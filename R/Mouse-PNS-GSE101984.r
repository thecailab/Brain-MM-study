library("Seurat")
library("patchwork")
library("Matrix")
library("plyr")
library("dplyr")

# GSE101984 data was downloaded from GEO database
setwd("C:\\Users\\xuanxuan\\Dropbox\\2021_Researches\\Xuanxuan\\Brain-scRNAseq\\GSE101984")
GSE101984<-read.csv("GSE101984_cells.csv",row.names = 2)
GSE101984<-GSE101984[,-1]

GSE101984.seurat<-CreateSeuratObject(GSE101984)

# First step: identify neurons
GSE101984.seurat <- NormalizeData(GSE101984.seurat, normalization.method = "LogNormalize", scale.factor = 1000000)
GSE101984.seurat <- FindVariableFeatures(GSE101984.seurat, selection.method = "vst")
GSE101984.seurat <- ScaleData(GSE101984.seurat)
GSE101984.seurat <- RunPCA(GSE101984.seurat, npcs=20, features = VariableFeatures(object = GSE101984.seurat))
GSE101984.seurat <- FindNeighbors(GSE101984.seurat, dims = 1:20)
GSE101984.seurat <- FindClusters(GSE101984.seurat, dim.use=1:20,resolution = 0.2)


# -------------------------------------------------------------------------------------# 
# firstly, the Trigeminal Neorons were identified using marker genes in the original study
arker<-c("Scn9a","Tubb3","Snap25","Trpv1", "Trpm8" , "Piezo2","Plp1", "Mbp" , "Epcam")
pdf("GSE101984level1.dotplot.pdf",width = 12,height = 8)
DotPlot(GSE101984.seurat, features = unique(marker), cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()
dev.off()


# according to dotplot, assign each cluster with cell type based on the expression of marker genes.
clusters<-as.numeric(GSE101984.seurat@active.ident)-1  # facotr 0 will be numeric 1, so that minus 1 to be consistant
table(clusters)
celltype<-clusters
celltype[which(clusters==0 | clusters==2 | clusters==4 | clusters==5 | clusters==6 
               | clusters==10 | clusters==12)]<-"Trigeminal Neorons"
celltype[which(clusters==1 | clusters==15 )]  <-"Epithelial"
celltype[which(clusters==8 | clusters==14)]   <-"Meylin"
celltype[which(clusters==3 | clusters==7 | clusters==9 | clusters==11 | clusters==13)]<-"Unknown"
Idents(GSE101984.seurat)<-celltype


# TSNE plot
GSE101984.seurat<- RunTSNE(GSE101984.seurat,dims=1:20, reduction = "pca")
pdf("GSE101984level1.tsne.pdf",width = 12,height = 8)
DimPlot(GSE101984.seurat, reduction = "tsne", label = T)
dev.off()


# To the interest of research, only keep cells of Trigeminal Neorons satisfying several criteria.
# criteria for inclusion: 500±7,500 genes, 0.2% mitochondrial transcripts; genes expressed in less
# than 6 neurons were also excluded leaving a dataset of 3580 neurons and more than 15,000 genes.
# 500±7,500 genes, 0.2% mitochondrial transcripts
GSE101984.puri<-GSE101984.seurat[,which(Idents(GSE101984.seurat)=="Trigeminal Neorons")]
GSE101984.puri[["percent.mt"]] <- PercentageFeatureSet(GSE101984.puri, pattern = "^MT-")
GSE101984.puri <- subset(GSE101984.puri, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 0.2)

# ------------------------------------------------------------------------- #
# 13 clusters in neurons
GSE101984.puri <- FindVariableFeatures(GSE101984.puri, selection.method = "vst")
GSE101984.puri <- ScaleData(GSE101984.puri)
GSE101984.puri <- RunPCA(GSE101984.puri, npcs=20, features = VariableFeatures(object = GSE101984.puri))
GSE101984.puri <- FindNeighbors(GSE101984.puri, dims = 1:20)
GSE101984.puri <- FindClusters(GSE101984.puri, dim.use=1:20,resolution = 0.3)

marker1<-read.csv("markers for mouse brain Fig5.csv")
marker2<-read.csv("markers for mouse brain FigS4.csv")
marker.gene<-c(marker1$markers,marker2$markers)

# m.g.l<-c()
# for (i in 1: length(marker.gene)){
#   m.g.l<-paste0(m.g.l,", ",marker.gene[i])
# }

pdf("GSE101984level1.dotplot.fig5.pdf",width = 30,height = 8)
DotPlot(GSE101984.puri, features = unique(marker.gene), cols = c("white","blue"),scale=F, dot.scale = 8) +
  RotatedAxis()
dev.off()

# assign cell type according to dotplot
clusters<-as.numeric(GSE101984.puri@active.ident)-1  # facotr 0 will be numeric 1, so that minus 1 to be consistant
table(clusters)
celltype<-clusters
celltype[which(clusters==0 | clusters==4)]<-"PEP"
celltype[which(clusters==1)]<-"Satellite"
celltype[which(clusters==2)]<-"cLTMR"
celltype[which(clusters==3 | clusters==7)]<-"NF"
celltype[which(clusters==5 | clusters==8)]<-"Endothelial"
celltype[which(clusters==9 | clusters==10)]<-"NP"
celltype[which(clusters==11)]<-"Schwann"
celltype[which(clusters==12)]<-"ENT.alike"
celltype[which(clusters==6)]<-"Unknown"
table(celltype)
Idents(GSE101984.puri)<-celltype
GSE101984.puri$cell.type<-celltype

GSE101984.puri<- RunTSNE(GSE101984.puri,dims=1:20, reduction = "pca")

pdf("GSE101984.filtered.tsne.pdf",width = 12,height = 8)
DimPlot(GSE101984.puri, reduction = "tsne", label = T)
dev.off()

# --------------------------------------------------------------------- #
# add meta ifo

GSE101984.puri@misc$DE$all <- FindAllMarkers(object = GSE101984.puri,
                                               only.pos = FALSE,
                                               min.pct = 0.1)
GSE101984.puri@misc$DE$all %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% as.data.frame() -> GSE101984.puri@misc$DE$top10
GSE101984.puri@misc$DE$all %>% group_by(cluster) %>% top_n(30, avg_log2FC) %>% as.data.frame() -> GSE101984.puri@misc$DE$top30
GSE101984.puri@misc$DE$all %>% group_by(cluster) %>% top_n(100, avg_log2FC) %>% as.data.frame() -> GSE101984.puri@misc$DE$top100
# -------------------------------------------------------------------------------#
# ------------add meta info -------------------- #
# data_loc   <- "C:/Users/YUXUANXUAN/Dropbox/2021_Researches/Xuanxuan/Brain-scRNAseq/Reference/"
# meta information template
meta.info.tmp <- read.csv(paste0(data_loc, "Tirosh-MM-2016-CD45P.csv"),
                          check.names = F,
                          row.names = 1,
                          na.strings = F,
                          stringsAsFactors = F)

meta1<-data.frame(Title="Diversity amongst trigeminal neurons revealed by high throughput single cell sequencing",
                  Authors="Minh Q. Nguyen, Youmei Wu, Lauren S. Bonilla, Lars J. von Buchholtz, Nicholas J.P. Ryba",
                  Publication= "PLOS ONE(2017)",
                  Summary="The trigeminal ganglion contains somatosensory neurons that detect a range of thermal, mechanical and chemical cues and innervate unique sensory compartments in the head and neck including the eyes, nose, mouth, meninges and vibrissae. We used single-cell sequencing and in situ hybridization to examine the cellular diversity of the trigeminal ganglion in mice, defining thirteen clusters of neurons. We show that clusters are well conserved in dorsal root ganglia suggesting they represent distinct functional classes of somatosensory neurons and not specialization associated with their sensory targets. Notably, functionally important genes (e.g. the mechanosensory channel Piezo2 and the capsaicin gated ion channel Trpv1) segregate into multiple clusters and often are expressed in subsets of cells within a cluster. Therefore, the 13 genetically-defined classes are likely to be physiologically heterogeneous rather than highly parallel (i.e., redundant) lines of sensory input. Our analysis harnesses the power of single-cell sequencing to provide a unique platform for in silico expression profiling that complements other approaches linking gene-expression with function and exposes unexpected diversity in the somatosensory system.",
                  Sample_Name="trigeminal neurons-PLOS ONE-2017",
                  Disease_Status="NA",
                  Tissue="NA",
                  Enrichment="NA",
                  Number_of_Cells="3943",
                  Platform="Illumina MiSeq",
                  Data_Source="GEO:GSE101984",
                  Data_POI="Minh Q Nguyen",
                  Analyst="NA",
                  Species="Mus musculus",
                  Reference_Genome="NA",
                  ENSEMBL="NA",
                  Project_ID="PRJNA396069",
                  Run_ID="NA",
                  Sample_ID="NA"
)
meta1<-as.data.frame(t(meta1))
rownames(meta1)<-rownames(meta.info.tmp)
colnames(meta1)<-colnames(meta.info.tmp)
meta.info <- as.list(subset(meta1, select = "Description", drop = T))
names(meta.info) <- rownames(meta1)
GSE101984.puri@misc$meta.info <- meta.info
# -------------------------------------------------------------------------------#
# ------------add demographic meta info -------------------- #
GSE101984.puri@meta.data$Age <- rep("3-6 wks",ncol(GSE101984.puri))
GSE101984.puri@misc$DataSegregation <- list("Age" = names(table(GSE101984.seurat@meta.data$Age)))