library("Seurat")
library("patchwork")
library("Matrix")
library("plyr")
library("dplyr")
library("R.utils") #Use gunzip function
library("GEOquery") # use getGEO function
# setwd("C:\\Users\\YUXUANXUAN\\Dropbox\\2021_Researches\\Cardiovascular Datasets\\GSE67835\\rawdata")

# untar files
untarfiles<-function(Dir){
  filenames<-list.files(Dir)
  full_pathname<-paste(Dir,"/",filenames,sep="")
  for (i in 1:length(full_pathname)) {
    gunzip(full_pathname[i],remove=TRUE)
  }
}
untarfiles("C:\\Users\\YUXUANXUAN\\Dropbox\\2021_Researches\\Cardiovascular Datasets\\GSE67825\\rawdata")

# read datasets, merge them by gene names
files<-list.files()
mat <- read.delim(file = files[1])
colnames(mat)<-c("gene",paste0("cell",1))
for (i in 2:length(files)){
  tmpmat <- read.delim(file = files[i])
  colnames(tmpmat)<-c("gene",paste0("cell",i))
  mat<-merge(mat,tmpmat,by="gene")
}
# modify rownames as gene names, then drop the gene names (first column) from the data matrix
rownames(mat)<-mat[,1]
mat<-mat[,-1]

# download demographic variables from GEO
TMP<-getGEO("GSE67835")
GSM.list<-c(TMP[["GSE67835-GPL15520_series_matrix.txt.gz"]]@phenoData@data[["geo_accession"]],
            TMP[["GSE67835-GPL18573_series_matrix.txt.gz"]]@phenoData@data[["geo_accession"]])
age<-c(TMP[["GSE67835-GPL15520_series_matrix.txt.gz"]]@phenoData@data$`age:ch1`,
       TMP[["GSE67835-GPL18573_series_matrix.txt.gz"]]@phenoData@data$`age:ch1`)
celltype<-c(TMP[["GSE67835-GPL15520_series_matrix.txt.gz"]]@phenoData@data$`cell type:ch1`,
            TMP[["GSE67835-GPL18573_series_matrix.txt.gz"]]@phenoData@data$`cell type:ch1`)
tissue<-c(TMP[["GSE67835-GPL15520_series_matrix.txt.gz"]]@phenoData@data$`tissue:ch1`,
          TMP[["GSE67835-GPL18573_series_matrix.txt.gz"]]@phenoData@data$`tissue:ch1`)

file.gsm<-substr(files,1,10)
orders<-match(file.gsm,GSM.list)
GSM.list<-GSM.list[orders]
age<-age[orders]
celltype<-celltype[orders]
tissue<-tissue[orders]

# create Seurate object
GSE67835<-CreateSeuratObject(mat)
# annotate the cells using the celltype identified in the original study
Idents(GSE67835)<-celltype

GSE67835 <- NormalizeData(GSE67835, normalization.method = "LogNormalize", scale.factor = 1000000)
GSE67835 <- FindVariableFeatures(GSE67835, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GSE67835)
GSE67835 <- ScaleData(GSE67835, features = all.genes)
GSE67835 <- RunPCA(GSE67835, features = VariableFeatures(object = GSE67835))
GSE67835 <- RunTSNE(GSE67835,dims = 1:10)


GSE67835@misc$DE$all <- FindAllMarkers(GSE67835, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GSE67835@misc$DE$all %>% group_by(cluster) %>% top_n(10, avg_logFC) %>% as.data.frame() -> GSE67835@misc$DE$top10
GSE67835@misc$DE$all %>% group_by(cluster) %>% top_n(30, avg_logFC) %>% as.data.frame() -> GSE67835@misc$DE$top30
GSE67835@misc$DE$all %>% group_by(cluster) %>% top_n(100, avg_logFC) %>% as.data.frame() -> GSE67835@misc$DE$top100

# setwd("C:\\Users\\YUXUANXUAN\\Dropbox\\2021_Researches\\Cardiovascular Datasets\\GSE67825")
# TSNE plot
pdf("GSE67835.tsne.pdf",width = 12,height = 8)
DimPlot(GSE67835,reduction = "tsne")
dev.off()

# -------------------------------------------------------------------------------#
# ------------add meta info -------------------- #
# data_loc   <- "C:/Users/YUXUANXUAN/Dropbox/2021_Researches/Xuanxuan/Brain-scRNAseq/Reference/"
# meta information template
meta.info.tmp <- read.csv(paste0(data_loc, "Tirosh-MM-2016-CD45P.csv"),
                          check.names = F,
                          row.names = 1,
                          na.strings = F,
                          stringsAsFactors = F)

meta1<-data.frame(Title="A survey of human brain transcriptome diversity at the single cell level",
                  Authors="Spyros Darmanis, Steven A. Sloanc, Ye Zhangc, Martin Engea, Christine Caneda, Lawrence M. Shuerd, Melanie G. Hayden Gephartd, Ben A. Barres, and Stephen R. Quake",
                  Publication= "PNAS(2015)",
                  Summary="We used single cell RNA sequencing on 466 cells to capture the cellular complexity of the adult and fetal human brain at a whole transcriptome level. Healthy adult temporal lobe tissue was obtained from epileptic patients during temporal lobectomy for medically refractory seizures. We were able to classify individual cells into all of the major neuronal, glial, and vascular cell types in the brain.",
                  Sample_Name="human brain-PNAS-2015",
                  Disease_Status="human brain",
                  Tissue="cortex",
                  Enrichment="NA",
                  Number_of_Cells="466",
                  Platform="Illumina MiSeq, Illumina NextSeq 500",
                  Data_Source="GEO:GSE67835",
                  Data_POI="Martin Enge[martin.enge@ki.se]",
                  Analyst="NA",
                  Species="Homo Sapiens",
                  Reference_Genome="hg19",
                  ENSEMBL="NA",
                  Project_ID="PRJNA281204",
                  Run_ID="NA",
                  Sample_ID="NA"
)
meta1<-as.data.frame(t(meta1))
rownames(meta1)<-rownames(meta.info.tmp)
colnames(meta1)<-colnames(meta.info.tmp)
meta.info <- as.list(subset(meta1, select = "Description", drop = T))
names(meta.info) <- rownames(meta1)
GSE67835@misc$meta.info <- meta.info
# -------------------------------------------------------------------------------#
# ------------add demographic meta info -------------------- #
GSE67835@meta.data$Age <- age
GSE67835@meta.data$Tissue <- tissue

GSE67835@misc$DataSegregation <- list("Age" = names(table(GSE67835@meta.data$Age)),
                                  "Tissue" = names(table(GSE67835@meta.data$Tissue)))

# saveRDS(GSE67835,"GSE67835.rds")

# ============================================================================== #
#45 marker analysis

genelist<-"marker45_M&D_group4"
marker.mat<-read.table(file.path("..",paste0(genelist,".txt")),header=TRUE,sep="\t")

SeuratObject.g<-SeuratObject[which(rownames(GSE67835) %in% marker.mat$HumanName),]

idents.ordered<-c("neurons","oligodendrocytes","microglia","astrocytes","endothelial")

SeuratObject.s<-SeuratObject.g[,which(Idents(SeuratObject) %in% idents.ordered)]
levels(SeuratObject.s) <- rev(idents.ordered)

#DotPlot
DotPlot(SeuratObject.s, features = unique(marker.mat$HumanName), cols = c("grey","blue"), dot.scale = 6) + RotatedAxis()

AveExp<-AverageExpression(SeuratObject.s, features = unique(marker.mat$HumanName), slot="data")
hm.dat<-t(AveExp$RNA)
hm.dat<-hm.dat[match(idents.ordered,rownames(hm.dat)),]

#Heatmap
pheatmap(hm.dat, scale = "column", cluster_rows = FALSE, cluster_cols = FALSE)

#DE analysis
DE <- FindAllMarkers(SeuratObject.s, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
