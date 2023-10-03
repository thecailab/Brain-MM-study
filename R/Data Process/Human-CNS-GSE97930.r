setwd("C:/Users/YUXUANXUAN/Dropbox/2021_Researches/Xuanxuan/Brain-scRNAseq/GSE97930_snDrop_seq")
setwd("C:/Users/xuanxuan/Dropbox/2021_Researches/Xuanxuan/Brain-scRNAseq/GSE97930_snDrop_seq")

celltypeinfo<-read.csv("celltypeinfo.csv")
geographyinfo<-read.csv("geographyinfo.csv")
geographyinfo$Library=tolower(geographyinfo$Library)
celltypeinfo<-merge(celltypeinfo,geographyinfo,by="Library")
# table(geographyinfo$Library.Barcode)
# table(geographyinfo$Experiment)
# table(geographyinfo$Library)
# table(celltypeinfo$Library)

excitory<-c("Ex1","Ex2","Ex3a","Ex3b","Ex3c","Ex3d","Ex3e","Ex4","Ex5a","Ex5b","Ex6a","Ex6b","Ex8","Gran")
inibitory<-c("In1a","In1b","In1c","In2","In3","In4a","In4b","In6a","In6b","In7","In8")
celltypeinfo$celltype=celltypeinfo$Identity
celltypeinfo$celltype[which(celltypeinfo$Identity %in% c(excitory,inibitory))]<-"Neurons"
celltypeinfo$celltype[which(celltypeinfo$Identity == "Oli")]<-"Oligodendrocyte"
celltypeinfo$celltype[which(celltypeinfo$Identity == "Mic")]<-"Microglia"
celltypeinfo$celltype[which(celltypeinfo$Identity %in% c("Ast","Ast_Cer"))]<-"Astrocyte"
celltypeinfo$celltype[which(celltypeinfo$Identity == "End")]<-"Endothelial"
celltypeinfo<-celltypeinfo[which(celltypeinfo$celltype %in% c("Neurons","Oligodendrocyte","Microglia","Astrocyte","Endothelial")),]
celltypeinfo$cellID<-paste0(celltypeinfo$Identity,"_",celltypeinfo$Sample.Names..Library_Barcode.)

table(celltypeinfo$celltype)
# match(names(table(celltypeinfo$Library)),tolower(geographyinfo$Library))


# ============================================================================== #
# FrontalCortex.seurat
DATA1<-read.delim("GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")
keep.list<-celltypeinfo$cellID[na.omit(match(colnames(DATA1),celltypeinfo$cellID))]
# na.omit(match(colnames(DATA1),keep.list))
celltypeinfo1<-celltypeinfo[match(keep.list,celltypeinfo$cellID),]
DATA1<-DATA1[,match(keep.list,colnames(DATA1))]

# create Seurate object
GSE97930_FrontalCortex<-CreateSeuratObject(DATA1)
# annotate the cells using the celltype identified in the original study
Idents(GSE97930_FrontalCortex)<-celltypeinfo1$celltype

GSE97930_FrontalCortex <- NormalizeData(GSE97930_FrontalCortex, normalization.method = "LogNormalize", scale.factor = 1000000)
GSE97930_FrontalCortex <- FindVariableFeatures(GSE97930_FrontalCortex, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GSE97930_FrontalCortex)
GSE97930_FrontalCortex <- ScaleData(GSE97930_FrontalCortex, features = all.genes)
GSE97930_FrontalCortex <- RunPCA(GSE97930_FrontalCortex, features = VariableFeatures(object = GSE97930_FrontalCortex))
GSE97930_FrontalCortex <- RunTSNE(GSE97930_FrontalCortex,dims = 1:10)


GSE97930_FrontalCortex@misc$DE$all <- FindAllMarkers(GSE97930_FrontalCortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GSE97930_FrontalCortex@misc$DE$all %>% group_by(cluster) %>% top_n(10, avg_logFC) %>% as.data.frame() -> GSE97930_FrontalCortex@misc$DE$top10
GSE97930_FrontalCortex@misc$DE$all %>% group_by(cluster) %>% top_n(30, avg_logFC) %>% as.data.frame() -> GSE97930_FrontalCortex@misc$DE$top30
GSE97930_FrontalCortex@misc$DE$all %>% group_by(cluster) %>% top_n(100, avg_logFC) %>% as.data.frame() -> GSE97930_FrontalCortex@misc$DE$top100

# setwd("C:\\Users\\YUXUANXUAN\\Dropbox\\2021_Researches\\Cardiovascular Datasets\\GSE67825")
# TSNE plot
pdf("GSE97930_FrontalCortex.tsne.pdf",width = 12,height = 8)
DimPlot(GSE97930_FrontalCortex,reduction = "tsne")
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

meta1<-data.frame(Title="Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain",
                  Authors="Blue B Lake, Song Chen, Brandon C Sos, Jean Fan, Gwendolyn E Kaeser, Yun C Yung, Thu E Duong, Derek Gao, Jerold Chun, Peter V Kharchenko, Kun Zhang",
                  Publication= "Nature Biotechnology(2017)",
                  Summary="Detailed characterization of the cell types in the human brain requires scalable experimental approaches to examine multiple aspects of the molecular state of individual cells, as well as computational integration of the data to produce unified cell-state annotations. Here we report improved high-throughput methods for single-nucleus droplet-based sequencing (snDrop-seq) and single-cell transposome hypersensitive site sequencing (scTHS-seq). We used each method to acquire nuclear transcriptomic and DNA accessibility maps for >60,000 single cells from human adult visual cortex, frontal cortex, and cerebellum. Integration of these data revealed regulatory elements and transcription factors that underlie cell-type distinctions, providing a basis for the study of complex processes in the brain, such as genetic programs that coordinate adult remyelination. We also mapped disease-associated risk variants to specific cellular populations, which provided insights into normal and pathogenic cellular processes in the human brain. This integrative multi-omics approach permits more detailed single-cell interrogation of complex organs and tissues.",
                  Sample_Name="human Frontal Cortex",
                  Disease_Status="human Frontal Cortex",
                  Tissue="Multiple",
                  Enrichment="NA",
                  Number_of_Cells="10319",
                  Platform="Illumina MiSeq (Homo sapiens), Illumina HiSeq 2500 (Homo sapiens)",
                  Data_Source="GEO:GSE97930",
                  Data_POI="Kun Zhang (UCSD)",
                  Analyst="NA",
                  Species="Homo Sapiens",
                  Reference_Genome="NA",
                  ENSEMBL="NA",
                  Project_ID="PRJNA383372",
                  Run_ID="NA",
                  Sample_ID="0006-YO,0007-OX,0008-CR,5342"
)
meta1<-as.data.frame(t(meta1))
rownames(meta1)<-rownames(meta.info.tmp)
colnames(meta1)<-colnames(meta.info.tmp)
meta.info <- as.list(subset(meta1, select = "Description", drop = T))
names(meta.info) <- rownames(meta1)
GSE97930_FrontalCortex@misc$meta.info <- meta.info
# -------------------------------------------------------------------------------#
# ------------add demographic meta info -------------------- #
GSE97930_FrontalCortex@meta.data$Age          <- celltypeinfo1$AGE
GSE97930_FrontalCortex@meta.data$SEX          <- celltypeinfo1$SEX
GSE97930_FrontalCortex@meta.data$RACE         <- celltypeinfo1$RACE
GSE97930_FrontalCortex@meta.data$Identity     <- celltypeinfo1$Identity
GSE97930_FrontalCortex@meta.data$Patient      <- celltypeinfo1$Patient.UMB..x
GSE97930_FrontalCortex@meta.data$Library      <- celltypeinfo1$Library
GSE97930_FrontalCortex@meta.data$Brain.Region <- celltypeinfo1$Brain.Region.x

GSE97930_FrontalCortex@misc$DataSegregation <- list("AGE" = names(table(GSE97930_FrontalCortex@meta.data$AGE)),
                                                    "SEX" = names(table(GSE97930_FrontalCortex@meta.data$SEX)),
                                                    "RACE" = names(table(GSE97930_FrontalCortex@meta.data$RACE)),
                                                    "Identity" = names(table(GSE97930_FrontalCortex@meta.data$Identity)),
                                                    "Patient" = names(table(GSE97930_FrontalCortex@meta.data$Patient)),
                                                    "Library" = names(table(GSE97930_FrontalCortex@meta.data$Library)),
                                                    "Brain.Region" = names(table(GSE97930_FrontalCortex@meta.data$Brain.Region))
                                      )

save(GSE97930_FrontalCortex,"GSE97930_FrontalCortex.RDS")

# ============================================================================== #
#45 marker analysis

genelist<-"marker45_M&D_group4"
marker.mat<-read.table(file.path("..",paste0(genelist,".txt")),header=TRUE,sep="\t")

SeuratObject.g<-GSE97930_FrontalCortex[which(rownames(GSE97930_FrontalCortex) %in% marker.mat$HumanName),]

SeuratObject.g <- RenameIdents(object = SeuratObject.g, `Excitatory Neurons` = "Neurons", `Inhibitory Neurons` = "Neurons")

idents.ordered<-c("Neurons","Oligodendrocyte","Microglia","Astrocyte","Endothelial")

SeuratObject.s<-SeuratObject.g[,which(Idents(SeuratObject.g) %in% idents.ordered)]
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
