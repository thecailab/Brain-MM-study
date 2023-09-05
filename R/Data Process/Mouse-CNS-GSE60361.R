library("Seurat")
library("patchwork")
library("Matrix")
library("plyr")
library("dplyr")

setwd("C:/Users/xuanxuan/Dropbox/2021_Researches/Xuanxuan/Brain-scRNAseq/GSE60361")
GSE603612<-readRDS("GSE60361.RDS")
GSE603612@misc$meta.info
expression_mRNA<-read.delim("expression_mRNA_17-Aug-2014.txt")

meta.data<-as.data.frame(expression_mRNA[1:9,])
meta.data[10,]<-colnames(meta.data)
meta.data<-as.data.frame(t(meta.data))
colnames(meta.data)<-meta.data[2,]
meta.data<-meta.data[-c(1,2),]
table(meta.data$level1class)

data.mat<-expression_mRNA[11:nrow(expression_mRNA),]
rownames(data.mat)<-data.mat[,1]
data.mat<-data.mat[,-c(1,2)]
colnames(data.mat)<-meta.data$cell_id

GSE60361<-CreateSeuratObject(data.mat)
Idents(GSE60361)<-meta.data$level1class

GSE60361 <- NormalizeData(GSE60361, normalization.method = "LogNormalize", scale.factor = 1000000)
GSE60361 <- FindVariableFeatures(GSE60361, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GSE60361)
GSE60361 <- ScaleData(GSE60361, features = all.genes)
GSE60361 <- RunPCA(GSE60361, features = VariableFeatures(object = GSE60361))
GSE60361 <- RunTSNE(GSE60361,dims = 1:10)


GSE60361@misc$DE$all <- FindAllMarkers(GSE60361, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GSE60361@misc$DE$all %>% group_by(cluster) %>% top_n(10, avg_logFC) %>% as.data.frame() -> GSE60361@misc$DE$top10
GSE60361@misc$DE$all %>% group_by(cluster) %>% top_n(30, avg_logFC) %>% as.data.frame() -> GSE60361@misc$DE$top30
GSE60361@misc$DE$all %>% group_by(cluster) %>% top_n(100, avg_logFC) %>% as.data.frame() -> GSE60361@misc$DE$top100

# TSNE plot
pdf("GSE60361.tsne.pdf",width = 12,height = 8)
DimPlot(GSE60361,reduction = "tsne")
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

meta1<-data.frame(Title="Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq",
                  Authors="Amit Zeisel, Ana B. Muñoz-Manchado, Simone Codeluppi, Peter Lönnerberg, Gioele La Manno, Anna Juréus, Sueli Marques, Hermany Munguba, Liqun He, Christer Betsholtz, Charlotte Rolny, Gonçalo Castelo-Branco, Jens Hjerling-Leffler, Sten Linnarsson",
                  Publication= "Science(2015)",
                  Summary="The mammalian cerebral cortex supports cognitive functions such as sensorimotor integration, memory, and social behaviors. Normal brain function relies on a diverse set of differentiated cell types, including neurons, glia, and vasculature. Here,we have used large-scale single-cell RNA sequencing (RNA-seq) to classify cells in the mouse somatosensory cortex and hippocampal CA1 region.We found 47 molecularly distinct subclasses, comprising all known major cell types in the cortex.We identified numerous marker genes, which allowed alignment with known cell types, morphology, and location.We found a layer I interneuron expressing Pax6 and a distinct postmitotic oligodendrocyte subclass marked by Itpr2. Across the diversity of cortical cell types, transcription factors formed a complex, layered regulatory code, suggesting a mechanism for the maintenance of adult cell type identity.",
                  Sample_Name="cortex and hippocampus_Science(2015)",
                  Disease_Status="NA",
                  Tissue="cortex and hippocampus",
                  Enrichment="NA",
                  Number_of_Cells="3005",
                  Platform="Illumina MiSeq, Illumina NextSeq 500",
                  Data_Source="GEO:GSE60361;http://linnarssonlab.org/cortex/",
                  Data_POI="Amit Zeisel",
                  Analyst="NA",
                  Species="mouse",
                  Reference_Genome="NA",
                  ENSEMBL="NA",
                  Project_ID="NA",
                  Run_ID="NA",
                  Sample_ID="NA"
)
meta1<-as.data.frame(t(meta1))
rownames(meta1)<-rownames(meta.info.tmp)
colnames(meta1)<-colnames(meta.info.tmp)
meta.info <- as.list(subset(meta1, select = "Description", drop = T))
names(meta.info) <- rownames(meta1)
GSE60361@misc$meta.info <- meta.info
# -------------------------------------------------------------------------------#
# ------------add demographic meta info -------------------- #
GSE60361@meta.data$Age <- meta.data$age
GSE60361@meta.data$Tissue <- meta.data$tissue

GSE60361@misc$DataSegregation <- list("Age" = names(table(GSE60361@meta.data$Age)),
                                      "Tissue" = names(table(GSE60361@meta.data$Tissue)))

# saveRDS(GSE60361,"GSE60361.rds")
