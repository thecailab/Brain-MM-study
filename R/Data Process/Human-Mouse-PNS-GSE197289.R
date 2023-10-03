library("Seurat")
library("dplyr")

# human and mouse data were downloaded from GEO database
# ===================================================================================== #
# human data
setwd("C:/Users/xuanxuan/Dropbox/2021_Researches/Xuanxuan/Brain-scRNAseq/GSE197289")
humanPNS.meta<-read.csv("GSE197289_snRNA-seq_human_barcode_meta.csv")
humanPNS<-readRDS("snRNA-seq_human_raw_counts.RDS")
humanPNS<-CreateSeuratObject(humanPNS)
humanPNS.meta<-humanPNS.meta[match(colnames(humanPNS),humanPNS.meta$V1),]

# use the cell type identified in the original study
Idents(humanPNS)<-humanPNS.meta$cellID

humanPNS <- NormalizeData(humanPNS, normalization.method = "LogNormalize", scale.factor = 1000000)
humanPNS <- FindVariableFeatures(humanPNS, selection.method = "vst")
humanPNS <- ScaleData(humanPNS)
humanPNS <- RunPCA(humanPNS, npcs=20, features = VariableFeatures(object = humanPNS))
humanPNS <- RunTSNE(humanPNS,dims=1:20, reduction = "pca")

pdf("human.TSNE_plot.pdf",width = 8,height = 7)
DimPlot(humanPNS, reduction = "tsne",label = T)
dev.off()

humanPNS@misc$DE$all <- FindAllMarkers(object = humanPNS,
                                        only.pos = FALSE,
                                        min.pct = 0.1)
humanPNS@misc$DE$all %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% as.data.frame() -> humanPNS@misc$DE$top10
humanPNS@misc$DE$all %>% group_by(cluster) %>% top_n(30, avg_log2FC) %>% as.data.frame() -> humanPNS@misc$DE$top30
humanPNS@misc$DE$all %>% group_by(cluster) %>% top_n(100, avg_log2FC) %>% as.data.frame() -> humanPNS@misc$DE$top100

# -------------------------------------------------------------------------------#
# ------------add meta info -------------------- #
# data_loc   <- "D:/Univeristy of South Carolina/2020FALL/ENHS793/Single cell datasets/20201204_SCANER_add_meta_info/20201203_add_meta_info/melanoma/"
data_loc   <- "C:/Users/xuanxuan/Desktop/SCANNER.add.meta.info/"
meta.info.tmp <- read.csv(paste0(data_loc, "Tirosh-MM-2016-CD45P.csv"),
                          check.names = F,
                          row.names = 1,
                          na.strings = F,
                          stringsAsFactors = F)

meta1<-data.frame(Title="Human and mouse trigeminal ganglia cell atlas implicates multiple cell types in migraine",
                  Authors="Lite Yang, Mengyi Xu, Shamsuddin A. Bhuiyan, Jia Li, Jun Zhao, Randall J. Cohrs, Justin T. Susterich, Sylvia Signorelli, Ursula Green, James R. Stone, Dan Levy, Jochen K. Lennerz, and William Renthal",
                  Publication= "Neuron (2022)",
                  Summary="The sensitization of trigeminal ganglion neurons contributes to primary headache disorders such as migraine, but the specific neuronal and non-neuronal trigeminal subtypes involved remain unclear. We thus developed a cell atlas in which human and mouse trigeminal ganglia are transcriptionally and epigenomically profiled at single-cell resolution. These data describe evolutionarily conserved and human-specific gene expression patterns within each trigeminal ganglion cell type, as well as the transcription factors and gene regulatory elements that contribute to cell-type-specific gene expression. We then leverage these data to identify trigeminal ganglion cell types that are implicated both by human genetic variation associated with migraine and two mouse models of headache. This trigeminal ganglion cell atlas improves our understanding of the cell types, genes, and epigenomic features involved in headache pathophysiology and establishes a rich resource of cell-type-specific molecular features to guide the development of more selective treatments for headache and facial pain.",
                  Sample_Name="Human",
                  Disease_Status=" trigeminal ganglia cell",
                  Tissue="trigeminal ganglia cell",
                  Enrichment="NA",
                  Number_of_Cells="38028",
                  Platform="NextSeq 550 (Homo sapiens)",
                  Data_Source="GEO:GSE197289",
                  Data_POI="yanglite1211@gmail.com",
                  Analyst="NA",
                  Species="Homo sapiens",
                  Reference_Genome="NA",
                  ENSEMBL="NA",
                  Project_ID="	PRJNA809683",
                  Run_ID="NA",
                  Sample_ID="1,2,3"
)
meta1<-as.data.frame(t(meta1))
rownames(meta1)<-rownames(meta.info.tmp)
colnames(meta1)<-colnames(meta.info.tmp)
# write.csv(meta1,"human.meta1.csv")
# meta1<-read.csv("human.meta1.csv",row.names = 1)
meta.info <- as.list(subset(meta1, select = "Description", drop = T))
names(meta.info) <- rownames(meta1)
humanPNS@misc$meta.info <- meta.info


# ------------add demographic meta info -------------------- #
humanPNS@meta.data$Tissue        <- humanPNS.meta$class
humanPNS@meta.data$SampleID      <- humanPNS.meta$donor
humanPNS@meta.data$Library       <- humanPNS.meta$library
humanPNS@meta.data$Ganglia       <- humanPNS.meta$ganglia
humanPNS@meta.data$sub.cell.type <- humanPNS.meta$subtype
humanPNS@misc$DataSegregation <- list("Tissue"   = names(table(humanPNS@meta.data$Tissue)),
                                       "SampleID" = names(table(humanPNS@meta.data$SampleID)),
                                       "Library"  = names(table(humanPNS@meta.data$Library)),
                                      "Ganglia"   = names(table(humanPNS@meta.data$Ganglia)),
                                       "Sub cell type" = names(table(humanPNS@meta.data$sub.cell.type))
)

saveRDS(humanPNS, "GSE197289_humanPNS.RDS")


# ===================================================================================== #
# Mouse data
setwd("C:/Users/xuanxuan/Dropbox/2021_Researches/Xuanxuan/Brain-scRNAseq/GSE197289")
mousePNS.meta<-read.csv("GSE197289_snRNA-seq_mouse_barcode_meta.csv")
mousePNS<-readRDS("snRNA-seq_mouse_raw_counts.RDS")

mousePNS<-CreateSeuratObject(mousePNS)
mousePNS.meta<-mousePNS.meta[match(colnames(mousePNS),mousePNS.meta$V1),]

# use the cell type identified in the original study
Idents(mousePNS)<-mousePNS.meta$cellID
mousePNS <- NormalizeData(mousePNS, normalization.method = "LogNormalize", scale.factor = 1000000)
mousePNS <- FindVariableFeatures(mousePNS, selection.method = "vst")
mousePNS <- ScaleData(mousePNS)
mousePNS <- RunPCA(mousePNS, npcs=20, features = VariableFeatures(object = mousePNS))
mousePNS <- RunTSNE(mousePNS,dims=1:20, reduction = "pca")


pdf("mouse.TSNE_plot.pdf",width = 8,height = 7)
DimPlot(mousePNS, reduction = "tsne",label = T)
dev.off()

mousePNS@misc$DE$all <- FindAllMarkers(object = mousePNS,
                                       only.pos = FALSE,
                                       min.pct = 0.1)
mousePNS@misc$DE$all %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% as.data.frame() -> mousePNS@misc$DE$top10
mousePNS@misc$DE$all %>% group_by(cluster) %>% top_n(30, avg_log2FC) %>% as.data.frame() -> mousePNS@misc$DE$top30
mousePNS@misc$DE$all %>% group_by(cluster) %>% top_n(100, avg_log2FC) %>% as.data.frame() -> mousePNS@misc$DE$top100

# -------------------------------------------------------------------------------#
# ------------add meta info -------------------- #
# data_loc   <- "D:/Univeristy of South Carolina/2020FALL/ENHS793/Single cell datasets/20201204_SCANER_add_meta_info/20201203_add_meta_info/melanoma/"
data_loc   <- "C:/Users/xuanxuan/Desktop/SCANNER.add.meta.info/"
meta.info.tmp <- read.csv(paste0(data_loc, "Tirosh-MM-2016-CD45P.csv"),
                          check.names = F,
                          row.names = 1,
                          na.strings = F,
                          stringsAsFactors = F)

meta1<-data.frame(Title="Human and mouse trigeminal ganglia cell atlas implicates multiple cell types in migraine",
                  Authors="Lite Yang, Mengyi Xu, Shamsuddin A. Bhuiyan, Jia Li, Jun Zhao, Randall J. Cohrs, Justin T. Susterich, Sylvia Signorelli, Ursula Green, James R. Stone, Dan Levy, Jochen K. Lennerz, and William Renthal",
                  Publication= "Neuron (2022)",
                  Summary="The sensitization of trigeminal ganglion neurons contributes to primary headache disorders such as migraine, but the specific neuronal and non-neuronal trigeminal subtypes involved remain unclear. We thus developed a cell atlas in which human and mouse trigeminal ganglia are transcriptionally and epigenomically profiled at single-cell resolution. These data describe evolutionarily conserved and human-specific gene expression patterns within each trigeminal ganglion cell type, as well as the transcription factors and gene regulatory elements that contribute to cell-type-specific gene expression. We then leverage these data to identify trigeminal ganglion cell types that are implicated both by human genetic variation associated with migraine and two mouse models of headache. This trigeminal ganglion cell atlas improves our understanding of the cell types, genes, and epigenomic features involved in headache pathophysiology and establishes a rich resource of cell-type-specific molecular features to guide the development of more selective treatments for headache and facial pain.",
                  Sample_Name="Mus musculus",
                  Disease_Status=" trigeminal ganglia cell",
                  Tissue="trigeminal ganglia cell",
                  Enrichment="NA",
                  Number_of_Cells="96933",
                  Platform="NextSeq 550 (Mus musculus)",
                  Data_Source="GEO:GSE197289",
                  Data_POI="yanglite1211@gmail.com",
                  Analyst="NA",
                  Species="Mus musculus",
                  Reference_Genome="NA",
                  ENSEMBL="NA",
                  Project_ID="PRJNA809683",
                  Run_ID="NA",
                  Sample_ID="CSD_1.5h_male_rep1, CSD_1.5h_male_rep1, CSD_1.5h_male_rep2, CSD_6h_male_rep1, CSD_6h_male_rep2, IS_1h_female_rep1, IS_1h_female_rep2, IS_1h_male_rep1, IS_1h_male_rep2, IS_1h_male_rep3, IS_24h_male_rep1, IS_24h_male_rep2, IS_6h_male_rep1, IS_6h_male_rep2, Naive_female_rep1, Naive_female_rep2, Naive_male_rep1, Naive_male_rep10, Naive_male_rep2, Naive_male_rep3, Naive_male_rep4, Naive_male_rep5, Naive_male_rep6, Naive_male_rep7, Naive_male_rep8, Naive_male_rep9, PBS_1h_female_rep1, PBS_1h_male_rep1"
)
meta1<-as.data.frame(t(meta1))
rownames(meta1)<-rownames(meta.info.tmp)
colnames(meta1)<-colnames(meta.info.tmp)
# write.csv(meta1,"mouse.meta1.csv")
# meta1<-read.csv("mouse.meta1.csv",row.names = 1)
meta.info               <- as.list(subset(meta1, select = "Description", drop = T))
names(meta.info)        <- rownames(meta1)
mousePNS@misc$meta.info <- meta.info


# ------------add demographic meta info -------------------- #
mousePNS@meta.data$Tissue        <- mousePNS.meta$class
mousePNS@meta.data$SampleID      <- mousePNS.meta$sample_name
mousePNS@meta.data$Library       <- mousePNS.meta$library
mousePNS@meta.data$Sex           <- mousePNS.meta$sex
mousePNS@meta.data$sub.cell.type <- mousePNS.meta$subtype
mousePNS@meta.data$Model         <- mousePNS.meta$model
mousePNS@meta.data$Realtime      <- mousePNS.meta$realtime
mousePNS@misc$DataSegregation <- list("Tissue"   = names(table(mousePNS@meta.data$Tissue)),
                                      "SampleID" = names(table(mousePNS@meta.data$SampleID)),
                                      "Library"  = names(table(mousePNS@meta.data$Library)),
                                      "Sub cell type" = names(table(mousePNS@meta.data$sub.cell.type)),
                                      "Sex"      = names(table(mousePNS@meta.data$Sex)),
                                      "Model"    = names(table(mousePNS@meta.data$Model)),
                                      "Real time"= names(table(mousePNS@meta.data$Realtime))
)

saveRDS(mousePNS, "GSE197289_mousePNS.RDS")


# ============================================================================== #
#45 marker analysis

data.l<-list("human"=humanPNS,"mouse"=mousePNS)

genelist<-"marker45_M&D_group4"
marker.mat<-read.table(file.path("..",paste0(genelist,".txt")),header=TRUE,sep="\t")

for(species.i in names(data.l)){
  if(species.i=="human"){
		marker.vec<-marker.mat$HumanName
	}
	if(species.i=="mouse"){
		marker.vec<-marker.mat$MouseName
	}
  
  SeuratObject<-data.l[[species.i]]
  SeuratObject.g<-SeuratObject[which(rownames(SeuratObject) %in% marker.vec),]

  SeuratObject.g<-SeuratObject.g[,which(SeuratObject.g$Model %in% c("Control","Naive"))]

  idents.ordered<-c("NF","PEP","NP","Schwann","Satglia","cLTMR")

  SeuratObject.s<-SeuratObject.g[,which(Idents(SeuratObject.g) %in% idents.ordered)]
  levels(SeuratObject.s) <- rev(idents.ordered)

  #DotPlot
  DotPlot(SeuratObject.s, features = unique(marker.vec), cols = c("grey","blue"), dot.scale = 6) + RotatedAxis()

  AveExp<-AverageExpression(SeuratObject.s, features = unique(marker.vec), slot="data")
  hm.dat<-t(AveExp$RNA)
  hm.dat<-hm.dat[match(idents.ordered,rownames(hm.dat)),]

  #Heatmap
  pheatmap(hm.dat, scale = "column", cluster_rows = FALSE, cluster_cols = FALSE)

  #DE analysis
  DE <- FindAllMarkers(SeuratObject.s, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
}

