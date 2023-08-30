library("LoomExperiment") # read data
library("Seurat")
library("patchwork")
library("Matrix")
library("plyr")
library("dplyr")
library("SeuratObject")

setwd("/home/xuanxuan/20210609_mouse brain scRNAseq/mouse_brian.org")
mouse.brain <- import("l5_all.loom", type="SingleCellLoomExperiment")

meta.data<-as.data.frame(mouse.brain@colData@listData)
ind<-which(meta.data$Taxonomy_group=="Astrocytes" |
             meta.data$Taxonomy_group=="Oligodendrocytes" |
             meta.data$Taxonomy_group=="Microglia" |
             meta.data$Taxonomy_group=="Di- and mesencephalon excitatory neurons" |
             meta.data$Taxonomy_group=="Di- and mesencephalon inhibitory neurons" |
             meta.data$Taxonomy_group=="Spinal cord excitatory neurons" |
             meta.data$Taxonomy_group=="Spinal cord inhibitory neurons" |
             meta.data$Taxonomy_group=="Telencephalon projecting excitatory neurons" |
             meta.data$Taxonomy_group=="Telencephalon projecting inhibitory neurons" |
             # meta.data$Taxonomy_group=="Peripheral sensory peptidergic neurons" |
             # meta.data$Taxonomy_group=="Peripheral sensory non-peptidergic neurons" |
             # meta.data$Taxonomy_group=="Peripheral sensory neurofilament neurons" |
             meta.data$Taxonomy_group=="Sympathetic cholinergic neurons" |
             meta.data$Taxonomy_group=="Sympathetic noradrenergic neurons" |
             meta.data$Taxonomy_group=="Schwann cells" |
             meta.data$Taxonomy_group=="Satellite glia" |
             meta.data$Taxonomy_group=="Vascular and leptomeningeal cells" |
             meta.data$Taxonomy_group=="Vascular endothelial cells" |
             meta.data$Taxonomy_group=="Vascular smooth muscle cells" |
             meta.data$Description=="Neurofilament (NF1), DRG" |
             meta.data$Description=="Neurofilament (NF2/3), DRG" |
             meta.data$Description=="Neurofilament (NF4/5), DRG" |
             meta.data$Description=="Non-peptidergic (NP1.1), DRG" |
             meta.data$Description=="Non-peptidergic (NP1.2), DRG" |
             meta.data$Description=="Non-peptidergic (NP2.1), DRG" |
             meta.data$Description=="Non-peptidergic (NP2.2), DRG" |
             meta.data$Description=="Non-peptidergic (NP3), DRG" |
             meta.data$Description=="Non-peptidergic (TH), DRG" |
             meta.data$Description=="Peptidergic (PEP1..4), DRG" |
             meta.data$Description=="Peptidergic (PEP1.1), DRG" |
             meta.data$Description=="Peptidergic (PEP1.2), DRG" |
             meta.data$Description=="Peptidergic (PEP1.3), DRG" |
             meta.data$Description=="Peptidergicv (PEP2), DRG" |
             meta.data$Description=="Peptidergic (TrpM8), DRG")

meta.data<-meta.data[ind,]
mat1<-mouse.brain@assays@data@listData[["matrix"]][,ind]

colnames(mat1)<-meta.data$CellID
rownames(mat1)<-mouse.brain@rowRanges@elementMetadata@listData[["Gene"]]

# assign cell type
meta.data$cell.type<-meta.data$Age
meta.data$cell.type[which(meta.data$Taxonomy_group=="Astrocytes")]<-"Astrocytes"
meta.data$cell.type[which(meta.data$Taxonomy_group=="Oligodendrocytes")]<-"Oligodendrocytes"
meta.data$cell.type[which(meta.data$Taxonomy_group=="Microglia")]<-"Microglia"
meta.data$cell.type[which(meta.data$Taxonomy_group=="Di- and mesencephalon excitatory neurons" |
                            meta.data$Taxonomy_group=="Spinal cord excitatory neurons" |
                            meta.data$Taxonomy_group=="Telencephalon projecting excitatory neurons")]<-"Excitatory"
meta.data$cell.type[which(meta.data$Taxonomy_group=="Di- and mesencephalon inhibitory neurons" |
                            meta.data$Taxonomy_group=="Spinal cord inhibitory neurons" |
                            meta.data$Taxonomy_group=="Telencephalon projecting inhibitory neurons")]<-"Inhibitory"
meta.data$cell.type[which(meta.data$Taxonomy_group=="Schwann cells")]<-"Schwann"
meta.data$cell.type[which(meta.data$Taxonomy_group=="Satellite glia")]<-"Satellite glia"
meta.data$cell.type[which(meta.data$Taxonomy_group=="Vascular and leptomeningeal cells")]<-"Vasc/leptomening"
meta.data$cell.type[which(meta.data$Taxonomy_group=="Vascular endothelial cells")]<-"Vasc/endoth"
meta.data$cell.type[which(meta.data$Taxonomy_group=="Vascular smooth muscle cells")]<-"Vasc/smooth mus"
meta.data$cell.type[which(meta.data$Taxonomy_group=="Sympathetic cholinergic neurons")]<-"Symp chol."
meta.data$cell.type[which(meta.data$Taxonomy_group=="Sympathetic noradrenergic neurons")]<-"Symp norad."
meta.data$cell.type[which(meta.data$Description=="Neurofilament (NF1), DRG" |
                            meta.data$Description=="Neurofilament (NF2/3), DRG" |
                            meta.data$Description=="Neurofilament (NF4/5), DRG")]<-"NF"
meta.data$cell.type[which(meta.data$Description=="Non-peptidergic (NP1.1), DRG" |
                            meta.data$Description=="Non-peptidergic (NP1.2), DRG" |
                            meta.data$Description=="Non-peptidergic (NP2.1), DRG" |
                            meta.data$Description=="Non-peptidergic (NP2.2), DRG" |
                            meta.data$Description=="Non-peptidergic (NP3), DRG" |
                            meta.data$Description=="Non-peptidergic (TH), DRG")]<-"NP"
meta.data$cell.type[which(meta.data$Description=="Peptidergic (TrpM8), DRG")]<-"cLTMR"
meta.data$cell.type[which(meta.data$Description=="Peptidergic (PEP1..4), DRG" |
                            meta.data$Description=="Peptidergic (PEP1.1), DRG" |
                            meta.data$Description=="Peptidergic (PEP1.2), DRG" |
                            meta.data$Description=="Peptidergic (PEP1.3), DRG" |
                            meta.data$Description=="Peptidergicv (PEP2), DRG")]<-"PEP"

mouse.brain.seurat<-CreateSeuratObject(mat1)
mouse.brain.seurat@meta.data<-cbind(mouse.brain.seurat@meta.data,meta.data)
Idents(mouse.brain.seurat)<-meta.data$cell.type

# add meta infos
mouse.brain.seurat <- NormalizeData(mouse.brain.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
mouse.brain.seurat <- FindVariableFeatures(mouse.brain.seurat, selection.method = "vst")
mouse.brain.seurat <- ScaleData(mouse.brain.seurat)
mouse.brain.seurat <- RunPCA(mouse.brain.seurat, npcs=20, features = VariableFeatures(object = mouse.brain.seurat))
mouse.brain.seurat<- RunTSNE(mouse.brain.seurat,dims=1:20, reduction = "pca")

pdf("mouse.brain.tsne.pdf",width = 12,height = 8)
DimPlot(mouse.brain.seurat, reduction = "tsne", label = T)
dev.off()

mouse.brain.seurat@misc$DE$all <- FindAllMarkers(object = mouse.brain.seurat,
                                                 only.pos = FALSE,
                                                 min.pct = 0.1)

mouse.brain.seurat@misc$DE$all %>% group_by(cluster) %>% top_n(10, avg_logFC) %>% as.data.frame() -> mouse.brain.seurat@misc$DE$top10
mouse.brain.seurat@misc$DE$all %>% group_by(cluster) %>% top_n(30, avg_logFC) %>% as.data.frame() -> mouse.brain.seurat@misc$DE$top30
mouse.brain.seurat@misc$DE$all %>% group_by(cluster) %>% top_n(100, avg_logFC) %>% as.data.frame() -> mouse.brain.seurat@misc$DE$top100

dim(mouse.brain.seurat)
names(mouse.brain.seurat@misc$DE$all)
# -------------------------------------------------------------------------------#
# ------------add meta info -------------------- #
# setwd("C:/Users/xuanxuan/Dropbox/2021_Researches/Xuanxuan/Brain-scRNAseq/mouse_Brain.org")
meta.info.tmp <- read.csv("Tirosh-MM-2016-CD45P.csv",check.names = F,row.names = 1,na.strings = F,stringsAsFactors = F)

meta1<-data.frame(Title="Molecular Architecture of the Mouse Nervous System",
                  Authors="Amit Zeisel, Hannah Hochgerner, Peter Lönnerberg, Anna Johnsson, Fatima Memic, Job van der Zwan, Martin Häring, Emelie Braun, Lars E Borm, Gioele La Manno, Simone Codeluppi, Alessandro Furlan, Kawai Lee, Nathan Skene, Kenneth D Harris, Jens Hjerling-Leffler, Ernest Arenas, Patrik Ernfors, Ulrika Marklund, Sten Linnarsson",
                  Publication= "Cell(2018)",
                  Summary="The mammalian nervous system executes complex behaviors controlled by specialized, precisely positioned, and interacting cell types. Here, we used RNA sequencing of half a million single cells to create a detailed census of cell types in the mouse nervous system. We mapped cell types spatially and derived a hierarchical, data-driven taxonomy. Neurons were the most diverse and were grouped by developmental anatomical units and by the expression of neurotransmitters and neuropeptides. Neuronal diversity was driven by genes encoding cell identity, synaptic connectivity, neurotransmission, and membrane conductance. We discovered seven distinct, regionally restricted astrocyte types that obeyed developmental boundaries and correlated with the spatial distribution of key glutamate and glycine neurotransmitters. In contrast, oligodendrocytes showed a loss of regional identity followed by a secondary diversification. The resource presented here lays a solid foundation for understanding the molecular architecture of the mammalian nervous system and enables genetic manipulation of specific cell types.",
                  Sample_Name="mouse brain-Cell-2018",
                  Disease_Status="NA",
                  Tissue="NA",
                  Enrichment="NA",
                  Number_of_Cells="102079",
                  Platform="10X Genomics Chromium Single Cell Kit Version 1/Version 2",
                  Data_Source="http://mousebrain.org/downloads.html",
                  Data_POI="NA",
                  Analyst="NA",
                  Species="mice",
                  Reference_Genome="NA",
                  ENSEMBL="NA",
                  Project_ID="NA",
                  Run_ID="NA",
                  Sample_ID="NA"
)
meta1<-as.data.frame(t(meta1))
rownames(meta1)<-rownames(meta.info.tmp)
colnames(meta1)<-colnames(meta.info.tmp)
# write.csv(meta1,"meta1.csv")
# meta1<-read.csv("meta1.csv")
meta.info <- as.list(subset(meta1, select = "Description", drop = T))
names(meta.info) <- rownames(meta1)
mouse.brain.seurat@misc$meta.info <- meta.info
# -------------------------------------------------------------------------------#
# ------------add demographic meta info -------------------- #
mouse.brain.seurat@misc$DataSegregation <- list("Age" = names(table(mouse.brain.seurat@meta.data$Age)))
mouse.brain.seurat@misc$DataSegregation <- list("Taxonomy_group" = names(table(mouse.brain.seurat@meta.data$Taxonomy_group)))
mouse.brain.seurat@misc$DataSegregation <- list("Description" = names(table(mouse.brain.seurat@meta.data$Description)))
mouse.brain.seurat@misc$DataSegregation <- list("Class" = names(table(mouse.brain.seurat@meta.data$Class)))
mouse.brain.seurat@misc$DataSegregation <- list("Subclass" = names(table(mouse.brain.seurat@meta.data$Subclass)))
mouse.brain.seurat@misc$DataSegregation <- list("Tissue" = names(table(mouse.brain.seurat@meta.data$Tissue)))
mouse.brain.seurat@misc$DataSegregation <- list("Clusters" = names(table(mouse.brain.seurat@meta.data$Clusters)))
mouse.brain.seurat@misc$DataSegregation <- list("TaxonomyRank1" = names(table(mouse.brain.seurat@meta.data$TaxonomyRank1)))
mouse.brain.seurat@misc$DataSegregation <- list("TaxonomyRank2" = names(table(mouse.brain.seurat@meta.data$TaxonomyRank2)))
mouse.brain.seurat@misc$DataSegregation <- list("TaxonomyRank3" = names(table(mouse.brain.seurat@meta.data$TaxonomyRank3)))
mouse.brain.seurat@misc$DataSegregation <- list("TaxonomyRank4" = names(table(mouse.brain.seurat@meta.data$TaxonomyRank4)))


saveRDS(mouse.brain.seurat, "mouse.brain.seurat.RDS")
