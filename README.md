# Brain-MM-study

## Datasets

Four CNS and four PNS datasets were analyzed in this study, including two human CNS datasets (GSE67835, 267 cells; GSE97930, 33,722 cells),  two mouse CNS datasets (the Mouse Brain Atlas, 95,753 cells; GSE60361, 1,667 cells), three mouse PNS datasets (the Mouse Brain Atlas, 2,311 cells; GSE101984, 3,688 cells; GSE197289, 67,991), and a human PNS dataset (GSE197289, 17,223 cells). The mouse PNS dataset, GSE197289, are from  mouse of normal models (PBS) and mouse of two headache models (CSD, IS), while all other datasets are from normal tissue. In CNS and PNS datasets, we studied cell types that considered to be involved in migraine pathophysiology, including CNS cell types (neurons, oligodendrocytes, microglia, astrocytes); PNS cell types (peptidergic nociceptors (PEP), non-peptidergic nociceptors (NP), large diameter neurofilament-positive mechanoreceptors (NF), C low threshold mechanoreceptors (cLTMR), Satellite cells, Schwann cells); and neurovascular endothelial cells (vasc/endoth).

## Data processing and cell type identification

The sequencing counts of single cell RNA-seq datasets were processed using Seurat. Cell types in the mouse Brain Atlas datasets, the mouse CNS dataset(GSE60361), the human CNS datasets (GSE67835, GSE97930), the PNS datasets of human and mouse cells (GSE197289) were identified by their original studies, respectively. For all data, log cpm (counts per million) were calculated by “NormalizeData” function and then centered and scaled by its standard deviation.  For the mouse dataset of trigeminal ganglion neurons (GSE101984), the principal component analysis (PCA) was performed on the top 2000 variable genes identified by “FindVariableFeatures” function. Cells were assigned into 13 clusters using “FindClusters” function with the top 20 PCs as the inputs and resolution was set as 0.3. The average expressions and expressed percentages of marker genes used by the Mouse Brain Atlas study were used to identify NF, satellite, endothelial, cLTMR, schwann, PEP, NP cell types. 

## Differential expression analysis

The Wilcox rank sum test implemented in “Seurat” was performed to compare gene expression in a certain cell type with that in all other cell types. Further, Bonferroni method was applied to adjust 𝑝 values to control false positives due to multiple comparison. A significant difference was detected when three Seurat default criteria were satisfied, including adjust 𝑝 value was less than 0.05, log fold change was larger than 0.25, and expression was detected in more than 10% of cells.
