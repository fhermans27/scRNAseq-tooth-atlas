# scRNAseq-tooth-atlas

This repository contains the R and Python code used perform the single-cell RNA-sequencing (scRNAseq) analysis used in :

*"Establishment of inclusive single-cell transcriptome atlases from mouse and human tooth as powerful resource for dental research", Hermans, F, Bueds , C, Hemeryck, L, Lambrichts, I, Bronckaers, A and Vankelecom, H. (2022) Front. Cell Dev. Biol. doi: 10.3389/fcell.2022.1021459*

Which can be found at: https://www.frontiersin.org/articles/10.3389/fcell.2022.1021459/full 

# Table of contents
- [Mouse Tooth Atlas](#mouse-tooth-atlas)
  - [Setup](#mta-setup)
  - [Quality Control](#mta-quality-control)
  - [Initial Integration on three groups](#initial-integration-on-three-groups)
  - [Removal of background or ambient RNA using SoupX](#removal-of-background-or-ambient-rna-using-soupx)
  - [rPCA integration of all data](#rpca-integration-of-all-data)
  - [Nebulosa and Module Score analysis](#nebulosa-and-module-score-analysis)
  - [DEG Analysis](#deg-analysis-mta)
  - [Analysis of GRNs using pySCENIC](#analysis-of-grns-using-pyscenic)
  - [Inference of LR interactions using CellPhoneDB](#inference-of-lr-interactions-using-cellphonedb)
  - [Subclustering of mouse epithelial clusters (focus on ameloblast lineage)](#subclustering-of-mouse-epithelium)

- [Human Tooth Atlas](#human-tooth-atlas)
  - [Healthy HTA Setup](#healthy-hta-setup)
  - [Healthy HTA Quality Control](#healthy-hta-quality-control)
  - [Healthy HTA Removal of background or ambient RNA using SoupX](#healthy-hta-removal-of-background-or-ambient-rna-using-soupx)
  - [Healthy HTA rPCA integration](#healthy-hta-rpca-integration)
  - [Healthy HTA subclustering of dental epithelium](#healthy-hta-subclustering-of-dental-epithelium)



# Mouse Tooth Atlas
## MTA setup

```
suppressMessages({
#set working directory
setwd("./MTA")
getwd()

#load in required software packages
suppressWarnings({
library(reticulate)
library(dplyr)
library(Seurat)
library(SingleCellExperiment, quietly = TRUE)
library(scater, quietly = TRUE)
library(tidyr)
library(purrr)
library(cowplot)
library(ggrepel)
library(viridis)
library(gridExtra)
library(biomaRt)
library(schex)
library(Nebulosa)
library(dittoSeq)
library(escape)
library(SoupX)
library(monocle3) 
})

# choose random seed for reproducibility
set.seed(27011995)

})
```

Import all datasets, and add relevant metadata obtained from the publications. 

```
#Import data
GSM4407907 <- Read10X(data.dir = "./Chiba")
Chiba <- CreateSeuratObject(counts = GSM4407907, project = "GSM4407907", min.cells = 3, min.features = 200)

GSM3393718 <- Read10X(data.dir = "./Takahashi")
Takahashi <- CreateSeuratObject(counts = GSM3393718, project = "GSM3393718", min.cells = 3, min.features = 200)

GSM3767568 <- Read10X(data.dir = "./Sharir")
Sharir <- CreateSeuratObject(counts = GSM3767568, project = "GSM3767568", min.cells = 3, min.features = 200)

FB1RKA6 <- Read10X(data.dir = "./Chen")
Chen <- CreateSeuratObject(counts = FB1RKA6, project = "FB1RKA6", min.cells = 3, min.features = 200)

FB1X3X8 <- Read10X(data.dir = "./Wen")
Wen <- CreateSeuratObject(counts = FB1X3X8, project = "FB1X3X8", min.cells = 3, min.features = 200)

GSM5121221 <- Read10X(data.dir = "./Chiba_2_epithelium")
Chiba_2_epithelium <- CreateSeuratObject(counts = GSM5121221, project = "GSM5121221", min.cells = 3, min.features = 200)

GSM5121222 <- Read10X(data.dir = "./Chiba_2_mesencyhme")
Chiba_2_mesencyhme <- CreateSeuratObject(counts = GSM5121222, project = "GSM5121222", min.cells = 3, min.features = 200)

GSM5140554 <- Read10X(data.dir = "./Nagata")
Nagata <- CreateSeuratObject(counts = GSM5140554, project = "GSM5140554", min.cells = 3, min.features = 200)
```

```
# Import of Krivanek mouse incisor data
data = read.table(gzfile("./GSM4365604_counts_incisor_10x.txt.gz"), sep=" ", header = TRUE, row.names = NULL)

names <- make.unique(data$row.names)
rownames(data) <- names
data <- data[,-1] # get rid of old names

#Fix cell names
colnames(data) <- gsub('X', '', colnames(data))
colnames(data) <- gsub('\\.', '-', colnames(data))

#Fix gene formatting
rownames(data) <- str_to_title(rownames(data))

#Create Seurat Object
Krivanek_incisor_10x <- CreateSeuratObject(counts = data, project = "GSM4365604", min.cells = 3, min.features = 200)
```

```
# Import of Krivanek mouse molar data 1
data = read.table(gzfile("./GSM4365605_counts_molar_10x.1.txt.gz"), sep=" ", header = TRUE, row.names = NULL)

names <- make.unique(data$row.names)
rownames(data) <- names
data <- data[,-1] # get rid of old names

#Fix cell names
colnames(data) <- gsub('X', '', colnames(data))
colnames(data) <- gsub('\\.', '-', colnames(data))

#Fix gene formatting
rownames(data) <- str_to_title(rownames(data))

#Create Seurat Object
Krivanek_molar_10x_1 <- CreateSeuratObject(counts = data, project = "GSM4365605", min.cells = 3, min.features = 200)
```

```
# Import of Krivanek mouse molar data 2
data = read.table(gzfile("./GSM4365611_counts_molar_10x.2.txt.gz"), sep=" ", header = TRUE, row.names = NULL)

names <- make.unique(data$row.names)
rownames(data) <- names
data <- data[,-1] # get rid of old names

#Fix cell names
colnames(data) <- gsub('X', '', colnames(data))
colnames(data) <- gsub('\\.', '-', colnames(data))

#Fix gene formatting
rownames(data) <- str_to_title(rownames(data))

#Create Seurat Object
Krivanek_molar_10x_2 <- CreateSeuratObject(counts = data, project = "GSM4365611", min.cells = 3, min.features = 200)
```

```
# Import of Zhao periodontal tissue 1
data = read.table(gzfile("./GSM4872144_Data1.txt.gz"), sep="\t", header = TRUE, row.names = NULL)

names <- make.unique(data$row.names)
rownames(data) <- names
data <- data[,-1] # get rid of old names

#Fix cell names
colnames(data) <- gsub('SHAM.CTRL_', '', colnames(data))

#Create Seurat Object
Zhao_perio_1 <- CreateSeuratObject(counts = data, project = "GSM4872144", min.cells = 3, min.features = 200)
```

```
# Import of Zhao periodontal tissue 2
data2 = read.table(gzfile("./GSM4872145_Data2.txt.gz"), sep="\t", header = TRUE, row.names = NULL)

names <- make.unique(data2$row.names)
rownames(data2) <- names
data2 <- data2[,-1] # get rid of old names

#Fix cell names
colnames(data2) <- gsub('L120b_', '', colnames(data2))

#Create Seurat Object
Zhao_perio_2 <- CreateSeuratObject(counts = data2, project = "GSM4872145", min.cells = 3, min.features = 200)
```

```
# Calculate % of mitochondrial genes per cell, and append this to the metadata
Chiba[["percent.mito"]] <- PercentageFeatureSet(Chiba, pattern = "^mt-")
Takahashi[["percent.mito"]] <- PercentageFeatureSet(Takahashi, pattern = "^mt-")
Sharir[["percent.mito"]] <- PercentageFeatureSet(Sharir, pattern = "^mt-")
Chen[["percent.mito"]] <- PercentageFeatureSet(Chen, pattern = "^mt-")
Wen[["percent.mito"]] <- PercentageFeatureSet(Wen, pattern = "^mt-")
Chiba_2_epithelium[["percent.mito"]] <- PercentageFeatureSet(Chiba_2_epithelium, pattern = "^mt-")
Chiba_2_mesencyhme[["percent.mito"]] <- PercentageFeatureSet(Chiba_2_mesencyhme, pattern = "^mt-")
Nagata[["percent.mito"]] <- PercentageFeatureSet(Nagata, pattern = "^mt-")
Krivanek_incisor_10x[["percent.mito"]] <- PercentageFeatureSet(Krivanek_incisor_10x, pattern = "^Mt-")
Krivanek_molar_10x_1[["percent.mito"]] <- PercentageFeatureSet(Krivanek_molar_10x_1, pattern = "^Mt-")
Krivanek_molar_10x_2[["percent.mito"]] <- PercentageFeatureSet(Krivanek_molar_10x_2, pattern = "^Mt-")
Zhao_perio_1[["percent.mito"]] <- PercentageFeatureSet(Zhao_perio_1, pattern = "^mt-")
Zhao_perio_2[["percent.mito"]] <- PercentageFeatureSet(Zhao_perio_2, pattern = "^mt-")
```

```
# Add extra metadata to Seurat object
Chiba@meta.data$Dataset <- "Chiba_1"
Chiba@meta.data$Age <- "PD7"
Chiba@meta.data$Tooth_type <- "Incisors"
Chiba@meta.data$Cell_source <- "Epithelial_and_mesenchyme"
Chiba@meta.data$Genotype <- "Krt14-RFP"
Chiba@meta.data$Technology <- "10X"

Takahashi@meta.data$Dataset <- "Takahashi"
Takahashi@meta.data$Age <- "PD6"
Takahashi@meta.data$Tooth_type <- "Molars"
Takahashi@meta.data$Cell_source <- "PTHrP-mCherry+_DF_cells"
Takahashi@meta.data$Genotype <- "PTHrP-mCherry"
Takahashi@meta.data$Technology <- "10X"

Sharir@meta.data$Dataset <- "Sharir"
Sharir@meta.data$Age <- "8-12weeks"
Sharir@meta.data$Tooth_type <- "Incisors"
Sharir@meta.data$Cell_source <- "Sorted_epithelial_cells"
Sharir@meta.data$Genotype <- "C57BL/6J"
Sharir@meta.data$Technology <- "10X"

Chen@meta.data$Dataset <- "Chen"
Chen@meta.data$Age <- "4weeks"
Chen@meta.data$Tooth_type <- "Incisors"
Chen@meta.data$Cell_source <- "Epithelial_and_mesenchyme"
Chen@meta.data$Genotype <- "C57BL/6J"
Chen@meta.data$Technology <- "10X"

Wen@meta.data$Dataset <- "Wen"
Wen@meta.data$Age <- "PD7.5"
Wen@meta.data$Tooth_type <- "Molars"
Wen@meta.data$Cell_source <- "Epithelial_and_mesenchyme"
Wen@meta.data$Genotype <- "Runx2_fl/fl_controls"
Wen@meta.data$Technology <- "10X"

Chiba_2_epithelium@meta.data$Dataset <- "Chiba_2_epithelium"
Chiba_2_epithelium@meta.data$Age <- "PD1"
Chiba_2_epithelium@meta.data$Tooth_type <- "Molars"
Chiba_2_epithelium@meta.data$Cell_source <- "Epithelium"
Chiba_2_epithelium@meta.data$Genotype <- "Krt14-RFP"
Chiba_2_epithelium@meta.data$Technology <- "10X"

Chiba_2_mesencyhme@meta.data$Dataset <- "Chiba_2_mesenchyme"
Chiba_2_mesencyhme@meta.data$Age <- "PD1"
Chiba_2_mesencyhme@meta.data$Tooth_type <- "Molars"
Chiba_2_mesencyhme@meta.data$Cell_source <- "Mesencyhme"
Chiba_2_mesencyhme@meta.data$Genotype <- "Krt14-RFP"
Chiba_2_mesencyhme@meta.data$Technology <- "10X"

Nagata@meta.data$Dataset <- "Nagata"
Nagata@meta.data$Age <- "PD25"
Nagata@meta.data$Tooth_type <- "Molars"
Nagata@meta.data$Cell_source <- "Periodontal"
Nagata@meta.data$Genotype <- "Col1a1(2.3kb)-GFP; PTHrP-creER; R26R(tdTomato/+)"
Nagata@meta.data$Technology <- "10X"

Krivanek_incisor_10x@meta.data$orig.ident <- NULL
Krivanek_incisor_10x@meta.data$orig.ident <- "GSM4365604"
Krivanek_incisor_10x@meta.data$Dataset <- "Krivanek_incisor"
Krivanek_incisor_10x@meta.data$Age <- "8-16weeks"
Krivanek_incisor_10x@meta.data$Tooth_type <- "Incisors"
Krivanek_incisor_10x@meta.data$Cell_source <- "Epithelial_and_mesenchyme"
Krivanek_incisor_10x@meta.data$Genotype <- "C57BL/6J_and_Sox2eGFP"
Krivanek_incisor_10x@meta.data$Technology <- "10X"

Krivanek_molar_10x_1@meta.data$orig.ident <- NULL
Krivanek_molar_10x_1@meta.data$orig.ident <- "GSM4365605"
Krivanek_molar_10x_1@meta.data$Dataset <- "Krivanek_molar"
Krivanek_molar_10x_1@meta.data$Age <- "8-16weeks"
Krivanek_molar_10x_1@meta.data$Tooth_type <- "Molars"
Krivanek_molar_10x_1@meta.data$Cell_source <- "Epithelial_and_mesenchyme"
Krivanek_molar_10x_1@meta.data$Genotype <- "C57BL/6J_and_Sox2eGFP"
Krivanek_molar_10x_1@meta.data$Technology <- "10X"

Krivanek_molar_10x_2@meta.data$orig.ident <- NULL
Krivanek_molar_10x_2@meta.data$orig.ident <- "GSM4365605"
Krivanek_molar_10x_2@meta.data$Dataset <- "Krivanek_molar"
Krivanek_molar_10x_2@meta.data$Age <- "8-16weeks"
Krivanek_molar_10x_2@meta.data$Tooth_type <- "Molars"
Krivanek_molar_10x_2@meta.data$Cell_source <- "Epithelial_and_mesenchyme"
Krivanek_molar_10x_2@meta.data$Genotype <- "C57BL/6J_and_Sox2eGFP"
Krivanek_molar_10x_2@meta.data$Technology <- "10X"

Zhao_perio_1@meta.data$Dataset <- "Zhao_perio_1"
Zhao_perio_1@meta.data$Age <- "Adult"
Zhao_perio_1@meta.data$Tooth_type <- "Molars"
Zhao_perio_1@meta.data$Cell_source <- "Periodontal"
Zhao_perio_1@meta.data$Genotype <- "CD1"
Zhao_perio_1@meta.data$Technology <- "10X"

Zhao_perio_2@meta.data$Dataset <- "Zhao_perio_2"
Zhao_perio_2@meta.data$Age <- "Adult"
Zhao_perio_2@meta.data$Tooth_type <- "Molars"
Zhao_perio_2@meta.data$Cell_source <- "Periodontal"
Zhao_perio_2@meta.data$Genotype <- "CD1"
Zhao_perio_2@meta.data$Technology <- "10X"
```

```
# Add Dataset of origin to cell name to avoid possible identical cell names between datasets
Chiba <- RenameCells(Chiba, new.names = paste0(Chiba$Dataset, "_", Cells(Chiba)))
Takahashi <- RenameCells(Takahashi, new.names = paste0(Takahashi$Dataset, "_", Cells(Takahashi)))
Sharir <- RenameCells(Sharir, new.names = paste0(Sharir$Dataset, "_", Cells(Sharir)))
Chen <- RenameCells(Chen, new.names = paste0(Chen$Dataset, "_", Cells(Chen)))
Wen <- RenameCells(Wen, new.names = paste0(Wen$Dataset, "_", Cells(Wen)))
Chiba_2_epithelium <- RenameCells(Chiba_2_epithelium, new.names = paste0(Chiba_2_epithelium$Dataset, "_", Cells(Chiba_2_epithelium)))
Chiba_2_mesencyhme <- RenameCells(Chiba_2_mesencyhme, new.names = paste0(Chiba_2_mesencyhme$Dataset, "_", Cells(Chiba_2_mesencyhme)))
Nagata <- RenameCells(Nagata, new.names = paste0(Nagata$Dataset, "_", Cells(Nagata)))
Krivanek_incisor_10x <- RenameCells(Krivanek_incisor_10x, new.names = paste0(Krivanek_incisor_10x$Dataset, "_", Cells(Krivanek_incisor_10x)))
Krivanek_molar_10x_1 <- RenameCells(Krivanek_molar_10x_1, new.names = paste0(Krivanek_molar_10x_1$Dataset, "_", Cells(Krivanek_molar_10x_1)))
Krivanek_molar_10x_2 <- RenameCells(Krivanek_molar_10x_2, new.names = paste0(Krivanek_molar_10x_2$Dataset, "_", Cells(Krivanek_molar_10x_2)))
Zhao_perio_1 <- RenameCells(Zhao_perio_1, new.names = paste0(Zhao_perio_1$Dataset, "_", Cells(Zhao_perio_1)))
Zhao_perio_2 <- RenameCells(Zhao_perio_2, new.names = paste0(Zhao_perio_2$Dataset, "_", Cells(Zhao_perio_2)))
```


## MTA quality control
```
Krivanek_molars <- merge(x = Krivanek_molars_1, y = Krivanek_molars_2) #merged now because the cell number is too small for downstream steps otherwise 

# Perform QC individually
Chiba <- subset(Chiba, subset = nFeature_RNA > 1500 & nFeature_RNA < 6000 & percent.mito < 5 & nCount_RNA < 40000)
Takahashi <- subset(Takahashi, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mito < 20 & nCount_RNA < 40000)
Sharir <- subset(Sharir, subset = nFeature_RNA > 1000 & nFeature_RNA < 4500 & percent.mito < 5 & nCount_RNA < 40000)
Krivanek_incisors <- subset(Krivanek_incisors, subset = nFeature_RNA > 750 & nFeature_RNA < 2000 & percent.mito < 15 & nCount_RNA < 40000)
Krivanek_molars <- subset(Krivanek_molars, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mito < 8 & nCount_RNA < 40000)
Chen <- subset(Chen, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mito < 20 & nCount_RNA < 40000)
Wen <- subset(Wen, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mito < 20 & nCount_RNA < 40000)
Chiba_2_epi <- subset(Chiba_2_epi, subset = nFeature_RNA > 1000 & nFeature_RNA < 4500 & percent.mito < 5 & nCount_RNA < 40000)
Chiba_2_mes <- subset(Chiba_2_mes, subset = nFeature_RNA > 1000 & nFeature_RNA < 4500 & percent.mito < 5 & nCount_RNA < 40000)
Nagata <- subset(Nagata, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mito < 20 & nCount_RNA < 40000)
Zhao_1 <- subset(Zhao_1, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mito < 8 & nCount_RNA < 40000) 
Zhao_2 <- subset(Zhao_2, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mito < 8 & nCount_RNA < 40000) 
```

```
# Merge datasets
merged_final_qc <- merge(x = Chiba, y = c(Takahashi, Sharir, Krivanek_incisors, Krivanek_molars, Chen, Wen, Chiba_2_epi, Chiba_2_mes, Nagata, Zhao_1, Zhao_2))
```

## Initial integration on three groups 
(incisors, molars, periodontal), using 'standard Seurat integration workflow'

```
# Merge into the new objects
incisors <- merge(x = Chiba, y = c(Sharir, Krivanek_incisor, Chen))
molars <- merge(x = Krivanek_molar, y = c(Wen, Chiba_2_epithelium, Chiba_2_mesenchyme))
periodontal <- merge(x = Takahashi, y = c(Nagata, Zhao_perio_1, Zhao_perio_2))
```

```
# After creating the objects for each group, we first perform standard normalization and variable feature selectio on the list of objects
incisors.list <- SplitObject(incisors, split.by = "Dataset")
incisors.list <- lapply(X = incisors.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

molars.list <- SplitObject(molars, split.by = "Dataset")
molars.list <- lapply(X = molars.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

periodontal.list <- SplitObject(periodontal, split.by = "Dataset")
periodontal.list <- lapply(X = periodontal.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})
```  
  
```
# Perform integration
# Perform integration
incisor.anchors <- FindIntegrationAnchors(object.list = incisors.list, dims = 1:30)
incisors_integrated <- IntegrateData(anchorset = incisor.anchors, dims = 1:30)

molars.anchors <- FindIntegrationAnchors(object.list = molars.list, dims = 1:30)
molars_integrated <- IntegrateData(anchorset = molars.anchors, dims = 1:30)

periodontal.anchors <- FindIntegrationAnchors(object.list = periodontal.list, dims = 1:30)
periodontal_integrated <- IntegrateData(anchorset = periodontal.anchors, dims = 1:30)
```

```
# switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(incisors_integrated) <- "integrated"
DefaultAssay(molars_integrated) <- "integrated"
DefaultAssay(periodontal_integrated) <- "integrated"
```

```
# Run the standard workflow for visualization and clustering
incisors_integrated <- ScaleData(incisors_integrated, verbose = FALSE)
incisors_integrated <- RunPCA(incisors_integrated, verbose = FALSE)

molars_integrated <- ScaleData(molars_integrated, verbose = FALSE)
molars_integrated <- RunPCA(molars_integrated, verbose = FALSE)

periodontal_integrated <- ScaleData(periodontal_integrated, verbose = FALSE)
periodontal_integrated <- RunPCA(periodontal_integrated, verbose = FALSE)
```

``` 
# Perform UMAP dimensional reduction using umap-learn
incisors_integrated <- RunUMAP(incisors_integrated, umap.method = "umap-learn", dims = 1:30)
molars_integrated <- RunUMAP(molars_integrated, umap.method = "umap-learn", dims = 1:30)
periodontal_integrated <- RunUMAP(periodontal_integrated, umap.method = "umap-learn", dims = 1:30)
```

```
# Perform subclustering on each group for rough annotation
DefaultAssay(incisors_integrated) <- "integrated" #needs to be set to integrated because this is what we ran PCA on
incisors_integrated <- FindNeighbors(incisors_integrated, dims = 1:30)
incisors_integrated <- FindClusters(incisors_integrated, resolution = c(0.1, 0.2, 0.3, 0.6, 0.8, 1))

DefaultAssay(molars_integrated) <- "integrated" #needs to be set to integrated because this is what we ran PCA on
molars_integrated <- FindNeighbors(molars_integrated, dims = 1:30)
molars_integrated <- FindClusters(molars_integrated, resolution = c(0.1, 0.2, 0.3, 0.6, 0.8, 1))

DefaultAssay(periodontal_integrated) <- "integrated" #needs to be set to integrated because this is what we ran PCA on
periodontal_integrated <- FindNeighbors(periodontal_integrated, dims = 1:30)
periodontal_integrated <- FindClusters(periodontal_integrated, resolution = c(0.1, 0.2, 0.3, 0.6, 0.8, 1))
```

```
# Check different resolutions for each group
resolution <- c("integrated_snn_res.0.1", "integrated_snn_res.0.2", "integrated_snn_res.0.3", "integrated_snn_res.0.6", "integrated_snn_res.0.8","integrated_snn_res.1")

suppressMessages({
plot.list <- list()
for (i in 1:length(resolution)) {
    options(repr.plot.width=7, repr.plot.height=5)
    plot.list[[i]] <- DimPlot(incisors_integrated, reduction = "umap", group.by = resolution[i], label = TRUE) +
        ggtitle(paste(resolution[i], sep = ""))
    print(plot.list[[i]])
}
    
for (i in 1:length(resolution)) {
    options(repr.plot.width=7, repr.plot.height=5)
    plot.list[[i]] <- DimPlot(molars_integrated, reduction = "umap", group.by = resolution[i], label = TRUE) +
        ggtitle(paste(resolution[i], sep = ""))
    print(plot.list[[i]])
}
    
for (i in 1:length(resolution)) {
    options(repr.plot.width=7, repr.plot.height=5)
    plot.list[[i]] <- DimPlot(periodontal_integrated, reduction = "umap", group.by = resolution[i], label = TRUE) +
        ggtitle(paste(resolution[i], sep = ""))
    print(plot.list[[i]])
}   
})
```

``` 
# Select resolution for each group, and annotate based on marker expression
Idents(incisors_integrated) <- incisors_integrated@meta.data$integrated_snn_res.0.3
Idents(molars_integrated) <- molars_integrated@meta.data$integrated_snn_res.0.6
Idents(periodontal_integrated) <- periodontal_integrated@meta.data$integrated_snn_res.0.8
```


```
incisors_integrated <- RenameIdents(incisors_integrated, 
                           `0` = "Epithelial", `1` = "Dental Follicle", `2` = "Pulp", `3` = "Secretory Ameloblasts", 
                           `4` = "Macrophages", `5` = "Cycling cells", `6` = "Pre-Ameloblasts", `7` = "Odontoblasts", `8` = "Endothelial", 
                           `9` = "Epithelial", `10` = "Immune cells", `11` = "Perivascular", `12` = "RBCs", 
                           `13` = "Neutrophils", `14` = "Schwann/Oligo")

molars_integrated <- RenameIdents(molars_integrated, 
                           `0` = "Pre-Ameloblasts", `1` = "Pulp", `2` = "Pulp", `3` = "Pre-Ameloblasts", 
                           `4` = "Epithelial", `5` = "Dental Follicle", `6` = "Cycling Epithelial", `7` = "Cycling Epithelial", `8` = "Epithelial", 
                           `9` = "Endothelial", `10` = "Macrophages", `11` = "Odontoblasts", `12` = "Epithelial", 
                           `13` = "Cycling Pulp", `14` = "Secretory Ameloblasts", `15` = "Macrophages", `16` = "Schwann/Oligo", 
                           `17` = "Perivascular", `18` = "Dental Follicle", `19` = "Glial", `20` = "Dental Follicle", 
                            `21` = "Immune cells", `22` = "Maturation-stage Ameloblasts", `23` = "Immune cells")

periodontal_integrated <- RenameIdents(periodontal_integrated, 
                       `0` = "Dental Follicle", `1` = "Dental Follicle", `2` = "Endothelial", `3` = "Epithelial", 
                        `4` = "Dental Follicle", `5` = "Macrophages", `6` = "Endothelial", `7` = "Pulp", `8` = "Macrophages", 
                        `9` = "Pulp", `10` = "Cementoblasts", `11` = "PDL", `12` = "Neutrophils", 
                        `13` = "Immune cells", `14` = "Immune cells", `15` = "Osteoblasts", `16` = "Macrophages", 
                        `17` = "Macrophages", `18` = "Pulp", `19` = "Neutrophils", `20` = "Immune cells", 
                        `21` = "Immune cells", `22` = "RBCs", `23` = "Perivascular", `24` = "RBCs",
                        `25` = "Schwann/Oligo", `26` = "Macrophages")   

incisors_integrated$Split_Clusters_2 <- Idents(incisors_integrated)
molars_integrated$Split_Clusters_2 <- Idents(molars_integrated)
periodontal_integrated$Split_Clusters_2 <- Idents(periodontal_integrated)
```

## Removal of background or ambient RNA using SoupX
```
#Change from integrated to RNA, otherwise next command will fail
DefaultAssay(incisors_integrated) <- "RNA"
DefaultAssay(molars_integrated) <- "RNA"
DefaultAssay(periodontal_integrated) <- "RNA"
```

```
all.raw.data_incisors <- as.matrix(GetAssayData(incisors_integrated, slot = "counts"))
list_cluster_incisor <- incisors_integrated@meta.data[[sprintf("Split_Clusters")]]
names(list_cluster_incisor) <- incisors_integrated@assays[["RNA"]]@data@Dimnames[[2]]
list_cluster_incisor <- as.matrix(list_cluster_incisor)
umap_incisor <- (Embeddings(incisors_integrated, reduction = "umap"))


all.raw.data_molars <- as.matrix(GetAssayData(molars_integrated, slot = "counts"))
list_cluster_molar <- molars_integrated@meta.data[[sprintf("Split_Clusters")]]
names(list_cluster_molar) <- molars_integrated@assays[["RNA"]]@data@Dimnames[[2]]
list_cluster_molar <- as.matrix(list_cluster_molar)
umap_molar <- (Embeddings(molars_integrated, reduction = "umap"))


all.raw.data_periodontal <- as.matrix(GetAssayData(periodontal_integrated, slot = "counts"))
list_cluster_periodontal <- periodontal_integrated@meta.data[[sprintf("Split_Clusters")]]
names(list_cluster_periodontal) <- periodontal_integrated@assays[["RNA"]]@data@Dimnames[[2]]
list_cluster_periodontal <- as.matrix(list_cluster_periodontal)
umap_periodontal <- (Embeddings(periodontal_integrated, reduction = "umap"))
```

```
# Create a SoupChannel object; because we are starting from our Seurat object we do not have a file with empty droplets. We will import the raw counts matrix from Seurat both as the table of counts (toc) and table of droplets (tod)

toc_incisor <- all.raw.data_incisors
tod_incisor = all.raw.data_incisors
sc_incisor = SoupChannel(tod_incisor, toc_incisor, calcSoupProfile = FALSE)

toc_molar <- all.raw.data_molars
tod_molar = all.raw.data_molars
sc_molar = SoupChannel(tod_molar, toc_molar, calcSoupProfile = FALSE)

toc_periodontal <- all.raw.data_periodontal
tod_periodontal = all.raw.data_periodontal
sc_periodontal = SoupChannel(tod_periodontal, toc_periodontal, calcSoupProfile = FALSE)
```

```
# Calculate the Soup profile (of ambient RNA)

toc = sc_incisor$toc
scNoDrops = SoupChannel(toc, toc, calcSoupProfile = FALSE)
soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), 
    counts = rowSums(toc))
scNoDrops_incisor = setSoupProfile(scNoDrops, soupProf)

toc = sc_molar$toc
scNoDrops = SoupChannel(toc, toc, calcSoupProfile = FALSE)
soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), 
    counts = rowSums(toc))
scNoDrops_molar = setSoupProfile(scNoDrops, soupProf)

toc = sc_periodontal$toc
scNoDrops = SoupChannel(toc, toc, calcSoupProfile = FALSE)
soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), 
    counts = rowSums(toc))
scNoDrops_periodontal = setSoupProfile(scNoDrops, soupProf)
```

```
# Add UMAP coordinates and cluster information from Seurat analysis to the SoupChannel object
scNoDrops_incisor = setClusters(scNoDrops_incisor, setNames(list_cluster_incisor[,1], rownames(list_cluster_incisor)))
scNoDrops_incisor = setDR(scNoDrops_incisor, umap_incisor[,c('UMAP_1','UMAP_2')])

scNoDrops_molar = setClusters(scNoDrops_molar, setNames(list_cluster_molar[,1], rownames(list_cluster_molar)))
scNoDrops_molar = setDR(scNoDrops_molar, umap_molar[,c('UMAP_1','UMAP_2')])

scNoDrops_periodontal = setClusters(scNoDrops_periodontal, setNames(list_cluster_periodontal[,1], rownames(list_cluster_periodontal)))
scNoDrops_periodontal = setDR(scNoDrops_periodontal, umap_periodontal[,c('UMAP_1','UMAP_2')])
```

```
# Estimate the contamination fraction
options(repr.plot.width=7, repr.plot.height=7)
sc_incisor = autoEstCont(scNoDrops_incisor, verbose = TRUE)
sc_molar = autoEstCont(scNoDrops_molar, verbose = TRUE)
sc_periodontal = autoEstCont(scNoDrops_periodontal, verbose = TRUE)
```

```
# Remove the calculated contamination fraction from the original counts matrix, and add back to the original Seurat object. 
out_incisor = adjustCounts(sc_incisor)
incisors_integrated[["SoupX"]] <- CreateAssayObject(counts = out_incisor)

out_molar = adjustCounts(sc_molar)
molars_integrated[["SoupX"]] <- CreateAssayObject(counts = out_molar)

out_periodontal = adjustCounts(sc_periodontal)
periodontal_integrated[["SoupX"]] <- CreateAssayObject(counts = out_periodontal)
```

```
# Plot the change in expression due to the correction, e.g. for Ambn and Amelx
options(repr.plot.width=7, repr.plot.height=7)
plotChangeMap(sc_molar, out_molar, "Amelx", pointSize = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

plotChangeMap(sc_incisor, out_incisor, "Amelx", pointSize = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

plotChangeMap(sc_periodontal, out_periodontal, "Amelx", pointSize = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

plotChangeMap(sc_molar, out_molar, "Ambn", pointSize = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

plotChangeMap(sc_incisor, out_incisor, "Ambn", pointSize = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

plotChangeMap(sc_periodontal, out_periodontal, "Ambn", pointSize = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
```

```
# Normalize corrected counts, find variable features and scale data
DefaultAssay(incisors_integrated) <- "SoupX"
incisors_integrated <- NormalizeData(incisors_integrated)
incisors_integrated <- FindVariableFeatures(incisors_integrated, selection.method = "vst", nfeatures = 2000)
incisors_integrated <- ScaleData(incisors_integrated, verbose = FALSE)

DefaultAssay(molars_integrated) <- "SoupX"
molars_integrated <- NormalizeData(molars_integrated)
molars_integrated <- FindVariableFeatures(molars_integrated, selection.method = "vst", nfeatures = 2000)
molars_integrated <- ScaleData(molars_integrated, verbose = FALSE)

DefaultAssay(periodontal_integrated) <- "SoupX"
periodontal_integrated <- NormalizeData(periodontal_integrated)
periodontal_integrated <- FindVariableFeatures(periodontal_integrated, selection.method = "vst", nfeatures = 2000)
periodontal_integrated <- ScaleData(periodontal_integrated, verbose = FALSE)
```

``` 
# Merge the three objects with correct counts matrices for further analysis
merged <- merge(x =incisors_integrated, y = c(molars_integrated, periodontal_integrated))
```

## rPCA integration of all data
Integration of 3 groups and cell cycle regression using Seurat's reciprocal PCA (rPCA) workflow 
```
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can segregate this list into markers of G2/M phase and markers of S phase. 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# These lists of S and G2M phase genes contain human gene names, so they must be converted to their mouse orthologues first. 
#Basic function to convert human to mouse gene names using the biomaRt package
convertHumanGeneList <- function(x){
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
    humanx <- unique(genesV2[, 2])
    # Print the first 6 genes found to the screen
    print(head(humanx))
    return(humanx)
    }

s.genes <- convertHumanGeneList(s.genes)
g2m.genes <- convertHumanGeneList(g2m.genes)
```

```
# After acquiring the data, we first perform standard normalization and variable feature selection on the list of objects
tooth.list <- SplitObject(merged, split.by = "Dataset")
tooth.list <- lapply(X = tooth.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
})
```

```
# select features that are repeatedly variable across datasets for integration run PCA on each dataset using these features
features <- SelectIntegrationFeatures(object.list = tooth.list)
tooth.list <- lapply(X = tooth.list, FUN = function(x) {
    x <- ScaleData(x, features = features, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
```

```
# Perform integration
tooth.anchors <- FindIntegrationAnchors(object.list = tooth.list, dims = 1:30, anchor.features = features, reduction = "rpca")
integrated <- IntegrateData(anchorset = tooth.anchors, dims = 1:30)
```

```
# switch to integrated assay. The variable features of this assay are automaticallyset during IntegrateData
DefaultAssay(integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, verbose = FALSE)
```

```
# Run UMAP for different dimensions, save to object separately to be able to go back
integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:15, min.dist = 0.5)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_15dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP15_", assay = DefaultAssay(integrated))

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:20, min.dist = 0.5)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_20dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP20_", assay = DefaultAssay(integrated))

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:30, min.dist = 0.5)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_30dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP30_", assay = DefaultAssay(integrated))

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:40, min.dist = 0.5)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_40dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP40_", assay = DefaultAssay(integrated))

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:50, min.dist = 0.5)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_50dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP50_", assay = DefaultAssay(integrated))
```

```
options(repr.plot.width=14, repr.plot.height=7)
p1 <- DimPlot(integrated, reduction = "umap_15dims", group.by = "Dataset") +ggtitle("UMAP 15dims")
p2 <- DimPlot(integrated, reduction = "umap_15dims", group.by = "Age", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_20dims", group.by = "Dataset") +ggtitle("UMAP 20dims")
p2 <- DimPlot(integrated, reduction = "umap_20dims", group.by = "Age", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_30dims", group.by = "Dataset") +ggtitle("UMAP 30dims")
p2 <- DimPlot(integrated, reduction = "umap_30dims", group.by = "Age", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_40dims", group.by = "Dataset") +ggtitle("UMAP 40dims")
p2 <- DimPlot(integrated, reduction = "umap_40dims", group.by = "Age", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_50dims", group.by = "Dataset") +ggtitle("UMAP 50dims")
p2 <- DimPlot(integrated, reduction = "umap_50dims", group.by = "Age", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)
```

```
options(repr.plot.width=14, repr.plot.height=7)
p1 <- DimPlot(integrated, reduction = "umap_15dims", group.by = "Phase") +ggtitle("UMAP 15dims")
p2 <- DimPlot(integrated, reduction = "umap_15dims", group.by = "Tooth_type", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_20dims", group.by = "Phase") +ggtitle("UMAP 20dims")
p2 <- DimPlot(integrated, reduction = "umap_20dims", group.by = "Tooth_type", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_30dims", group.by = "Phase") +ggtitle("UMAP 30dims")
p2 <- DimPlot(integrated, reduction = "umap_30dims", group.by = "Tooth_type", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_40dims", group.by = "Phase") +ggtitle("UMAP 40dims")
p2 <- DimPlot(integrated, reduction = "umap_40dims", group.by = "Tooth_type", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_50dims", group.by = "Phase") +ggtitle("UMAP 50dims")
p2 <- DimPlot(integrated, reduction = "umap_50dims", group.by = "Tooth_type", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)
```

We continue with 50 dimensions. 

```
# Perform subclustering, assess different resolutions 
DefaultAssay(integrated) <- "integrated" #needs to be set to integrated because this is what we ran PCA on
integrated <- FindNeighbors(integrated, dims = 1:50) # set to same as used for dimensonal reduction
integrated <- FindClusters(integrated, resolution = c(0.1, 0.2, 0.3, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 5, 10))
```

```
# Visualize clustering at different resolutions
resolution <- c("integrated_snn_res.0.1", "integrated_snn_res.0.2", "integrated_snn_res.0.3", "integrated_snn_res.0.6", "integrated_snn_res.0.8", 
                "integrated_snn_res.1", "integrated_snn_res.1.2", "integrated_snn_res.1.4", "integrated_snn_res.1.6", "integrated_snn_res.1.8", "integrated_snn_res.2",
               "integrated_snn_res.5", "integrated_snn_res.10")

suppressMessages({
plot.list <- list()
for (i in 1:length(resolution)) {
    options(repr.plot.width=7, repr.plot.height=5)
    plot.list[[i]] <- DimPlot(integrated, reduction = "umap_50dims", group.by = resolution[i], label = TRUE) +
        ggtitle(paste(resolution[i], sep = ""))
    print(plot.list[[i]])
}  
})
```

We continue with resolution = 5, and annotate based on expression of marker genes

```
Idents(integrated) <- integrated@meta.data$integrated_snn_res.5

integrated <- RenameIdents(integrated,
                           `0` = "Apical Pulp", `1` = "Apical Pulp",
                           `2` = "preAB", `3` = "VEE/OEE",
                           `4` = "DEP", `5` = "Endothelial",
                           `6` = "DF", `7` = "Distal Pulp",
                           `8` = "Early OBs", `9` = "Macrophage",
                           `10` = "SI", `11` = "DEP",
                           `12` = "sAB", `13` = "Macrophage",
                           `14` = "Early OBs", `15` = "Distal Pulp",
                           `16` = "Cycling Pulp", `17` = "VEE/OEE",
                           `18` = "Cycling DEP", `19` = "preAB",
                           `20` = "DEP", `21` = "Cycling DEP",
                           `22` = "DF", `23` = "DF",
                           `24` = "Endothelial", `25` = "DF",
                           `26` = "Macrophage", `27` = "PDL",
                           `28` = "DEP", `29` = "sAB",
                           `30` = "DF", `31` = "Cycling DEP",
                           `32` = "PDL", `33` = "mOB",
                           `34` = "IEE/OEE", `35` = "Endothelial",
                           `36` = "Cycling DF", `37` = "DEP",
                           `38` = "DEP", `39` = "Cycling Pulp",
                           `40` = "IEE/OEE", `41` = "Cycling DEP",
                           `42` = "Distal Pulp", `43` = "B",
                           `44` = "Macrophage", `45` = "Low Quality",
                           `46` = "Cycling Pulp", `47` = "DF",
                           `48` = "Perivascular", `49` = "Macrophage",
                           `50` = "Macrophage", `51` = "Cycling DEP",
                           `52` = "DF", `53` = "NK",
                           `54` = "Macrophage", `55` = "sAB",
                           `56` = "Fibroblasts", `57` = "Cycling DEP",
                           `58` = "Neutrophil", `59` = "CB",
                           `60` = "PDL", `61` = "MSCs",
                           `62` = "sAB", `63` = "Endothelial",
                           `64` = "SR", `65` = "T",
                           `66` = "Perivascular", `67` = "Distal Pulp",
                           `68` = "Neutrophil", `69` = "Neutrophil",
                           `70` = "Cycling DEP", `71` = "Schwann/Oligo",
                           `72` = "RBC", `73` = "Cycling SI",
                           `74` = "Marrow Stromal", `75` = "Neutrophil",
                           `76` = "Endothelial", `77` = "Glial",
                           `78` = "PDL", `79` = "DEP",
                           `80` = "Osteoblast", `81` = "mAB",
                           `82` = "RBC", `83` = "Macrophage",
                           `84` = "pDC", `85` = "Neutrophil",
                           `86` = "Fibroblasts", `87` = "Distal Pulp",
                           `88` = "sAB", `89` = "Endothelial")
```

## Nebulosa and Module Score analysis
```
DefaultAssay(integrated) <- "SoupX"
options(repr.plot.width=12, repr.plot.height=12)
```

```
# Gene expression of DESC genes
DESC_features <- list(c('Sox2','Lgr5','Gli1','Lrig1','Bmi1','Ptch1','Pknox2','Zfp273','Spock1','Pcp4','Sfrp5'))

plot.list <- list()

for (i in 1:length(DESC_features)) {
    options(repr.plot.width=12, repr.plot.height=12)
    plot.list[[i]] <- FeaturePlot(integrated, reduction = "umap_50dims", features = DESC_features[i], pt.size = 1, cols = (c("lightgrey", "steelblue", "blue4")))
    print(plot.list[[i]])
}
```

```
# DESC Module with AddModuleScore (Seurat)

DESC_features <- list(c('Sox2','Lgr5','Gli1','Lrig1','Bmi1','Ptch1','Pknox2','Zfp273','Spock1','Pcp4','Sfrp5'))

integrated <- AddModuleScore(
  object = integrated,
  features = DESC_features,
  ctrl = 5,
  name = 'DESC_features',
    assay = "SoupX"
)

FeaturePlot(integrated, reduction = "umap_50dims", features = "DESC_features1", pt.size = 1, cols = (c("lightgrey", "steelblue", "blue4")))
```

```
# DESC (joint) densities with Nebulosa
p_list <- plot_density(seu, c('Sox2','Lgr5','Gli1','Lrig1','Bmi1','Ptch1','Pknox2','Zfp273','Spock1','Pcp4','Sfrp5'),
  joint = TRUE, combine = FALSE)
  
p_list[[length(p_list)]]
```

## DEG Analysis MTA
```
DefaultAssay(integrated) <- "SoupX"
```

```
all.markers <- FindAllMarkers(object = integrated,
                             logfc.threshold = 0.5, #default = 0.25
                             min.pct = 0.25) # default = 0.1
```                                                

## Analysis of GRNs using pySCENIC
First we extract the required data from our integrated Seurat object in our R environment. Afterwards we switch to our Python Conda environment containing pySCENIC. After setup, analysis can be divided into several steps: 
  (1) inference of co-expression modules
  (2) generation of regulons
  (3) Calculate AUC (i.e. cellular regulon enrichment)
  (4) Create a loom file (to share data via SCope, and to calculate RSS scores)
  (5) Calculate RSS scores

Afterwards, we move back to our R environment for further analyses using the pySCENIC output:
  (6) Create a regulon-based Seurat object and perform regulon-based integration
  (7) Calculate MRA

#### Extract required data
```
normCounts <- integrated[["SoupX"]]@data
write.csv(normCounts, "./normCounts_4scenic.csv")

annotation <- integrated@meta.data
write.csv(annotation, "./annotation_4scenic.csv")

umap_dims <- integrated@reductions$umap_50dims@cell.embeddings
write.csv(umap_dims, "./SCENIC/resources_folder")
```

#### Switch to the pySCENIC Conda environment and set up workspace
```
import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import seaborn as sns
```

```
cd /MTA/SCENIC/
```

```
DATA_FOLDER="/MTA/SCENIC/data_folder"
RESOURCES_FOLDER="/MTA/SCENIC/resources_folder"
DATABASE_FOLDER = "/MTA/database_folder"
SCHEDULER="123.122.8.24:8786"
DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "mm10*.mc9nr.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'mm_mgi_tfs.txt')
SC_EXP_FNAME = os.path.join(RESOURCES_FOLDER, "normCounts_4scenic.csv")
REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.p")
MOTIFS_FNAME = os.path.join(DATA_FOLDER, "motifs.csv")
```
NB: The motifs and feather files were obtained from https://resources.aertslab.org/cistarget/ 

NB2: The list of TFs can be found from our resources folder in the repository. xxx

#### Load in counts, TF and ranking databases. 
```
exp_mtx = pd.read_csv(os.path.join(RESOURCES_FOLDER, "normCounts_4scenic.csv"), index_col=0, sep=',').T
exp_mtx.shape #check file
exp_mtx.iloc[0:10,100:110] #check file
```
```
tf_names = load_tf_names(MM_TFS_FNAME)
```
```
db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
dbs
```

#### (1) inference of co-expression modules
```
from distributed import Client, LocalCluster
local_cluster = LocalCluster(n_workers=10, threads_per_worker=1)
client = Client(local_cluster)
```
```
adjacencies = grnboost2(
    expression_data=exp_mtx,
    tf_names=tf_names,
    verbose=True,
    client_or_address=client)
    
adjacencies.to_csv("adjacencies_norm.csv")
adjacencies.head()
```
```
adjacencies = pd.read_csv('adjacencies_norm.csv', index_col=0)
adjacencies.head()
```

#### (2) generation of regulons
```
# Generate modules based on calculated adjacencies
modules = list(modules_from_adjacencies(adjacencies, exp_mtx))
```
```
# Prune modules using cic-regulatory information (RcisTarget)
with ProgressBar():
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME, client_or_address=client)
    
regulons = df2regulons(df)
len(regulons)

df.to_csv(MOTIFS_FNAME)
```
```
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)
```

#### (3) Calculate AUC
```
regulons = df2regulons(df)
```
```
auc_mtx = aucell(exp_mtx, regulons, num_workers=8)
auc_mtx.to_csv("auc_mtx.csv")
```
```
auc_mtx = pd.read_csv('auc_mtx.csv', index_col=0)
auc_mtx.head()
```
```
auc_mtx.index = auc_mtx.index.str.replace('.', '-') #fix formatting of cell names
auc_mtx.index = auc_mtx.index.str.replace('X', '') #fix formatting of cell names, somehow and X was added to the cell names
auc_mtx.to_csv("auc_mtx_corrected_index.csv")
auc_mtx.head()
auc_mtx.shape
```

#### (4) Loom file generation
Set up environment
```
from pyscenic.export import export2loom
from pyscenic.utils import load_motifs, Sequence
from pyscenic.transform import df2regulons
```
```
# Import gene expression matrix, fix formatting of cell names.
exp_mtx = pd.read_csv(os.path.join(RESOURCES_FOLDER, "normCounts_4scenic.csv"), index_col=0, sep=',').T
exp_mtx.index = exp_mtx.index.str.replace('.', '-')
exp_mtx.head()
```
```
# Check if expression matrix has the correct format
def is_valid_exp_matrix(mtx):
    return (all(isinstance(idx, str) for idx in mtx.index) 
            and all(isinstance(idx, str) for idx in mtx.columns)
            and (mtx.index.nlevels == 1)
            and (mtx.columns.nlevels == 1))
is_valid_exp_matrix(exp_mtx)
```
Import all required data, and ensure correct formatting. 
```
# Import motifs
motifs = df
```
```
# Import metadata from Seurat
annotation_loom = pd.read_csv("/MTA/SCENIC/resources_folder/annotation_4scenic.csv", index_col=0, sep=',')
annotation_loom = annotation_loom.loc[:,"Annotation"] # Extract cell annotation
annotation_loom.index = annotation_loom.index.str.replace('.', '-') # fix formatting of cell names

annotation_loom.to_csv("_annotation_loom.csv")
!mv /MTA/SCENIC/annotation_loom.csv /MTA/SCENIC/resources_folder/
```
```
ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "FH_atlas_annotation_loom.csv")
with open(ANNOTATIONS_FNAME, "rt") as f:
     annotations = dict(line.strip().replace("\"", "").split(",") for idx, line in enumerate(f) if idx > 0)
```
```
#Check if annotation has the correct format
def is_valid_annotation_mapping(m):
    return (all(isinstance(k, str) for k in m.keys()) 
            and all(isinstance(v, str) for v in m.values()))
is_valid_annotation_mapping(annotations)
dict(list(annotations.items())[:5])
```
```
# Import UMAP dimensional reduction coordinates
embeddings = pd.read_csv("/MTA/SCENIC/resources_folder/UMAP_dims_4scenic.csv", index_col=0, sep=',')
embeddings.index = embeddings.index.str.replace('.', '-') # fix formatting of cell names
```
```
embeddings = embeddings.iloc[:,[0,1]]
embeddings.columns=['_X', '_Y']
```
```
embeddings = { "UMAP (default)" : pd.DataFrame(data=embeddings,
                                      index=exp_mtx.index, columns=['_X', '_Y']) } # (n_cells, 2)
                                      
```
```
from collections import OrderedDict
id2name = OrderedDict()
embeddings_X = pd.DataFrame(index=exp_mtx.index)
embeddings_Y = pd.DataFrame(index=exp_mtx.index)
for idx, (name, df_embedding) in enumerate(embeddings.items()):
    if(len(df_embedding.columns)!=2):
        raise Exception('The embedding should have two columns.')
```
```
embedding_id = idx - 1 # Default embedding must have id == -1 for SCope.
id2name[embedding_id] = name
```
```
embedding = df_embedding.copy()
embedding.columns=['_X', '_Y']
embeddings_X = pd.merge(embeddings_X, embedding['_X'].to_frame().rename(columns={'_X': str(embedding_id)}), left_index=True, right_index=True)
embeddings_Y = pd.merge(embeddings_Y, embedding['_Y'].to_frame().rename(columns={'_Y': str(embedding_id)}), left_index=True, right_index=True)
```
```
# Calculate the number of genes per cell
binary_mtx = exp_mtx.copy()
binary_mtx[binary_mtx != 0] = 1.0
ngenes = binary_mtx.sum(axis=1).astype(int)
```
```
# Encode genes in regulons as a binary membership matrix
from operator import attrgetter
genes = np.array(exp_mtx.columns)
n_genes = len(genes)
n_regulons = len(regulons)
data = np.zeros(shape=(n_genes, n_regulons), dtype=int)
for idx, regulon in enumerate(regulons):
    data[:, idx] = np.isin(genes, regulon.genes).astype(int)
regulon_assignment = pd.DataFrame(data=data,
                                    index=exp_mtx.columns,
                                    columns=list(map(attrgetter('name'), regulons)))
```
```
# Encode cell annotations
name2idx = dict(map(reversed, enumerate(sorted(set(annotations.values())))))

clusterings = pd.DataFrame(data=exp_mtx.index.values,
                               index=exp_mtx.index,
                               columns=['0']).replace(annotations).replace(name2idx)

clusterings.head() # Data for first cell always gets mixed up; add correct cluster manually (below)
clusterings.iloc[[0],[0]] = 34
clusterings.head()
```
```
# Encode cluster_labels (for RSS score calculation)
cluster_labels = pd.DataFrame(data=exp_mtx.index.values,
                               index=exp_mtx.index,
                               columns=['0']).replace(annotations)

cluster_labels.head() # Data for first cell always gets mixed up; add correct cluster manually (below)
cluster_labels.iloc[[0],[0]] = "sAB"
```
Create metadata structure of loom file
```
def create_structure_array(df):
    return np.array([tuple(row) for row in df.values],
                    dtype=np.dtype(list(zip(df.columns, df.dtypes))))
```
```
    default_embedding = next(iter(embeddings.values())).copy()
    default_embedding.columns=['_X', '_Y']
    column_attrs = {
        "CellID": exp_mtx.index.values.astype('str'),
        "nGene": ngenes.values,
        "Embedding": create_structure_array(default_embedding),
        "RegulonsAUC": create_structure_array(auc_mtx),
        "Clusterings": create_structure_array(clusterings),
        "ClusterID": clusterings.values,
        "Cluster_labels":cluster_labels.values,
        'Embeddings_X': create_structure_array(embeddings_X),
        'Embeddings_Y': create_structure_array(embeddings_Y),
        }
    row_attrs = {
        "Gene": exp_mtx.columns.values.astype('str'),
        "Regulons": create_structure_array(regulon_assignment),
        }
```
```
def fetch_logo(context):
    for elem in context:
        if elem.endswith('.png'):
            return elem
    return ""
```
```
def fetch_logo(context):
    for elem in context:
        if elem.endswith('.png'):
            return elem
    return ""
name2logo = {reg.name: fetch_logo(reg.context) for reg in regulons}
regulon_thresholds = [{"regulon": name,
                        "defaultThresholdValue":(threshold if isinstance(threshold, float) else threshold[0]),
                        "defaultThresholdName": "guassian_mixture_split",
                        "allThresholds": {"guassian_mixture_split": (threshold if isinstance(threshold, float) else threshold[0])},
                        "motifData": name2logo.get(name, "")} for name, threshold in auc_mtx.iteritems()] 
```
```
import json
general_attrs = {
    "title": "RegulonsMouseToothAtlas",
    "MetaData": json.dumps({
        "embeddings": [{'id': identifier, 'name': name} for identifier, name in id2name.items()],
        "annotations": [{
            "name": "",
            "values": []
        }],
        "clusterings": [{
            "id": 0,
            "group": "celltype",
            "name": "Cell Type",
            "clusters": [{"id": idx, "description": name} for name, idx in name2idx.items()]
        }],
        "regulonThresholds": regulon_thresholds
    }),
    "Genome": "mm10"}
```
Add tree structure
```
tree_structure: Sequence[str] = ()
```
```
from itertools import islice
import itertools
assert len(tree_structure) <= 3, ""
general_attrs.update(("SCopeTreeL{}".format(idx+1), category)
                        for idx, category in enumerate(list(islice(itertools.chain(tree_structure, itertools.repeat("")), 3))))
```
```
def compress_encode(value):
    '''
    Compress using ZLIB algorithm and encode the given value in base64.
    From: https://github.com/aertslab/SCopeLoomPy/blob/5438da52c4bcf48f483a1cf378b1eaa788adefcb/src/scopeloompy/utils/__init__.py#L7
    '''
    return base64.b64encode(zlib.compress(value.encode('ascii'))).decode('ascii')
```
```
import base64
import zlib
general_attrs["MetaData"] = compress_encode(value=general_attrs["MetaData"])
```
Create loom file:
```
import loompy as lp
lp.create(filename="MTA_Regulons_loom.loom",
              layers=exp_mtx.T.values,
              row_attrs=row_attrs,
              col_attrs=column_attrs,
              file_attrs=general_attrs)
```
```
from pyscenic.genesig import Regulon
from typing import Union
```
```
def add_scenic_metadata(adata: 'sc.AnnData',
                        auc_mtx: pd.DataFrame,
                        regulons: Union[None, Sequence[Regulon]] = None,
                        bin_rep: bool = False,
                        copy: bool = False) -> 'sc.AnnData':
    """
    Add AUCell values and regulon metadata to AnnData object.
    :param adata: The AnnData object.
    :param auc_mtx: The dataframe containing the AUCell values (#observations x #regulons).
    :param bin_rep: Also add binarized version of AUCell values as separate representation. This representation
    is stored as `adata.obsm['X_aucell_bin']`.
    :param copy: Return a copy instead of writing to adata.
    :
    """
    # To avoid dependency with scanpy package the type hinting intentionally uses string literals.
    # In addition, the assert statement to assess runtime type is also commented out.
    #assert isinstance(adata, sc.AnnData)
    assert isinstance(auc_mtx, pd.DataFrame)
    assert len(auc_mtx) == adata.n_obs

    REGULON_SUFFIX_PATTERN = 'Regulon({})'

    result = adata.copy() if copy else adata

    # Add AUCell values as new representation (similar to a PCA). This facilitates the usage of
    # AUCell as initial dimensional reduction.
    result.obsm['X_aucell'] = auc_mtx.values.copy()
    if bin_rep:
        bin_mtx, _ = binarize(auc_mtx)
        result.obsm['X_aucell_bin'] = bin_mtx.values

    # Encode genes in regulons as "binary" membership matrix.
    if regulons is not None:
        genes = np.array(adata.var_names)
        data = np.zeros(shape=(adata.n_vars, len(regulons)), dtype=bool)
        for idx, regulon in enumerate(regulons):
            data[:, idx] = np.isin(genes, regulon.genes).astype(bool)
        regulon_assignment = pd.DataFrame(data=data, index=genes,
                                          columns=list(map(lambda r: REGULON_SUFFIX_PATTERN.format(r.name), regulons1)))
        result.var = pd.merge(result.var, regulon_assignment, left_index=True, right_index=True, how='left')

    # Add additional meta-data/information on the regulons.
    def fetch_logo(context):
        for elem in context:
            if elem.endswith('.png'):
                return elem
        return ""
    result.uns['aucell'] = {
        'regulon_names': auc_mtx.columns.map(lambda s: REGULON_SUFFIX_PATTERN.format(s)).values,
        'regulon_motifs': np.array([fetch_logo(reg.context) for reg in regulons] if regulons is not None else [])
    }

    # Add the AUCell values also as annotations of observations. This way regulon activity can be
    # depicted on cellular scatterplots.
    mtx = auc_mtx.copy()
    mtx.columns = result.uns['aucell']['regulon_names']
    result.obs = pd.merge(result.obs, mtx, left_index=True, right_index=True, how='left')

    return result
```
```
def export_regulons(regulons: Sequence[Regulon], fname: str) -> None:
    """
    Export regulons as GraphML.
    :param regulons: The sequence of regulons to export.
    :param fname: The name of the file to create.
    """
    graph = nx.DiGraph()
    for regulon in regulons:
        src_name = regulon.transcription_factor
        graph.add_node(src_name, group='transcription_factor')
        edge_type = 'activating' if 'activating' in regulon.context else 'inhibiting'
        node_type = 'activated_target' if 'activating' in regulon.context else 'inhibited_target'
        for dst_name, edge_strength in regulon.gene2weight.items():
            graph.add_node(dst_name, group=node_type, **regulon.context)
            graph.add_edge(src_name, dst_name, weight=edge_strength, interaction=edge_type, **regulon.context)
    nx.readwrite.write_graphml(graph, fname)
```


#### (5) Calculate RSS scores
```
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from pyscenic.binarization import binarize
```
```
f_final_loom = 'MTA_Regulons_loom.loom'
```
```
lf = lp.connect( f_final_loom, mode='r', validate=False )
auc_mtx_loom = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
cellAnnot = pd.concat([pd.DataFrame( lf.ca.Cluster_labels, index=lf.ca.CellID )], axis=1)
lf.close()
```
```
rss_cellType = regulon_specificity_scores(auc_mtx, cellAnnot[0])
rss_cellType
```
```
rss_cellType.to_csv("rss_cellType.csv")
```
```
#RSS panel plot with all cell types

cats = sorted(list(set(cellAnnot[0])))

fig = plt.figure(figsize=(40, 40))
for c,num in zip(cats, range(1,len(cats)+1)):
    x=rss_clusters.T[c]
    ax = fig.add_subplot(5,7,num)
    plot_rss(rss_clusters, c, top_n=10, max_n=None, ax=ax)
    ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
 
fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
        'figure.titlesize': 'large' ,
        'axes.labelsize': 'medium',
        'axes.titlesize':'large',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'
        })
plt.savefig("clusters-RSS-top5.pdf", dpi=600, bbox_inches = "tight")
plt.show()
```

#### (6) Calculate MRA (R)
Move back to our R environment, as done before. 
```
integrated
```
```
auc_mtx <- read.csv("/MTA/SCENIC/auc_mtx_corrected_index.csv")
```
```
UMAP_dims <- integrated@reductions$umap_50dims@cell.embeddings
annotation <- integrated@meta.data
```
```
anno <- data.frame(UMAP_dims[,-1],
                   Cell = UMAP_dims$X,
                   clusters = annotation$Annotation,
                   dataset = annotation$Dataset,
                   age = annotation$Age,
                   tooth_type = annotation$Tooth_type,
                  row.names = UMAP_dims$X)
```
```
# Ensure cell IDs match
anno$Cell <- gsub('\\.', '-', anno$Cell) 
rownames(anno) <- anno$Cell
auc_mtx$Cell <- gsub('\\.', '-', auc_mtx$Cell)
UMAP_dims$X <- gsub('\\.', '-', UMAP_dims$X )
rownames(annotation) <- gsub('\\.', '-', rownames(annotation))
```
```
# Join the data
regulon_anno <- inner_join(auc_mtx, anno, by = "Cell")
regulon_anno_long <- gather(regulon_anno, regulon, activity, -clusters, -dataset, -age, -tooth_type, -Cell, -UMAP50_1, -UMAP50_2)
```
```
# Remove 3 dots at the end of regulon's name
regulon_anno_long$regulon <- gsub('.{0,3}$', '', regulon_anno_long$regulon)
```
Calculate MRA:
```
meanRegPerCluster <- regulon_anno_long %>%
                        dplyr::select(-c(Cell, UMAP50_1, UMAP50_2)) %>%
                        group_by(clusters, regulon) %>%
                        summarize(mean_activity = mean(activity))
write.table(meanRegPerCluster, "./MTA_SCENIC_MRA_Clusters.txt", col.names=NA, sep="\t")
```

#### (7) Create a regulon-based Seurat object and perform regulon-based integration (R)
```
auc_mtx <- auc_mtx %>%
                dplyr::filter(Cell %in% regulon_anno$Cell)
```
```
# Remove 3 dots at the end of regulon's name
colnames(auc_mtx) <- gsub('\\...','', colnames(auc_mtx))
```
```
auc_mtx <- data.frame(auc_mtx[,-1], row.names = auc_mtx[,1])
```
```
# Transpose matrix
auc_mtx_T <- t(auc_mtx)
```
```
meta_data <- data.frame(annotation)
```
```
seurat <- CreateSeuratObject(counts = auc_mtx_T, meta.data = meta_data, min.cells = 0, min.features = 0, project = "AUC")
```
Perform standard integration
```
tooth.list <- SplitObject(seurat, split.by = "Dataset")
tooth.list <- lapply(X = tooth.list, FUN = function(x) {
    x <- FindVariableFeatures(x, verbose = FALSE)
})
```
```
tooth.anchors <- FindIntegrationAnchors(object.list = tooth.list, dims = 1:30)
integrated_regulons <- IntegrateData(anchorset = tooth.anchors, dims = 1:30)
```
```
DefaultAssay(integrated_regulons) <- "integrated_regulons"
integrated_regulons <- ScaleData(integrated_regulons, verbose = FALSE, vars.to.regress = c("S.Score", "G2M.Score"))
integrated_regulons <- RunPCA(integrated_regulons, npcs = 100, verbose = FALSE)
```
```
integrated_regulons <- RunUMAP(integrated_regulons, umap.method = "umap-learn", dims = 1:15, min.dist = 0.5)
umap_dims <- integrated_regulons@reductions$umap@cell.embeddings
integrated_regulons[["umap_15dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP15_", assay = DefaultAssay(integrated_regulons))
```


## Inference of LR interactions using CellPhoneDB
First we divide the established MTA into our three groups (incisors, molars, periodontal). Then we convert the mouse gene names to human gene names (required for CellphoneDB) using the Homologene table from [Hart, R., BioRxiv (2019)](https://www.biorxiv.org/content/10.1101/671115v1)). Next, we extract the required data to run CellPhoneDB. Finally, we move to our Python Conda environment to run CellPhoneDB. 

```
# Split object into 3 groups
split <- SplitObject(integrated, split.by = "Group")

incisors <- split$Incisors
molars <- split$Molar
periodontal <- split$Periodontal
```

```
# Load in Homologene table
geneTrans=read.table("./geneTrans.txt",sep=",",header=T,stringsAsFactors = F,row.names = 1)

geneTrans$hg19=gsub("hg19_","",geneTrans$hg19)
geneTrans$mm10=gsub("mm10_","",geneTrans$mm10)
row.names(geneTrans)=gsub("_","-",row.names(geneTrans))
```
Below is the code for the incisors. Replace with the other Seurat objects for the other groups. 
```
#extract counts
mm.raw=GetAssayData(incisors,slot="counts") # replace with molars/periodontal to perform on the other groups.

#subset rows in geneTrans (not necessary but simpler)
mm.raw=mm.raw[row.names(mm.raw) %in% geneTrans$mm10,]

#translate mouse symbols to human
mm.trans=merge(x=mm.raw,y=geneTrans[,c(5,6)],by.x=0,by.y="mm10",all.x=T)
rownames(mm.trans)=mm.trans$hg19
mm.trans=mm.trans[,!(names(mm.trans) %in% c("Row.names","hg19"))]
```

```
# Save converted gene names + normalized counts
#write.table(mm.trans, "./incisor_cpdb_count.txt", sep="\t", quote=F)

# generating meta file
meta_data <- cbind(rownames(incisors@meta.data), incisors@meta.data[,"Annotation", drop=F])
#write.table(meta_data, "./incisor_cpdb_meta.txt", sep="\t", quote=F, row.names=F)
```

For the next part of the CellPhoneDB analysis of Ligand-Receptor interactions in the MTA we use our Python Conda environment with the CellPhoneDB tool installed, and perform the analysis using the command line. 

``` 
conda activate py_cpdb
```

```
# Example for incisor group, same done for molar/periodontal groups
cellphonedb method statistical_analysis 
  --counts-data=gene_name incisor_cpdb_meta.txt incisor_cpdb_count.txt 
  --iterations=1000 
  --threshold=0.2 
  --threads=8 
  --output-path=CellPhoneDB 
  --project-name=MTA_incisor 
  --quiet
```

Next, make a rows.txt and columns.txt file with the L/R pairs and cluster-cluster pairs for further analysis. Run the CellPhoneDB plotting function for each group.

```
# Example for incisor group, same done for molar/periodontal groups
cellphonedb plot dot_plot 
  --means-path=./CellPhoneDB/MTA_incisor/means.txt 
  --pvalues-path=./CellPhoneDB/MTA_incisor/pvalues.txt 
  --rows ./CellPhoneDB/MTA_incisor/rows.txt 
  --columns ./CellPhoneDB/MTA_incisor/columns.txt 
  --output-path=./CellPhoneDB/MTA_incisor/ 
  --output-name=cpdb_MTA_incisor_dotplot.pdf
```

## Subclustering of Mouse Epithelium
Here we focus on the ameloblast lineage trajectory, starting from the IEE/OEE cluster. 

### Setup
```
mouse_epithelial <- subset(x = mta, idents = c("sAB", "mAB", "preAB", "DEP", "Cycling DEP", "IEE/OEE"))

DefaultAssay(mouse_epithelial) <- "SoupX"

mouse_epithelial[['integrated']] <- NULL
mouse_epithelial[['RNA']] <- NULL
mouse_epithelial <- RenameAssays(object = mouse_epithelial, SoupX = 'RNA')
```

We divide the datasets into 3 groups (if we would integrate on the datasets, as we did before, it wouldn't work due to too few cells for some datasets). By dividing the datasets into three groups (where we tried to group them based on biological similarity) we circumvent this problem and are able to integrate the data. 

```
# Divide epithelial data into 3 groups to enable integration
Idents(mouse_epithelial) <- mouse_epithelial@meta.data$Dataset

mouse_epithelial <- RenameIdents(mouse_epithelial, `Krivanek_molar` = "Group3", `Zhao_perio_1` = "Group3", 
                    `Zhao_perio_2` = "Group3", `Krivanek_incisor` = "Group3", `Nagata` = "Group3",
                    `Chen` = "Group4", `Chiba` = "Group5", `Chiba_2_epithelium` = "Group6", `Chiba_2_mesenchyme` = "Group7", 
                    `Sharir` = "Group8", `Takahashi` = "Group3", `Wen` = "Group3")

mouse_epithelial@meta.data$Integration_Groups <- Idents(mouse_epithelial)
```

### Integration
```
# As done before using rPCA
tooth.list <- SplitObject(mouse_epithelial, split.by = "Integration_Groups")
tooth.list <- lapply(X = tooth.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = tooth.list)
tooth.list <- lapply(X = tooth.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE) # approx = false https://github.com/satijalab/seurat/issues/1963
})

tooth.anchors <- FindIntegrationAnchors(object.list = tooth.list, dims = 1:30, anchor.features = features, reduction = "rpca")
epi_integrated <- IntegrateData(anchorset = tooth.anchors, dims = 1:30)

DefaultAssay(epi_integrated) <- "epi_integrated"

epi_integrated <- ScaleData(epi_integrated, verbose = FALSE)
epi_integrated <- RunPCA(epi_integrated, verbose = FALSE)
```

```
# Run UMAP for different dimensions, save to object separately to be able to compare

epi_integrated <- RunUMAP(epi_integrated, umap.method = "umap-learn", dims = 1:10)
umap_dims <- epi_integrated@reductions$umap@cell.embeddings
epi_integrated[["umap_10dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP10_", assay = DefaultAssay(epi_integrated))

epi_integrated <- RunUMAP(epi_integrated, umap.method = "umap-learn", dims = 1:20)
umap_dims <- epi_integrated@reductions$umap@cell.embeddings
epi_integrated[["umap_20dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP20_", assay = DefaultAssay(epi_integrated))

epi_integrated <- RunUMAP(epi_integrated, umap.method = "umap-learn", dims = 1:30)
umap_dims <- epi_integrated@reductions$umap@cell.embeddings
epi_integrated[["umap_30dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP30_", assay = DefaultAssay(epi_integrated))

epi_integrated <- RunUMAP(epi_integrated, umap.method = "umap-learn", dims = 1:40)
umap_dims <- epi_integrated@reductions$umap@cell.embeddings
epi_integrated[["umap_40dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP40_", assay = DefaultAssay(epi_integrated))

epi_integrated <- RunUMAP(epi_integrated, umap.method = "umap-learn", dims = 1:50)
umap_dims <- epi_integrated@reductions$umap@cell.embeddings
epi_integrated[["umap_50dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP50_", assay = DefaultAssay(epi_integrated))
```

By plotting the metadata, prior cluster annotation and marker genes we selected umap_40dims for further downstream analyses.

### Redo Nebulosa and Module Score analysis on subclustered epithelium
```
DefaultAssay(epi_integrated) <- "SoupX"
options(repr.plot.width=12, repr.plot.height=12)
```

```
# Gene expression of DESC genes
DESC_features <- list(c('Sox2','Lgr5','Gli1','Lrig1','Bmi1','Ptch1','Pknox2','Zfp273','Spock1','Pcp4','Sfrp5'))

plot.list <- list()

for (i in 1:length(DESC_features)) {
    options(repr.plot.width=12, repr.plot.height=12)
    plot.list[[i]] <- FeaturePlot(epi_integrated, reduction = "umap_40dims", features = DESC_features[i], pt.size = 1, cols = (c("lightgrey", "steelblue", "blue4")))
    print(plot.list[[i]])
}
```

```
# DESC Module with AddModuleScore (Seurat)

DESC_features <- list(c('Sox2','Lgr5','Gli1','Lrig1','Bmi1','Ptch1','Pknox2','Zfp273','Spock1','Pcp4','Sfrp5'))

epi_integrated <- AddModuleScore(
  object = epi_integrated,
  features = DESC_features,
  ctrl = 5,
  name = 'DESC_features',
    assay = "SoupX"
)

FeaturePlot(epi_integrated, reduction = "umap_40dims", features = "DESC_features1", pt.size = 1, cols = (c("lightgrey", "steelblue", "blue4")))
```

```
# DESC (joint) densities with Nebulosa
p_list <- plot_density(epi_integrated, c('Sox2','Lgr5','Gli1','Lrig1','Bmi1','Ptch1','Pknox2','Zfp273','Spock1','Pcp4','Sfrp5'),
  joint = TRUE, combine = FALSE)
  
p_list[[length(p_list)]]
```

### Pseudotime analysis on subclustered epithelium using Monocle3
First we manually create the cds object required for Monocle3, from our Seurat object

```
### Building the necessary parts for a basic cds
# part one, gene annotations
gene_annotation <- as.data.frame(rownames(epi_integrated@reductions[["pca"]]@feature.loadings), row.names = rownames(epi_integrated@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"
# part two, cell information
cell_metadata <- as.data.frame(epi_integrated@assays[["RNA"]]@counts@Dimnames[[2]], row.names = epi_integrated@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
# part three, counts sparse matrix
New_matrix <- epi_integrated@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(epi_integrated@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix
#Cluster information
list_cluster <- epi_integrated@meta.data[[sprintf("Annotation")]]
names(list_cluster) <- epi_integrated@assays[["RNA"]]@data@Dimnames[[2]]
#UMAP coordinates
umap <- epi_integrated@reductions[["umap_40dims"]]@cell.embeddings
#Feature loadings
feature_loadings <- epi_integrated@reductions[["pca"]]@feature.loadings


# Make the CDS object
cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

### Construct and assign the made up partition

recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

### Assign the cluster info
cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

### Could be a space-holder, but essentially fills out louvain parameters
cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

# Assign UMAP coordinates
cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <- umap

#Assign feature loading for downstream module analysis
cds_from_seurat@preprocess_aux$gene_loadings <- feature_loadings

#Learn graph
cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = T)
```

```
colData(cds_from_seurat)$new_clusters <- as.character(cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]])

# Select a cluster and collect unique vertex information for cells in that cluster to start pseudotime manually
root_group = colnames(cds_from_seurat)[pData(cds_from_seurat)$new_clusters == "IEE/OEE"] # We start from IEE/OEE, the putative stem cell region
closest_vertex = cds_from_seurat@principal_graph_aux[['UMAP']]$pr_graph_cell_proj_closest_vertex
root_pr_nodes = paste("Y_", closest_vertex[root_group,], sep="")
unique_nodes = unique(root_pr_nodes)
unique_nodes
unique_nodes <- as.vector(unique_nodes) # this manually prints a list of nodes, by trial and error we identify the node at the root of our trajectory
```

```
# Input the manually identified root node to order cells and calculate pseudotime
cds_from_seurat = order_cells(cds_from_seurat, root_pr_nodes = c('Y_52'))
```

```
# Plot results
options(repr.plot.width=20, repr.plot.height=7)
p1 <- plot_cells(cds_from_seurat, 
                color_cells_by = 'cluster',
                label_groups_by_cluster=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                cell_size = 1.25,
                 label_roots = FALSE,
                group_label_size = 5)
p2 <- plot_cells(cds_from_seurat,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           cell_size = 1.25,
           graph_label_size=5)
p1+p2
```

```
# Extract pseudotime values and add to Seurat object
epi_integrated <- AddMetaData(
  object = epi_integrated,
  metadata = cds_from_seurat@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Epi_pseudotime"
)
```

Next, we calculate DEGs along pseudotime using the method of [Jevitt, A. et al., PLOS Biology (2020) ](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000538), as described on their [GitHub page](https://github.com/crickbabs/ZebrafishDevelopingHindbrainAtlas#trajectory-analysis)
```
#Calculate DEGs along pseudotime
pr_graph_test_res <- graph_test(cds_from_seurat, neighbor_graph="principal_graph", cores=1, verbose = FALSE) # Group relative to pseudotime
hits      <- as.data.frame(pr_graph_test_res[pr_graph_test_res$gene_short_name,])
hits$pass <- hits$morans_I > 0.1 & hits$q_value < 0.01

write.table(hits,file="./genes_changing_with_pseudotime.txt",col.names=T,row.names=F,sep="\t",quote=F)
genes <- read.delim("./changing_with_pseudotime.txt",header=T,sep="\t",stringsAsFactors=F)[,5] #extract list of genes identified
```

```
# Intersect identified genes with matrix, z-score normalization of expression
pt.matrix <- as.matrix(cds_from_seurat@assays@data$counts[match(genes,rowData(cds_from_seurat)[,1]),order(pseudotime(cds_from_seurat))])
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
save(pt.matrix, file="matrix.Rdata")
```


We have a separate Conda environment containing the ComplexHeatmap package (R_heatmap xxx)
```
#Setup environment
#set working directory
setwd("./MTA")

#load in required software packages
suppressPackageStartupMessages({
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(magick)
})

load("matrix.Rdata")
```

```
# Create pseudotime heatmap
options(repr.plot.width=7, repr.plot.height=7)
ht <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  split                        = split,
  row_title_rot                = 0,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE
  raster_by_magick             = TRUE)
print(ht)
```

# Human Tooth Atlas

## Healthy HTA setup


```
suppressMessages({
#set working directory
setwd("./HTA")
getwd()

#load in required software packages
suppressWarnings({
library(reticulate)
library(dplyr)
library(Seurat)
library(SingleCellExperiment, quietly = TRUE)
library(scater, quietly = TRUE)
library(tidyr)
library(purrr)
library(cowplot)
library(ggrepel)
library(viridis)
library(gridExtra)
library(biomaRt)
library(schex)
library(Nebulosa)
library(dittoSeq)
library(escape)
library(SoupX)
library(monocle3) 
})

# choose random seed for reproducibility
set.seed(13101999)

})
```

Import all datasets, and add relevant metadata obtained from the publications. 

```
#Import data
GSM4998457 <- Read10X(data.dir = "./GSM4998457_Pagella_Pulp_1")
Pagella_Pulp_1 <- CreateSeuratObject(counts = GSM4998457, project = "GSM4998457", min.cells = 3, min.features = 200)

GSM4998458 <- Read10X(data.dir = "./GSM4998458_Pagella_Pulp_2")
Pagella_Pulp_2 <- CreateSeuratObject(counts = GSM4998458, project = "GSM4998458", min.cells = 3, min.features = 200)

GSM4998459 <- Read10X(data.dir = "./GSM4998459_Pagella_Pulp_3")
Pagella_Pulp_3 <- CreateSeuratObject(counts = GSM4998459, project = "GSM4998459", min.cells = 3, min.features = 200)

GSM4998460 <- Read10X(data.dir = "./GSM4998460_Pagella_Pulp_4")
Pagella_Pulp_4 <- CreateSeuratObject(counts = GSM4998460, project = "GSM4998460", min.cells = 3, min.features = 200)

GSM4998461 <- Read10X(data.dir = "./GSM4998461_Pagella_Pulp_5")
Pagella_Pulp_5 <- CreateSeuratObject(counts = GSM4998461, project = "GSM4998461", min.cells = 3, min.features = 200)

GSM4904134 <- Read10X(data.dir = "./GSM4904134_Pagella_Perio_1")
Pagella_Perio_1 <- CreateSeuratObject(counts = GSM4904134, project = "GSM4904134", min.cells = 3, min.features = 200)

GSM4904135 <- Read10X(data.dir = "./GSM4904135_Pagella_Perio_2")
Pagella_Perio_2 <- CreateSeuratObject(counts = GSM4904135, project = "GSM4904135", min.cells = 3, min.features = 200)

GSM4904136 <- Read10X(data.dir = "./GSM4904136_Pagella_Perio_3")
Pagella_Perio_3 <- CreateSeuratObject(counts = GSM4904136, project = "GSM4904136", min.cells = 3, min.features = 200)

GSM4904137 <- Read10X(data.dir = "./GSM4904137_Pagella_Perio_4")
Pagella_Perio_4 <- CreateSeuratObject(counts = GSM4904137, project = "GSM4904137", min.cells = 3, min.features = 200)

GSM4904138 <- Read10X(data.dir = "./GSM4904138_Pagella_Perio_5")
Pagella_Perio_5 <- CreateSeuratObject(counts = GSM4904138, project = "GSM4904138", min.cells = 3, min.features = 200)

xxx Hemeryck/Yin/Opasawatchai
```

```
# Import of Krivanek apical papilla 1 (human germinectomy)
data = read.table(gzfile("./GSM4365601_counts_human_germinectomy.txt.gz"), sep=" ", header = TRUE, row.names = NULL)

names <- make.unique(data$row.names)
rownames(data) <- names
data <- data[,-1] # get rid of old names

#Fix cell names
colnames(data) <- gsub('\\.', '-', colnames(data))
colnames(data) <- gsub('human_germinectomy_1_one', 'Krivanek_Apical_Papilla_1', colnames(data))

#Create Seurat Object
Krivanek_Apical_Papilla_1 <- CreateSeuratObject(counts = data, project = "GSM4365601", min.cells = 3, min.features = 200)
```

```
# Import of Krivanek apical papilla 2 (apical_papilla_female_15yo)
data = read.table(gzfile("./GSM4365609_counts_human_germ_molar_apical_papilla_female_15yo.txt.gz"), sep=" ", header = TRUE, row.names = NULL)

names <- make.unique(data$row.names)
rownames(data) <- names
data <- data[,-1] # get rid of old names

#Fix cell names
colnames(data) <- gsub('\\.', '-', colnames(data))
colnames(data) <- gsub('human_germ_molar_apical_papilla_female_15yo_one', 'Krivanek_Apical_Papilla_2', colnames(data))

#Create Seurat Object
Krivanek_Apical_Papilla_2 <- CreateSeuratObject(counts = data, project = "GSM4365609", min.cells = 3, min.features = 200)
```

```
# Import of Krivanek dental pulp (human_molar_pulp)
data = read.table(gzfile("./GSM4365602_counts_human_molar_pulp.txt.gz"), sep=" ", header = TRUE, row.names = NULL)

names <- make.unique(data$row.names)
rownames(data) <- names
data <- data[,-1] # get rid of old names

#Fix cell names
colnames(data) <- gsub('\\.', '-', colnames(data))
colnames(data) <- gsub('human_molar_pulp_4_one', 'Krivanek_Dental_Pulp', colnames(data))

#Create Seurat Object
Krivanek_Dental_Pulp <- CreateSeuratObject(counts = data3, project = "GSM4365602", min.cells = 3, min.features = 200)
```

```
# Import of Krivanek whole molar 1 (molar_healthy_1)
data = read.table(gzfile("./GSM4365608_counts_human_molar_healthy_1.txt.gz"), sep=" ", header = TRUE, row.names = NULL)

names <- make.unique(data$row.names)
rownames(data) <- names
data <- data[,-1] # get rid of old names

#Fix cell names
colnames(data) <- gsub('\\.', '-', colnames(data))
colnames(data) <- gsub('human_molar_healthy_1_one', 'Krivanek_Whole_Molar_1', colnames(data))

#Create Seurat Object
Krivanek_Whole_Molar_1 <- CreateSeuratObject(counts = data4, project = "GSM4365608", min.cells = 3, min.features = 200)
```

```
# Import of Krivanek whole molar 2 (molar_healthy_2)
data = read.table(gzfile("./GSM4365607_counts_human_molar_healthy_2.txt.gz"), sep=" ", header = TRUE, row.names = NULL)

names <- make.unique(data$row.names)
rownames(data) <- names
data <- data[,-1] # get rid of old names

#Fix cell names
colnames(data) <- gsub('\\.', '-', colnames(data))
colnames(data) <- gsub('human_molar_healthy_2_one', 'Krivanek_Whole_Molar_2', colnames(data))

#Create Seurat Object
Krivanek_Whole_Molar_2 <- CreateSeuratObject(counts = data5, project = "GSM4365607", min.cells = 3, min.features = 200)
```

```
# Import of Krivanek whole molar 3 (No5_24_yo_healthy_retained)
data = read.table(gzfile("./GSM4365610_counts_No5_24_yo_healthy_retained.txt.gz"), sep=" ", header = TRUE, row.names = NULL)

names <- make.unique(data$row.names)
rownames(data) <- names
data <- data[,-1] # get rid of old names

#Fix cell names
colnames(data) <- gsub('\\.', '-', colnames(data))
colnames(data) <- gsub('No5_24_yo_healthy_retained_one', 'Krivanek_Whole_Molar_3', colnames(data))

#Create Seurat Object
Krivanek_Whole_Molar_3 <- CreateSeuratObject(counts = data6, project = "GSM4365610", min.cells = 3, min.features = 200)
```

```
# Import of Shi tooth germ dataset (from Mendeley Data, Robject)
load("/staging/leuven/stg_00073/SC_Florian/Human/Song_human_tooth_germ.RData")
Song_Tooth_Germ <- all
```

```
# Calculate % of mitochondrial genes per cell, and append this to the metadata
Pagella_Pulp_1[["percent.mito"]] <- PercentageFeatureSet(Pagella_Pulp_1, pattern = "^MT-")
Pagella_Pulp_2[["percent.mito"]] <- PercentageFeatureSet(Pagella_Pulp_2, pattern = "^MT-")
Pagella_Pulp_3[["percent.mito"]] <- PercentageFeatureSet(Pagella_Pulp_3, pattern = "^MT-")
Pagella_Pulp_4[["percent.mito"]] <- PercentageFeatureSet(Pagella_Pulp_4, pattern = "^MT-")
Pagella_Pulp_5[["percent.mito"]] <- PercentageFeatureSet(Pagella_Pulp_5, pattern = "^MT-")
Pagella_Perio_1[["percent.mito"]] <- PercentageFeatureSet(Pagella_Perio_1, pattern = "^MT-")
Pagella_Perio_2[["percent.mito"]] <- PercentageFeatureSet(Pagella_Perio_2, pattern = "^MT-")
Pagella_Perio_3[["percent.mito"]] <- PercentageFeatureSet(Pagella_Perio_3, pattern = "^MT-")
Pagella_Perio_4[["percent.mito"]] <- PercentageFeatureSet(Pagella_Perio_4, pattern = "^MT-")
Pagella_Perio_5[["percent.mito"]] <- PercentageFeatureSet(Pagella_Perio_5, pattern = "^MT-")
Krivanek_Apical_Papilla_1[["percent.mito"]] <- PercentageFeatureSet(Krivanek_Apical_Papilla_1, pattern = "^MT-")
Krivanek_Apical_Papilla_2[["percent.mito"]] <- PercentageFeatureSet(Krivanek_Apical_Papilla_2, pattern = "^MT-")
Krivanek_Dental_Pulp[["percent.mito"]] <- PercentageFeatureSet(Krivanek_Dental_Pulp, pattern = "^MT-")
Krivanek_Whole_Molar_1[["percent.mito"]] <- PercentageFeatureSet(Krivanek_Whole_Molar_1, pattern = "^MT-")
Krivanek_Whole_Molar_2[["percent.mito"]] <- PercentageFeatureSet(Krivanek_Whole_Molar_2, pattern = "^MT-")
Krivanek_Whole_Molar_3[["percent.mito"]] <- PercentageFeatureSet(Krivanek_Whole_Molar_3, pattern = "^MT-")
#Song_Tooth_Germ dataset already contains percent.mito metadata
```

```
# Add extra metadata to Seurat object
Pagella_Pulp_1@meta.data$Dataset <- "Pagella Pulp 1"
Pagella_Pulp_1@meta.data$Tissue <- "Dental Pulp"
Pagella_Pulp_1@meta.data$Extraction <- "Third Molar"
Pagella_Pulp_1@meta.data$Species <- "Human"
Pagella_Pulp_1@meta.data$Age <- "18-35"
Pagella_Pulp_1@meta.data$Conditon <- "Healthy"
Pagella_Pulp_1@meta.data$Technology <- "10X"

Pagella_Pulp_2@meta.data$Dataset <- "Pagella Pulp 2"
Pagella_Pulp_2@meta.data$Tissue <- "Dental Pulp"
Pagella_Pulp_2@meta.data$Extraction <- "Third Molar"
Pagella_Pulp_2@meta.data$Species <- "Human"
Pagella_Pulp_2@meta.data$Age <- "18-35"
Pagella_Pulp_2@meta.data$Conditon <- "Healthy"
Pagella_Pulp_2@meta.data$Technology <- "10X"

Pagella_Pulp_3@meta.data$Dataset <- "Pagella Pulp 3"
Pagella_Pulp_3@meta.data$Tissue <- "Dental Pulp"
Pagella_Pulp_3@meta.data$Extraction <- "Third Molar"
Pagella_Pulp_3@meta.data$Species <- "Human"
Pagella_Pulp_3@meta.data$Age <- "18-35"
Pagella_Pulp_3@meta.data$Conditon <- "Healthy"
Pagella_Pulp_3@meta.data$Technology <- "10X"

Pagella_Pulp_4@meta.data$Dataset <- "Pagella Pulp 4"
Pagella_Pulp_4@meta.data$Tissue <- "Dental Pulp"
Pagella_Pulp_4@meta.data$Extraction <- "Third Molar"
Pagella_Pulp_4@meta.data$Species <- "Human"
Pagella_Pulp_4@meta.data$Age <- "18-35"
Pagella_Pulp_4@meta.data$Conditon <- "Healthy"
Pagella_Pulp_4@meta.data$Technology <- "10X"

Pagella_Pulp_5@meta.data$Dataset <- "Pagella Pulp 5"
Pagella_Pulp_5@meta.data$Tissue <- "Dental Pulp"
Pagella_Pulp_5@meta.data$Extraction <- "Third Molar"
Pagella_Pulp_5@meta.data$Species <- "Human"
Pagella_Pulp_5@meta.data$Age <- "18-35"
Pagella_Pulp_5@meta.data$Conditon <- "Healthy"
Pagella_Pulp_5@meta.data$Technology <- "10X"

Pagella_Perio_1@meta.data$Dataset <- "Pagella Perio 1"
Pagella_Perio_1@meta.data$Tissue <- "Periodontium"
Pagella_Perio_1@meta.data$Extraction <- "Third Molar"
Pagella_Perio_1@meta.data$Species <- "Human"
Pagella_Perio_1@meta.data$Age <- "18-35"
Pagella_Perio_1@meta.data$Conditon <- "Healthy"
Pagella_Perio_1@meta.data$Technology <- "10X"

Pagella_Perio_1@meta.data$Dataset <- "Pagella Perio 1"
Pagella_Perio_1@meta.data$Tissue <- "Periodontium"
Pagella_Perio_1@meta.data$Extraction <- "Third Molar"
Pagella_Perio_1@meta.data$Species <- "Human"
Pagella_Perio_1@meta.data$Age <- "18-35"
Pagella_Perio_1@meta.data$Conditon <- "Healthy"
Pagella_Perio_1@meta.data$Technology <- "10X"

Pagella_Perio_3@meta.data$Dataset <- "Pagella Perio 3"
Pagella_Perio_3@meta.data$Tissue <- "Periodontium"
Pagella_Perio_3@meta.data$Extraction <- "Third Molar"
Pagella_Perio_3@meta.data$Species <- "Human"
Pagella_Perio_3@meta.data$Age <- "18-35"
Pagella_Perio_3@meta.data$Conditon <- "Healthy"
Pagella_Perio_3@meta.data$Technology <- "10X"

Pagella_Perio_4@meta.data$Dataset <- "Pagella Perio 4"
Pagella_Perio_4@meta.data$Tissue <- "Periodontium"
Pagella_Perio_4@meta.data$Extraction <- "Third Molar"
Pagella_Perio_4@meta.data$Species <- "Human"
Pagella_Perio_4@meta.data$Age <- "18-35"
Pagella_Perio_4@meta.data$Conditon <- "Healthy"
Pagella_Perio_4@meta.data$Technology <- "10X"

Pagella_Perio_5@meta.data$Dataset <- "Pagella Perio 5"
Pagella_Perio_5@meta.data$Tissue <- "Periodontium"
Pagella_Perio_5@meta.data$Extraction <- "Third Molar"
Pagella_Perio_5@meta.data$Species <- "Human"
Pagella_Perio_5@meta.data$Age <- "18-35"
Pagella_Perio_5@meta.data$Conditon <- "Healthy"
Pagella_Perio_5@meta.data$Technology <- "10X"

Krivanek_Apical_Papilla_1@meta.data$orig.ident <- NULL
Krivanek_Apical_Papilla_1@meta.data$orig.ident <- "GSM4365601" #metadata not correct --> adapted 
Krivanek_Apical_Papilla_1@meta.data$Dataset <- "Apical_Papilla_1"
Krivanek_Apical_Papilla_1@meta.data$Tissue <- "Apical Papilla"
Krivanek_Apical_Papilla_1@meta.data$Extraction <- "Wisdom Tooth"
Krivanek_Apical_Papilla_1@meta.data$Species <- "Human"
Krivanek_Apical_Papilla_1@meta.data$Age <- "18-31"
Krivanek_Apical_Papilla_1@meta.data$Conditon <- "Healthy"
Krivanek_Apical_Papilla_1@meta.data$Technology <- "10X"

Krivanek_Apical_Papilla_2@meta.data$orig.ident <- NULL
Krivanek_Apical_Papilla_2@meta.data$orig.ident <- "GSM4365609" #metadata not correct --> adapted 
Krivanek_Apical_Papilla_2@meta.data$Dataset <- "Apical_Papilla_2"
Krivanek_Apical_Papilla_2@meta.data$Tissue <- "Apical Papilla"
Krivanek_Apical_Papilla_2@meta.data$Extraction <- "Wisdom Tooth"
Krivanek_Apical_Papilla_2@meta.data$Species <- "Human"
Krivanek_Apical_Papilla_2@meta.data$Age <- "15"
Krivanek_Apical_Papilla_2@meta.data$Conditon <- "Healthy"
Krivanek_Apical_Papilla_2@meta.data$Technology <- "10X"

Krivanek_Dental_Pulp@meta.data$orig.ident <- NULL
Krivanek_Dental_Pulp@meta.data$orig.ident <- "GSM4365602" #metadata not correct --> adapted
Krivanek_Dental_Pulp@meta.data$Dataset <- "Dental_Pulp"
Krivanek_Dental_Pulp@meta.data$Tissue <- "Dental Pulp"
Krivanek_Dental_Pulp@meta.data$Extraction <- "Wisdom Tooth"
Krivanek_Dental_Pulp@meta.data$Species <- "Human"
Krivanek_Dental_Pulp@meta.data$Age <- "18-31"
Krivanek_Dental_Pulp@meta.data$Conditon <- "Healthy"
Krivanek_Dental_Pulp@meta.data$Technology <- "10X"

Krivanek_Whole_Molar_1@meta.data$orig.ident <- NULL
Krivanek_Whole_Molar_1@meta.data$orig.ident <- "GSM4365608" #metadata not correct --> adapted
Krivanek_Whole_Molar_1@meta.data$Dataset <- "Whole_Molar_1"
Krivanek_Whole_Molar_1@meta.data$Tissue <- "Whole Molar"
Krivanek_Whole_Molar_1@meta.data$Extraction <- "Wisdom Tooth"
Krivanek_Whole_Molar_1@meta.data$Species <- "Human"
Krivanek_Whole_Molar_1@meta.data$Age <- "18-31"
Krivanek_Whole_Molar_1@meta.data$Conditon <- "Healthy"
Krivanek_Whole_Molar_1@meta.data$Technology <- "10X"

Krivanek_Whole_Molar_2@meta.data$orig.ident <- NULL
Krivanek_Whole_Molar_2@meta.data$orig.ident <- "GSM4365607" #metadata not correct --> adapted
Krivanek_Whole_Molar_2@meta.data$Dataset <- "Whole_Molar_2"
Krivanek_Whole_Molar_2@meta.data$Tissue <- "Whole Molar"
Krivanek_Whole_Molar_2@meta.data$Extraction <- "Wisdom Tooth"
Krivanek_Whole_Molar_2@meta.data$Species <- "Human"
Krivanek_Whole_Molar_2@meta.data$Age <- "18-31"
Krivanek_Whole_Molar_2@meta.data$Conditon <- "Healthy"
Krivanek_Whole_Molar_2@meta.data$Technology <- "10X"

Krivanek_Whole_Molar_3@meta.data$orig.ident <- NULL
Krivanek_Whole_Molar_3@meta.data$orig.ident <- "GSM4365610" #metadata not correct --> adapted
Krivanek_Whole_Molar_3@meta.data$Dataset <- "Whole_Molar_3"
Krivanek_Whole_Molar_3@meta.data$Tissue <- "Whole Molar"
Krivanek_Whole_Molar_3@meta.data$Extraction <- "Wisdom Tooth"
Krivanek_Whole_Molar_3@meta.data$Species <- "Human"
Krivanek_Whole_Molar_3@meta.data$Age <- "24"
Krivanek_Whole_Molar_3@meta.data$Conditon <- "Healthy"
Krivanek_Whole_Molar_3@meta.data$Technology <- "10X"

Song_Tooth_Germ@meta.data$S.Score <- NULL
Song_Tooth_Germ@meta.data$G2M.Score <- NULL
Song_Tooth_Germ@meta.data$Phase <- NULL
Song_Tooth_Germ@meta.data$CC.Difference <- NULL
Song_Tooth_Germ@meta.data$RNA_snn_res.0.8 <- NULL
Song_Tooth_Germ@meta.data$Dataset <- "Human_Tooth_Germ"
Song_Tooth_Germ@meta.data$Tissue <- "Tooth Germ"
Song_Tooth_Germ@meta.data$Extraction <- "Third molars"
Song_Tooth_Germ@meta.data$Species <- "Human"
Song_Tooth_Germ@meta.data$Age <- "Unkown"
Song_Tooth_Germ@meta.data$Conditon <- "Healthy"
Song_Tooth_Germ@meta.data$Technology <- "BD Rhapsody"

xxx Hemeryck/Yin/Opasawatchai 
```

```
# Add Dataset of origin to cell name to avoid possible identical cell names between datasets
Pagella_Pulp_1 <- RenameCells(Pagella_Pulp_1, new.names = paste0(Pagella_Pulp_1$Dataset, "_", Cells(Pagella_Pulp_1)))
Pagella_Pulp_2 <- RenameCells(Pagella_Pulp_2, new.names = paste0(Pagella_Pulp_2$Dataset, "_", Cells(Pagella_Pulp_2)))
Pagella_Pulp_3 <- RenameCells(Pagella_Pulp_3, new.names = paste0(Pagella_Pulp_3$Dataset, "_", Cells(Pagella_Pulp_3)))
Pagella_Pulp_4 <- RenameCells(Pagella_Pulp_4, new.names = paste0(Pagella_Pulp_4$Dataset, "_", Cells(Pagella_Pulp_4)))
Pagella_Pulp_5 <- RenameCells(Pagella_Pulp_5, new.names = paste0(Pagella_Pulp_5$Dataset, "_", Cells(Pagella_Pulp_5)))

Pagella_Perio_1 <- RenameCells(Pagella_Perio_1, new.names = paste0(Pagella_Perio_1$Dataset, "_", Cells(Pagella_Perio_1)))
Pagella_Perio_2 <- RenameCells(Pagella_Perio_2, new.names = paste0(Pagella_Perio_2$Dataset, "_", Cells(Pagella_Perio_2)))
Pagella_Perio_3 <- RenameCells(Pagella_Perio_3, new.names = paste0(Pagella_Perio_3$Dataset, "_", Cells(Pagella_Perio_3)))
Pagella_Perio_4 <- RenameCells(Pagella_Perio_4, new.names = paste0(Pagella_Perio_4$Dataset, "_", Cells(Pagella_Perio_4)))
Pagella_Perio_5 <- RenameCells(Pagella_Perio_5, new.names = paste0(Pagella_Perio_5$Dataset, "_", Cells(Pagella_Perio_5)))

#Krivanek and Song cell names were already adjusted.

xxx Hemeryck/Yin/Opasawatchai
```

```
# Group datasets together for QC and integration
Pagella_Pulp_1@meta.data$Dataset_merged <- "Pagella_Pulp"
Pagella_Pulp_2@meta.data$Dataset_merged <- "Pagella_Pulp"
Pagella_Pulp_3@meta.data$Dataset_merged <- "Pagella_Pulp"
Pagella_Pulp_4@meta.data$Dataset_merged <- "Pagella_Pulp"
Pagella_Pulp_5@meta.data$Dataset_merged <- "Pagella_Pulp"

Pagella_Perio_1@meta.data$Dataset_merged <- "Pagella_Periodontal"
Pagella_Perio_2@meta.data$Dataset_merged <- "Pagella_Periodontal"
Pagella_Perio_3@meta.data$Dataset_merged <- "Pagella_Periodontal"
Pagella_Perio_4@meta.data$Dataset_merged <- "Pagella_Periodontal"
Pagella_Perio_5@meta.data$Dataset_merged <- "Pagella_Periodontal"

Krivanek_Apical_Papilla_1@meta.data$Dataset_merged <- "Krivanek_Apical_Papilla"
Krivanek_Apical_Papilla_2@meta.data$Dataset_merged <- "Krivanek_Apical_Papilla"
Krivanek_Dental_Pulp@meta.data$Dataset_merged <- "Krivanek_Pulp"
Krivanek_Whole_Molar_1@meta.data$Dataset_merged <- "Krivanek_Whole_Molar"
Krivanek_Whole_Molar_2@meta.data$Dataset_merged <- "Krivanek_Whole_Molar"
Krivanek_Whole_Molar_3@meta.data$Dataset_merged <- "Krivanek_Whole_Molar"

Song_Tooth_Germ@meta.data$Dataset_merged <- "Song_Tooth_Germ"

Hemeryck_follicle_1@meta.data$Dataset_merged <- "Hemeryck_Dental_Follicle"
Hemeryck_follicle_2@meta.data$Dataset_merged <- "Hemeryck_Dental_Follicle"

Yin_Healthy@meta.data$Dataset_merged <- "Yin_Healthy"

Opasawatchai_Healthy@meta.data$Dataset_merged <- "Opasawatchai_Healthy"
```

## Healthy HTA quality control

```
merged <- merge(x = Pagella_Pulp_1, y = c(Pagella_Pulp_2, Pagella_Pulp_3, Pagella_Pulp_4,Pagella_Pulp_5,
                                          Pagella_Perio_1, Pagella_Perio_2, Pagella_Perio_3, Pagella_Perio_4, Pagella_Perio_5,
                                          Krivanek_Apical_Papilla_1, Krivanek_Apical_Papilla_2, Krivanek_Dental_Pulp, Krivanek_Whole_Molar_1, Krivanek_Whole_Molar_2, Krivanek_Whole_Molar_3, Song_Tooth_Germ, Hemeryck_follicle_1, Hemeryck_follicle_2, Opasawatchai_Healthy, Yin_Healthy))
                                          
tooth.list <- SplitObject(merged, split.by = "Dataset_merged")    

Pagella_Pulp <- tooth.list$Pagella_Pulp
Pagella_Periodontal <- tooth.list$Pagella_Periodontal
Krivanek_Apical_Papilla <- tooth.list$Krivanek_Apical_Papilla
Krivanek_Pulp <- tooth.list$Krivanek_Pulp
Krivanek_Whole_Molar <- tooth.list$Krivanek_Whole_Molar
Song_Tooth_Germ <- tooth.list$Song_Tooth_Germ
Hemeryck_Dental_Follicle <- tooth.list$Hemeryck_Dental_Follicle
Opasawatchai_Healthy <- tooth.list$Opasawatchai_Healthy
Yin_Healthy <- tooth.list$Yin_Healthy
```

```
#Perform QC individually
Pagella_Pulp <- subset(Pagella_Pulp, subset = nFeature_RNA > 800 & nFeature_RNA < 3000 & percent.mito < 15 & nCount_RNA < 40000)
Pagella_Periodontal <- subset(Pagella_Periodontal, subset = nFeature_RNA > 750 & nFeature_RNA < 4000 & percent.mito < 10 & nCount_RNA < 40000)
Krivanek_Apical_Papilla <- subset(Krivanek_Apical_Papilla, subset = nFeature_RNA > 500 & nFeature_RNA < 1500 & percent.mito < 6 & nCount_RNA < 40000)
Krivanek_Pulp <- subset(Krivanek_Pulp, subset = nFeature_RNA > 400 & nFeature_RNA < 1500 & percent.mito < 5 & nCount_RNA < 40000)
Krivanek_Whole_Molar <- subset(Krivanek_Whole_Molar, subset = nFeature_RNA > 1000 & nFeature_RNA < 3000 & percent.mito < 10 & nCount_RNA < 40000)
Song_Tooth_Germ <- subset(Song_Tooth_Germ, subset = nFeature_RNA > 750 & nFeature_RNA < 2500 & percent.mito < 20 & nCount_RNA < 40000)
Hemeryck_Dental_Follicle <- subset(Hemeryck_Dental_Follicle, subset = nFeature_RNA > 450 & nFeature_RNA < 5500 & percent.mito < 10 & nCount_RNA < 40000)
Opasawatchai_Healthy <- subset(Opasawatchai_Healthy, subset = nFeature_RNA > 900 & nFeature_RNA < 3000 & percent.mito < 3 & nCount_RNA < 40000)
Yin_Healthy <- subset(Yin_Healthy, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mito < 15 & nCount_RNA < 40000)
```

```
#Merge datasets
healthy_merged_final_qc <- merge(x = Pagella_Pulp, y = c(Pagella_Periodontal, Krivanek_Apical_Papilla, Krivanek_Pulp, Krivanek_Whole_Molar,
                                                         Song_Tooth_Germ, Hemeryck_Dental_Follicle, Opasawatchai_Healthy, Yin_Healthy))
```

## Healthy HTA initial integration and cell cycle regression

```
# Cell cycle gene sets rom the Seurat package:
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
```

```
# After acquiring the data, we first perform standard normalization and variable feature selection on the list of objects

tooth.list <- SplitObject(healthy_merged_final_qc, split.by = "Dataset_merged")
tooth.list <- lapply(X = tooth.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
})
```

```
# select features that are repeatedly variable across datasets for integration run PCA on each dataset using these features

features <- SelectIntegrationFeatures(object.list = tooth.list)
tooth.list <- lapply(X = tooth.list, FUN = function(x) {
    x <- ScaleData(x, features = features, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
```

```
# Perform integration
tooth.anchors <- FindIntegrationAnchors(object.list = tooth.list, dims = 1:30, anchor.features = features, reduction = "rpca")

healthy_integrated <- IntegrateData(anchorset = tooth.anchors, dims = 1:30)
```

```
# switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(healthy_integrated) <- "integrated"
```

```
# Run the standard workflow for visualization and clustering
healthy_integrated <- ScaleData(healthy_integrated, verbose = FALSE)
healthy_integrated <- RunPCA(healthy_integrated, verbose = FALSE, npcs = 100)
```

```
# Perform UMAP dimensional reduction using umap-learn
integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:40)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_40dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP40_", assay = DefaultAssay(integrated))

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:50)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_50dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP50_", assay = DefaultAssay(integrated))

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:60)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_60dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP60_", assay = DefaultAssay(integrated))

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:80)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_80dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP80_", assay = DefaultAssay(integrated))

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:100)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_100dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP100_", assay = DefaultAssay(integrated))
```

```
options(repr.plot.width=20, repr.plot.height=7)

p1 <- DimPlot(integrated, reduction = "umap_40dims", group.by = "Dataset_merged") +ggtitle("UMAP 40dims")
p2 <- DimPlot(integrated, reduction = "umap_40dims", group.by = "Tissue", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_50dims", group.by = "Dataset_merged") +ggtitle("UMAP 50dims")
p2 <- DimPlot(integrated, reduction = "umap_50dims", group.by = "Tissue", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

options(repr.plot.width=20, repr.plot.height=7)
p1 <- DimPlot(integrated, reduction = "umap_60dims", group.by = "Dataset_merged") +ggtitle("UMAP 60dims")
p2 <- DimPlot(integrated, reduction = "umap_60dims", group.by = "Tissue", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_80dims", group.by = "Dataset_merged") +ggtitle("UMAP 80dims")
p2 <- DimPlot(integrated, reduction = "umap_80dims", group.by = "Tissue", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_100dims", group.by = "Dataset_merged") +ggtitle("UMAP 100dims")
p2 <- DimPlot(integrated, reduction = "umap_100dims", group.by = "Tissue", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)
```

We continue with 80 dimensions

```
# Perform subclustering 
DefaultAssay(integrated) <- "integrated" #needs to be set to integrated because this is what we ran PCA on
integrated <- FindNeighbors(integrated, dims = 1:80)
integrated <- FindClusters(integrated, resolution = c(1.6, 1.8, 2, 5))
```

```
options(repr.plot.width=15, repr.plot.height=10)

DimPlot(integrated, reduction = "umap_80dims", group.by = "integrated_snn_res.1.6", label = TRUE) +
    labs(title = "80 dims UMAP with 1.6 res CCA")
# res 1.8
DimPlot(integrated, reduction = "umap_80dims", group.by = "integrated_snn_res.1.8", label = TRUE) +
    labs(title = "80 dims UMAP with 1.8 res CCA")
# res 2
DimPlot(integrated, reduction = "umap_80dims", group.by = "integrated_snn_res.2", label = TRUE) +
    labs(title = "80 dims UMAP with 2 res CCA")
# res 5
DimPlot(integrated, reduction = "umap_80dims", group.by = "integrated_snn_res.5", label = TRUE) +
    labs(title = "80 dims UMAP with 5 res CCA")
```

We continue with resolution = 5 

```
# Rough, initial annotation
Idents(integrated) <- integrated@meta.data$integrated_snn_res.5
integrated <- RenameIdents(integrated, `54` = "Naive_B", `58` = "B", `79` = "pDC",
                              `40` = "Endothelial", `55` = "Endothelial", `51` = "Endothelial",
                              `4` = "Endothelial", `12` = "Endothelial", `75` = "Endothelial",
                              `59` = "Endothelial", `68` = "Endothelial", `66` = "Endothelial",
                              `32` = "Endothelial", `17` = "Endothelial", `7` = "Endothelial",
                              `38` = "Endothelial", `52` = "Endothelial", `80` = "Schwann/Oligo",
                              `23` = "Schwann", `69` = "Schwann", `20` = "Glial",
                              `34` = "Glial", `83` = "Glia cells", `72` = "Glia cells",
                              `33` = "Glial", `41` = "SMC", `21` = "SMC", `0` = "SMC",
                              `53` = "SMC", `61` = "SMC", `65` = "Perivascular", `49` = "Perivascular",
                              `5` = "Perivascular", `2` = "Perivascular", `24` = "PDL",
                              `81` = "Dental follicle", `22` = "Dental follicle", `71` = "Dental follicle",
                              `10` = "Dental follicle", `35` = "Dental follicle", `11` = "Distal pulp",
                              `14` = "Distal pulp", `70` = "Distal pulp", `3` = "Distal pulp",
                              `13` = "Distal pulp", `36` = "Distal pulp", `76` = "Distal pulp",
                              `15` = "Distal pulp", `43` = "Distal pulp", `44` = "Distal pulp",
                              `30` = "Distal pulp", `6` = "Distal pulp", `48` = "Distal pulp",
                              `46` = "Odontoblasts", `18` = "Apical papilla", `19` = "Apical papilla",
                              `45` = "Apical papilla", `27` = "Apical papilla", `37` = "Apical papilla",
                              `50` = "Apical pulp", `82` = "Apical pulp", `16` = "Apical pulp",
                              `62` = "NK cells", `42` = "NK cells", `73` = "NK cells", `1` = "T cells",
                              `56` = "T cells", `57` = "T cells", `9` = "T cells", `26` = "T cells",
                              `63` = "T cells", `78` = "Unknown 1", `64` = "Unknown 2", `77` = "Naive T cells",
                              `74` = "Cycling cells", `47` = "ERM",  `84` = "Mast cells", `28` = "Macrophages",
                              `31` = "Macrophages", `8` = "Macrophages", `60` = "Macrophages",
                              `67` = "Macrophages", `25` = "Macrophages", `29` = "Macrophages", `39` = "Macrophages")
integrated$Annotation <- Idents(integrated)
```

## Healthy HTA removal of background or ambient RNA using SoupX

```
# Adjust for SoupX
integrated$Annotation <- gsub(' ', '_', integrated$Annotation)
integrated$Annotation <- gsub('/', '_', integrated$Annotation)
```

```
#Change from integrated to RNA, otherwise next command will fail
DefaultAssay(integrated) <- "RNA"
```

```
all.raw.data_integrated <- as.matrix(GetAssayData(integrated, slot = "counts"))

list_cluster_integrated <- integrated@meta.data[[sprintf("Annotation")]]
names(list_cluster_integrated) <- integrated@assays[["RNA"]]@data@Dimnames[[2]]
list_cluster_integrated <- as.matrix(list_cluster_integrated)

umap_integrated <- (Embeddings(integrated, reduction = "umap_80dims"))
```

```
# Create a SoupChannel object; because we are starting from our Seurat object we do not have a file with empty droplets. We will import the raw counts matrix from Seurat both as the table of counts (toc) and table of droplets (tod)

toc <- all.raw.data_integrated
tod = all.raw.data_integrated
sc = SoupChannel(tod, toc, calcSoupProfile = FALSE)
```

```
# Calculate the Soup profile (of ambient RNA)
toc = sc$toc
scNoDrops = SoupChannel(toc, toc, calcSoupProfile = FALSE)
soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), 
    counts = rowSums(toc))
scNoDrops = setSoupProfile(scNoDrops, soupProf)
```

```
# Add UMAP coordinates and cluster information from Seurat analysis to the SoupChannel object
scNoDrops = setClusters(scNoDrops, setNames(list_cluster_integrated[,1], rownames(list_cluster_integrated)))
scNoDrops = setDR(scNoDrops, umap_integrated[,c('UMAP80_1','UMAP80_2')])
```

```
# Estimate the contamination fraction
options(repr.plot.width=7, repr.plot.height=7)
sc = autoEstCont(scNoDrops, verbose = TRUE)=
```

```
# Remove the calculated contamination fraction from the original counts matrix, and add back to the original Seurat object. 
out = adjustCounts(sc)
integrated[["SoupX"]] <- CreateAssayObject(counts = out)
```

```
# Plot the change in expression due to the correction, e.g. for Ambn and Amelx
options(repr.plot.width=15, repr.plot.height=9)
x1 <- plotChangeMap(sc, out, "CD4")
x2 <- plotChangeMap(sc, out, "CD40")
x3 <- plotChangeMap(sc, out, "VIM")
x4 <- plotChangeMap(sc, out, "HBB")
(x1+x2)/(x3+x4)
```

```
# Normalize corrected counts, find variable features and scale data
DefaultAssay(integrated) <- "SoupX"
integrated <- NormalizeData(integrated)
integrated <- FindVariableFeatures(integrated, selection.method = "vst", nfeatures = 2000)
integrated <- ScaleData(integrated, verbose = FALSE)
```

## Healthy HTA rPCA integration

```
# Cell cycle gene sets rom the Seurat package:
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
```

```
# After acquiring the data, we first perform standard normalization and variable feature selection on the list of objects

DefaultAssay(integrated) <- "SoupX"
tooth.list <- SplitObject(integrated, split.by = "Dataset_merged")
tooth.list <- lapply(X = tooth.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
})
```

```
# select features that are repeatedly variable across datasets for integration run PCA on each dataset using these features

features <- SelectIntegrationFeatures(object.list = tooth.list)
tooth.list <- lapply(X = tooth.list, FUN = function(x) {
    x <- ScaleData(x, features = features, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
```

```
# Perform integration
tooth.anchors <- FindIntegrationAnchors(object.list = tooth.list, dims = 1:30, anchor.features = features, reduction = "rpca")

healthy_integrated <- IntegrateData(anchorset = tooth.anchors, dims = 1:30)
```

```
# switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(healthy_integrated) <- "integrated"
```

```
# Run the standard workflow for visualization and clustering
healthy_integrated <- ScaleData(healthy_integrated, verbose = FALSE)
healthy_integrated <- RunPCA(healthy_integrated, verbose = FALSE, npcs = 100)
```

```
# Perform UMAP dimensional reduction using umap-learn
integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:40)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_40dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP40_", assay = DefaultAssay(integrated))

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:50)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_50dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP50_", assay = DefaultAssay(integrated))

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:60)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_60dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP60_", assay = DefaultAssay(integrated))

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:80)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_80dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP80_", assay = DefaultAssay(integrated))

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:100)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_100dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP100_", assay = DefaultAssay(integrated))
```

```
options(repr.plot.width=20, repr.plot.height=7)

p1 <- DimPlot(integrated, reduction = "umap_40dims", group.by = "Dataset_merged") +ggtitle("UMAP 40dims")
p2 <- DimPlot(integrated, reduction = "umap_40dims", group.by = "Tissue", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_50dims", group.by = "Dataset_merged") +ggtitle("UMAP 50dims")
p2 <- DimPlot(integrated, reduction = "umap_50dims", group.by = "Tissue", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

options(repr.plot.width=20, repr.plot.height=7)
p1 <- DimPlot(integrated, reduction = "umap_60dims", group.by = "Dataset_merged") +ggtitle("UMAP 60dims")
p2 <- DimPlot(integrated, reduction = "umap_60dims", group.by = "Tissue", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_80dims", group.by = "Dataset_merged") +ggtitle("UMAP 80dims")
p2 <- DimPlot(integrated, reduction = "umap_80dims", group.by = "Tissue", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_100dims", group.by = "Dataset_merged") +ggtitle("UMAP 100dims")
p2 <- DimPlot(integrated, reduction = "umap_100dims", group.by = "Tissue", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)
```

We continue with 80 dimensions

```
# Perform subclustering 
DefaultAssay(integrated) <- "integrated" #needs to be set to integrated because this is what we ran PCA on
integrated <- FindNeighbors(integrated, dims = 1:80)
integrated <- FindClusters(integrated, resolution = c(1.6, 1.8, 2, 5))
```

```
options(repr.plot.width=15, repr.plot.height=10)

DimPlot(integrated, reduction = "umap_80dims", group.by = "integrated_snn_res.1.6", label = TRUE) +
    labs(title = "80 dims UMAP with 1.6 res CCA")
# res 1.8
DimPlot(integrated, reduction = "umap_80dims", group.by = "integrated_snn_res.1.8", label = TRUE) +
    labs(title = "80 dims UMAP with 1.8 res CCA")
# res 2
DimPlot(integrated, reduction = "umap_80dims", group.by = "integrated_snn_res.2", label = TRUE) +
    labs(title = "80 dims UMAP with 2 res CCA")
# res 5
DimPlot(integrated, reduction = "umap_80dims", group.by = "integrated_snn_res.5", label = TRUE) +
    labs(title = "80 dims UMAP with 5 res CCA")
```

We continues with res = 5 

```
Idents(integrated) <- integrated@meta.data$integrated_snn_res.5


integrated <- RenameIdents(integrated, `57` = "Naive B cells", `58` = "B cells", `78` = "PDC",
                              `15` = "Endothelial", `55` = "Endothelial", `61` = "Endothelial",
                              `47` = "Endothelial", `75` = "Endothelial", `1` = "Endothelial",
                              `39` = "Endothelial", `67` = "Endothelial", `52` = "Endothelial",
                              `11` = "Endothelial", `32` = "Endothelial", `12` = "Endothelial",
                              `56` = "Endothelial", `68` = "Endothelial", `69` = "Schwann/Oligo",
                              `22` = "Schwann/Oligo", `36` = "Glia cells",
                              `19` = "Glia cells", `79` = "Glia cells", `30` = "Glia cells",
                              `72` = "Glia cells", `40` = "SMC", `23` = "SMC", `16` = "SMC",
                              `50` = "SMC", `45` = "SMC", `60` = "SMC", `49` = "Perivascular",
                              `6` = "Perivascular", `3` = "Perivascular", `59` = "Perivascular", `34` = "PDL",
                              `71` = "PDL",`81` = "Dental follicle", `25` = "Dental follicle", `53` = "Dental follicle",
                              `13` = "Dental follicle", `44` = "Dental follicle", `10` = "Distal pulp",
                              `4` = "Distal pulp", `70` = "Distal pulp", `14` = "Distal pulp",
                              `8` = "Distal pulp", `27` = "Distal pulp", `33` = "Distal pulp",
                              `18` = "Distal pulp", `17` = "Distal pulp", `9` = "Distal pulp",
                              `54` = "Distal pulp",
                              `46` = "Odontoblasts", `37` = "Apical papilla", `5` = "Apical papilla",
                              `43` = "Apical papilla", `48` = "Apical papilla", `21` = "Apical papilla",
                              `80` = "Apical pulp", `31` = "Apical pulp", `35` = "Apical pulp",
                              `62` = "NK cells", `42` = "NK cells", `73` = "NK cells", `64` = "NK cells", `0` = "T cells",
                              `51` = "T cells", `7` = "T cells", `66` = "T cells", `29` = "T cells",
                              `78` = "PDC", `76` = "Unknown 1", `77` = "Unknown 2", `65` = "Unknown 3",
                              `74` = "Cycling cells", `41` = "ERM",  `82` = "Mast cells", `24` = "Macrophages",
                              `26` = "Macrophages", `20` = "Macrophages", `2` = "Macrophages",
                              `63` = "Macrophages", `38` = "Macrophages", `28` = "Macrophages")
seu$Annotation <- Idents(seu)
```

## Healthy HTA subclustering of dental epithelium

```
epithelium <- subset(x = seu, idents = c("Epithelial", "Cycling cells"))
DefaultAssay(epithelium) <- "SoupX"
```

Integration doesn't work because there are too few cells per dataset. Instead we will merge several together. Combined as follows:
- Pagella_Periodontal
- All others

```
ents(epithelium) <- epithelium@meta.data$Dataset_merged

epithelium <- RenameIdents(epithelium, `Hemeryck_Dental_Follicle` = "Group1", `Krivanek_Apical_Papilla` = "Group1", 
                    `Krivanek_Whole_Molar` = "Group1", `Pagella_Periodontal` = "Group2",
                    `Pagella_Pulp` = "Group1", `Yin_Healthy` = "Group1", `Song_Tooth_Germ` = "Group1")

epithelium@meta.data$Integration_Groups <- Idents(epithelium)
```

```
tooth.list <- SplitObject(epithelium, split.by = "Integration_Groups")
tooth.list <- lapply(X = tooth.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})
```

```
features <- SelectIntegrationFeatures(object.list = tooth.list)
tooth.list <- lapply(X = tooth.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE, npcs = 30, approx = FALSE) # approx = false for small datasets; see: https://github.com/satijalab/seurat/issues/1963
})
```

```
# Perform integration
tooth.anchors <- FindIntegrationAnchors(object.list = tooth.list, dims = 1:30, anchor.features = features, reduction = "rpca")
integrated <- IntegrateData(anchorset = tooth.anchors, dims = 1:30)
```

```
# switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(integrated) <- "integrated"
```

```
# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, verbose = FALSE)
```

```
# Run UMAP for different dimensions, save to object separately to be able to go back

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:10)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_10dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP10_", assay = DefaultAssay(integrated))

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:20)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_20dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP20_", assay = DefaultAssay(integrated))

integrated <- RunUMAP(integrated, umap.method = "umap-learn", dims = 1:30)
umap_dims <- integrated@reductions$umap@cell.embeddings
integrated[["umap_30dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP30_", assay = DefaultAssay(integrated))
```

```
options(repr.plot.width=20, repr.plot.height=7)

p1 <- DimPlot(integrated, reduction = "umap_10dims", group.by = "Dataset_merged") +ggtitle("umap_10dims")
p2 <- DimPlot(integrated, reduction = "umap_10dims", group.by = "Tissue", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_20dims", group.by = "Dataset_merged") +ggtitle("umap_20dims")
p2 <- DimPlot(integrated, reduction = "umap_20dims", group.by = "Tissue", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

options(repr.plot.width=20, repr.plot.height=7)
p1 <- DimPlot(integrated, reduction = "umap_30dims", group.by = "Dataset_merged") +ggtitle("umap_30dims")
p2 <- DimPlot(integrated, reduction = "umap_30dims", group.by = "Tissue", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)
```

```
options(repr.plot.width=20, repr.plot.height=7)

p1 <- DimPlot(integrated, reduction = "umap_10dims", group.by = "Integration_Groups") +ggtitle("umap_10dims")
p2 <- DimPlot(integrated, reduction = "umap_10dims", group.by = "integrated_snn_res.5", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(integrated, reduction = "umap_20dims", group.by = "Integration_Groups") +ggtitle("umap_20dims")
p2 <- DimPlot(integrated, reduction = "umap_20dims", group.by = "integrated_snn_res.5", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)

options(repr.plot.width=20, repr.plot.height=7)
p1 <- DimPlot(integrated, reduction = "umap_30dims", group.by = "Integration_Groups") +ggtitle("umap_30dims")
p2 <- DimPlot(integrated, reduction = "umap_30dims", group.by = "integrated_snn_res.5", label = TRUE, 
    repel = TRUE)
plot_grid(p1, p2)
```

We continue with dims = 30 

```
# Perform subclustering 
DefaultAssay(integrated) <- "integrated" #needs to be set to integrated because this is what we ran PCA on
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = c(0.2, 0.3, 0.6, 0.9, 1.2))
```

```
resolution <- c('integrated_snn_res.0.2', "integrated_snn_res.0.3", "integrated_snn_res.0.6", "integrated_snn_res.0.9", "integrated_snn_res.1.2")
suppressMessages({
plot.list <- list()
for (i in 1:length(resolution)) {
    options(repr.plot.width=7, repr.plot.height=5)
    plot.list[[i]] <- DimPlot(integrated, reduction = "umap_30dims", group.by = resolution[i], label = TRUE) +
        ggtitle(paste(resolution[i], sep = ""))
    print(plot.list[[i]])
}
    })
```

We continue with res = 1.2

```
Idents(integrated) <- integrated@meta.data$integrated_snn_res.1.2
integrated <- RenameIdents(integrated, `0` = "JE", `1` = "Pulp-Epi", `2` = "Perio-ERM",
                              `3` = "JE", `4` = "DF-ERM", `5` = "Perio-ERM",
                              `6` = "AP-EMT", `7` = "AP-Epi")
integrated@meta.data$epithelial <- Idents(integrated)
```

### DEG analysis of dental epithelial subclusters

```
DefaultAssay(integrated) <- "SoupX"
all.markers <- FindAllMarkers(object = integrated, logfc.threshold = 0.25, min.pct = 0.25)
```

### Analysis of GRNs using pySCENIC

First we extract the required data from our integrated Seurat object in our R environment. Afterwards we switch to our Python Conda environment containing pySCENIC. After setup, analysis can be divided into several steps: 
  (1) inference of co-expression modules
  (2) generation of regulons
  (3) Calculate AUC (i.e. cellular regulon enrichment)
  (4) Create a loom file (to share data via SCope, and to calculate RSS scores)
  (5) Calculate RSS scores

Afterwards, we move back to our R environment for further analyses using the pySCENIC output:
  (6) Create a regulon-based Seurat object and perform regulon-based integration
  (7) Calculate MRA

#### Extract required data
```
normCounts <- integrated[["SoupX"]]@data
write.csv(normCounts, "./HTA_SCENIC/ERM/resources_folder/normCounts_4scenic.csv")

annotation <- integrated@meta.data
write.csv(annotation, "./HTA_SCENIC/ERM/resources_folder/annotation_4scenic.csv")

umap_dims <- integrated@reductions$umap_50dims@cell.embeddings
write.csv(umap_dims, "./HTA_SCENIC/ERM/resources_folder/dims_4scenic.csv")
```

#### Switch to the pySCENIC Conda environment and set up workspace
```
import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import seaborn as sns
```

```
os.chdir("/HTA_SCENIC/ERM")
```

```
DATA_FOLDER="/HTA_SCENIC/ERM/data_folder"
RESOURCES_FOLDER="/HTA_SCENIC/ERM/resources_folder"
DATABASE_FOLDER = "/HTA_SCENIC/ERM/database_folder"
SCHEDULER="123.122.8.24:8786"
DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "hg38*.mc9nr.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, TF_names_v_1.01.txt')
SC_EXP_FNAME = os.path.join(RESOURCES_FOLDER, "normCounts_4scenic.csv")
REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.p")
MOTIFS_FNAME = os.path.join(DATA_FOLDER, "motifs.csv")
```
NB: The motifs and feather files were obtained from https://resources.aertslab.org/cistarget/ 

NB2: The list of TFs can be found from our resources folder in the repository. xxx

#### Load in counts, TF and ranking databases. 
```
exp_mtx = pd.read_csv(os.path.join(RESOURCES_FOLDER, "normCounts_4scenic.csv"), index_col=0, sep=',').T
exp_mtx.shape #check file
exp_mtx.iloc[0:10,100:110] #check file
```
```
tf_names = load_tf_names(MM_TFS_FNAME)
```
```
db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
dbs
```

#### (1) inference of co-expression modules
```
from distributed import Client, LocalCluster
local_cluster = LocalCluster(n_workers=10, threads_per_worker=1)
client = Client(local_cluster)
```
```
adjacencies = grnboost2(
    expression_data=exp_mtx,
    tf_names=tf_names,
    verbose=True,
    client_or_address=client)
    
adjacencies.to_csv("adjacencies_norm.csv")
adjacencies.head()
```
```
adjacencies = pd.read_csv('adjacencies_norm.csv', index_col=0)
adjacencies.head()
```

#### (2) generation of regulons
```
# Generate modules based on calculated adjacencies
modules = list(modules_from_adjacencies(adjacencies, exp_mtx))
```
```
# Prune modules using cic-regulatory information (RcisTarget)
with ProgressBar():
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME, client_or_address=client)
    
regulons = df2regulons(df)
len(regulons)

df.to_csv(MOTIFS_FNAME)
```
```
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)
```

#### (3) Calculate AUC
```
regulons = df2regulons(df)
```
```
auc_mtx = aucell(exp_mtx, regulons, num_workers=8)
auc_mtx.to_csv("auc_mtx.csv")
```
```
auc_mtx = pd.read_csv('auc_mtx.csv', index_col=0)
auc_mtx.head()
```
```
auc_mtx.index = auc_mtx.index.str.replace('.', '-') #fix formatting of cell names
auc_mtx.index = auc_mtx.index.str.replace('X', '') #fix formatting of cell names, somehow and X was added to the cell names
auc_mtx.to_csv("auc_mtx_corrected_index.csv")
auc_mtx.head()
auc_mtx.shape
```

#### (4) Loom file generation
Set up environment
```
from pyscenic.export import export2loom
from pyscenic.utils import load_motifs, Sequence
from pyscenic.transform import df2regulons
```
```
# Import gene expression matrix, fix formatting of cell names.
exp_mtx = pd.read_csv(os.path.join(RESOURCES_FOLDER, "normCounts_4scenic.csv"), index_col=0, sep=',').T
exp_mtx.index = exp_mtx.index.str.replace('.', '-')
exp_mtx.head()
```
```
# Check if expression matrix has the correct format
def is_valid_exp_matrix(mtx):
    return (all(isinstance(idx, str) for idx in mtx.index) 
            and all(isinstance(idx, str) for idx in mtx.columns)
            and (mtx.index.nlevels == 1)
            and (mtx.columns.nlevels == 1))
is_valid_exp_matrix(exp_mtx)
```
Import all required data, and ensure correct formatting. 
```
# Import motifs
motifs = df
```
```
# Import metadata from Seurat
annotation_loom = pd.read_csv("HTA_SCENIC/ERM/resources_folder/annotation_4scenic.csv", index_col=0, sep=',')
annotation_loom = annotation_loom.loc[:,"epithelial"] # Extract cell annotation
annotation_loom.index = annotation_loom.index.str.replace('.', '-') # fix formatting of cell names

annotation_loom.to_csv("_annotation_loom.csv")
!mv /MTA/SCENIC/annotation_loom.csv /HTA_SCENIC/ERM/resources_folder/
```
```
ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "FH_atlas_annotation_loom.csv")
with open(ANNOTATIONS_FNAME, "rt") as f:
     annotations = dict(line.strip().replace("\"", "").split(",") for idx, line in enumerate(f) if idx > 0)
```
```
#Check if annotation has the correct format
def is_valid_annotation_mapping(m):
    return (all(isinstance(k, str) for k in m.keys()) 
            and all(isinstance(v, str) for v in m.values()))
is_valid_annotation_mapping(annotations)
dict(list(annotations.items())[:5])
```
```
# Import UMAP dimensional reduction coordinates
embeddings = pd.read_csv("/HTA_SCENIC/ERM/resources_folder/UMAP_dims_4scenic.csv", index_col=0, sep=',')
embeddings.index = embeddings.index.str.replace('.', '-') # fix formatting of cell names
```
```
embeddings = embeddings.iloc[:,[0,1]]
embeddings.columns=['_X', '_Y']
```
```
embeddings = { "UMAP (default)" : pd.DataFrame(data=embeddings,
                                      index=exp_mtx.index, columns=['_X', '_Y']) } # (n_cells, 2)
                                      
```
```
from collections import OrderedDict
id2name = OrderedDict()
embeddings_X = pd.DataFrame(index=exp_mtx.index)
embeddings_Y = pd.DataFrame(index=exp_mtx.index)
for idx, (name, df_embedding) in enumerate(embeddings.items()):
    if(len(df_embedding.columns)!=2):
        raise Exception('The embedding should have two columns.')
```
```
embedding_id = idx - 1 # Default embedding must have id == -1 for SCope.
id2name[embedding_id] = name
```
```
embedding = df_embedding.copy()
embedding.columns=['_X', '_Y']
embeddings_X = pd.merge(embeddings_X, embedding['_X'].to_frame().rename(columns={'_X': str(embedding_id)}), left_index=True, right_index=True)
embeddings_Y = pd.merge(embeddings_Y, embedding['_Y'].to_frame().rename(columns={'_Y': str(embedding_id)}), left_index=True, right_index=True)
```
```
# Calculate the number of genes per cell
binary_mtx = exp_mtx.copy()
binary_mtx[binary_mtx != 0] = 1.0
ngenes = binary_mtx.sum(axis=1).astype(int)
```
```
# Encode genes in regulons as a binary membership matrix
from operator import attrgetter
genes = np.array(exp_mtx.columns)
n_genes = len(genes)
n_regulons = len(regulons)
data = np.zeros(shape=(n_genes, n_regulons), dtype=int)
for idx, regulon in enumerate(regulons):
    data[:, idx] = np.isin(genes, regulon.genes).astype(int)
regulon_assignment = pd.DataFrame(data=data,
                                    index=exp_mtx.columns,
                                    columns=list(map(attrgetter('name'), regulons)))
```
```
# Encode cell annotations
name2idx = dict(map(reversed, enumerate(sorted(set(annotations.values())))))

clusterings = pd.DataFrame(data=exp_mtx.index.values,
                               index=exp_mtx.index,
                               columns=['0']).replace(annotations).replace(name2idx)

clusterings.head() # Data for first cell always gets mixed up; add correct cluster manually (below)
clusterings.iloc[[0],[0]] = 34
clusterings.head()
```
```
# Encode cluster_labels (for RSS score calculation)
cluster_labels = pd.DataFrame(data=exp_mtx.index.values,
                               index=exp_mtx.index,
                               columns=['0']).replace(annotations)

cluster_labels.head() # Data for first cell always gets mixed up; add correct cluster manually (below)
cluster_labels.iloc[[0],[0]] = "Pulp-Epi"
```
Create metadata structure of loom file
```
def create_structure_array(df):
    return np.array([tuple(row) for row in df.values],
                    dtype=np.dtype(list(zip(df.columns, df.dtypes))))
```
```
    default_embedding = next(iter(embeddings.values())).copy()
    default_embedding.columns=['_X', '_Y']
    column_attrs = {
        "CellID": exp_mtx.index.values.astype('str'),
        "nGene": ngenes.values,
        "Embedding": create_structure_array(default_embedding),
        "RegulonsAUC": create_structure_array(auc_mtx),
        "Clusterings": create_structure_array(clusterings),
        "ClusterID": clusterings.values,
        "Cluster_labels":cluster_labels.values,
        'Embeddings_X': create_structure_array(embeddings_X),
        'Embeddings_Y': create_structure_array(embeddings_Y),
        }
    row_attrs = {
        "Gene": exp_mtx.columns.values.astype('str'),
        "Regulons": create_structure_array(regulon_assignment),
        }
```
```
def fetch_logo(context):
    for elem in context:
        if elem.endswith('.png'):
            return elem
    return ""
```
```
def fetch_logo(context):
    for elem in context:
        if elem.endswith('.png'):
            return elem
    return ""
name2logo = {reg.name: fetch_logo(reg.context) for reg in regulons}
regulon_thresholds = [{"regulon": name,
                        "defaultThresholdValue":(threshold if isinstance(threshold, float) else threshold[0]),
                        "defaultThresholdName": "guassian_mixture_split",
                        "allThresholds": {"guassian_mixture_split": (threshold if isinstance(threshold, float) else threshold[0])},
                        "motifData": name2logo.get(name, "")} for name, threshold in auc_mtx.iteritems()] 
```
```
import json
general_attrs = {
    "title": "RegulonsEpiHTA",
    "MetaData": json.dumps({
        "embeddings": [{'id': identifier, 'name': name} for identifier, name in id2name.items()],
        "annotations": [{
            "name": "",
            "values": []
        }],
        "clusterings": [{
            "id": 0,
            "group": "celltype",
            "name": "Cell Type",
            "clusters": [{"id": idx, "description": name} for name, idx in name2idx.items()]
        }],
        "regulonThresholds": regulon_thresholds
    }),
    "Genome": "mm10"}
```
Add tree structure
```
tree_structure: Sequence[str] = ()
```
```
from itertools import islice
import itertools
assert len(tree_structure) <= 3, ""
general_attrs.update(("SCopeTreeL{}".format(idx+1), category)
                        for idx, category in enumerate(list(islice(itertools.chain(tree_structure, itertools.repeat("")), 3))))
```
```
def compress_encode(value):
    '''
    Compress using ZLIB algorithm and encode the given value in base64.
    From: https://github.com/aertslab/SCopeLoomPy/blob/5438da52c4bcf48f483a1cf378b1eaa788adefcb/src/scopeloompy/utils/__init__.py#L7
    '''
    return base64.b64encode(zlib.compress(value.encode('ascii'))).decode('ascii')
```
```
import base64
import zlib
general_attrs["MetaData"] = compress_encode(value=general_attrs["MetaData"])
```
Create loom file:
```
import loompy as lp
lp.create(filename="HTA_Epi_Regulons_loom.loom",
              layers=exp_mtx.T.values,
              row_attrs=row_attrs,
              col_attrs=column_attrs,
              file_attrs=general_attrs)
```
```
from pyscenic.genesig import Regulon
from typing import Union
```
```
def add_scenic_metadata(adata: 'sc.AnnData',
                        auc_mtx: pd.DataFrame,
                        regulons: Union[None, Sequence[Regulon]] = None,
                        bin_rep: bool = False,
                        copy: bool = False) -> 'sc.AnnData':
    """
    Add AUCell values and regulon metadata to AnnData object.
    :param adata: The AnnData object.
    :param auc_mtx: The dataframe containing the AUCell values (#observations x #regulons).
    :param bin_rep: Also add binarized version of AUCell values as separate representation. This representation
    is stored as `adata.obsm['X_aucell_bin']`.
    :param copy: Return a copy instead of writing to adata.
    :
    """
    # To avoid dependency with scanpy package the type hinting intentionally uses string literals.
    # In addition, the assert statement to assess runtime type is also commented out.
    #assert isinstance(adata, sc.AnnData)
    assert isinstance(auc_mtx, pd.DataFrame)
    assert len(auc_mtx) == adata.n_obs

    REGULON_SUFFIX_PATTERN = 'Regulon({})'

    result = adata.copy() if copy else adata

    # Add AUCell values as new representation (similar to a PCA). This facilitates the usage of
    # AUCell as initial dimensional reduction.
    result.obsm['X_aucell'] = auc_mtx.values.copy()
    if bin_rep:
        bin_mtx, _ = binarize(auc_mtx)
        result.obsm['X_aucell_bin'] = bin_mtx.values

    # Encode genes in regulons as "binary" membership matrix.
    if regulons is not None:
        genes = np.array(adata.var_names)
        data = np.zeros(shape=(adata.n_vars, len(regulons)), dtype=bool)
        for idx, regulon in enumerate(regulons):
            data[:, idx] = np.isin(genes, regulon.genes).astype(bool)
        regulon_assignment = pd.DataFrame(data=data, index=genes,
                                          columns=list(map(lambda r: REGULON_SUFFIX_PATTERN.format(r.name), regulons1)))
        result.var = pd.merge(result.var, regulon_assignment, left_index=True, right_index=True, how='left')

    # Add additional meta-data/information on the regulons.
    def fetch_logo(context):
        for elem in context:
            if elem.endswith('.png'):
                return elem
        return ""
    result.uns['aucell'] = {
        'regulon_names': auc_mtx.columns.map(lambda s: REGULON_SUFFIX_PATTERN.format(s)).values,
        'regulon_motifs': np.array([fetch_logo(reg.context) for reg in regulons] if regulons is not None else [])
    }

    # Add the AUCell values also as annotations of observations. This way regulon activity can be
    # depicted on cellular scatterplots.
    mtx = auc_mtx.copy()
    mtx.columns = result.uns['aucell']['regulon_names']
    result.obs = pd.merge(result.obs, mtx, left_index=True, right_index=True, how='left')

    return result
```
```
def export_regulons(regulons: Sequence[Regulon], fname: str) -> None:
    """
    Export regulons as GraphML.
    :param regulons: The sequence of regulons to export.
    :param fname: The name of the file to create.
    """
    graph = nx.DiGraph()
    for regulon in regulons:
        src_name = regulon.transcription_factor
        graph.add_node(src_name, group='transcription_factor')
        edge_type = 'activating' if 'activating' in regulon.context else 'inhibiting'
        node_type = 'activated_target' if 'activating' in regulon.context else 'inhibited_target'
        for dst_name, edge_strength in regulon.gene2weight.items():
            graph.add_node(dst_name, group=node_type, **regulon.context)
            graph.add_edge(src_name, dst_name, weight=edge_strength, interaction=edge_type, **regulon.context)
    nx.readwrite.write_graphml(graph, fname)
```


#### (5) Calculate RSS scores
```
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from pyscenic.binarization import binarize
```
```
f_final_loom = 'HTA_Epi_Regulons_loom.loom'
```
```
lf = lp.connect( f_final_loom, mode='r', validate=False )
auc_mtx_loom = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
cellAnnot = pd.concat([pd.DataFrame( lf.ca.Cluster_labels, index=lf.ca.CellID )], axis=1)
lf.close()
```
```
rss_cellType = regulon_specificity_scores(auc_mtx, cellAnnot[0])
rss_cellType
```
```
rss_cellType.to_csv("HTA_EPI_rss_clusters.csv")
```
```
#RSS panel plot with all cell types

cats = sorted(list(set(cellAnnot[0])))

fig = plt.figure(figsize=(40, 40))
for c,num in zip(cats, range(1,len(cats)+1)):
    x=rss_clusters.T[c]
    ax = fig.add_subplot(5,7,num)
    plot_rss(rss_clusters, c, top_n=10, max_n=None, ax=ax)
    ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
 
fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
        'figure.titlesize': 'large' ,
        'axes.labelsize': 'medium',
        'axes.titlesize':'large',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'
        })
plt.savefig("clusters-RSS-top5.pdf", dpi=600, bbox_inches = "tight")
plt.show()
```

#### (6) Calculate MRA (R)
Move back to our R environment, as done before. 
```
integrated
```
```
auc_mtx <- read.csv("/HTA_SCENIC/ERM/auc_mtx.csv")
auc_mtx$Cell <- gsub('\\.', '-', auc_mtx$Cell)
```
```
UMAP_dims <- integrated@reductions$umap_30dims@cell.embeddings
annotation <- integrated@meta.data
```
```
anno <- data.frame(UMAP_dims[,-1],
                   Cell = UMAP_dims$X,
                   clusters = annotation$epithelial,
                   dataset = annotation$Dataset_merged,
                   extraction = annotation$Extraction,
                   tissue = annotation$Tissue,
                  row.names = UMAP_dims$X)
```
```
# Join the data
regulon_anno <- inner_join(auc_mtx, anno, by = "Cell")
regulon_anno_long <- gather(regulon_anno, regulon, activity, -clusters, -dataset, -tissue, -extraction, -Cell, -UMAP30_1, -UMAP30_2)
```

```
# Remove 3 dots at the end of regulon's name
regulon_anno_long$regulon <- gsub('.{0,3}$', '', regulon_anno_long$regulon)
```
Calculate MRA:
```
meanRegPerCluster <- regulon_anno_long %>%
                        dplyr::select(-c(Cell, UMAP30_1, UMAP30_2)) %>%
                        group_by(clusters, regulon) %>%
                        summarize(mean_activity = mean(activity))
```

#### (7) Create a regulon-based Seurat object and perform regulon-based integration (R)
```
auc_mtx <- auc_mtx %>%
                dplyr::filter(Cell %in% regulon_anno$Cell)
```
```
# Remove 3 dots at the end of regulon's name
colnames(auc_mtx) <- gsub('\\...','', colnames(auc_mtx))
```
```
auc_mtx <- data.frame(auc_mtx[,-1], row.names = auc_mtx[,1])
```
```
# Transpose matrix
auc_mtx_T <- t(auc_mtx)
```
```
meta_data <- data.frame(annotation)
```
```
seurat <- CreateSeuratObject(counts = auc_mtx_T, meta.data = meta_data, min.cells = 0, min.features = 0, project = "AUC")
```
Perform standard integration
```
tooth.list <- SplitObject(seurat, split.by = "Dataset")
tooth.list <- lapply(X = tooth.list, FUN = function(x) {
    x <- FindVariableFeatures(x, verbose = FALSE)
})
```
```
tooth.anchors <- FindIntegrationAnchors(object.list = tooth.list, dims = 1:30)
integrated_regulons <- IntegrateData(anchorset = tooth.anchors, dims = 1:30)
```
```
DefaultAssay(integrated_regulons) <- "integrated_regulons"
integrated_regulons <- ScaleData(integrated_regulons, verbose = FALSE, vars.to.regress = c("S.Score", "G2M.Score"))
integrated_regulons <- RunPCA(integrated_regulons, npcs = 100, verbose = FALSE)
```
```
# Run UMAP for different dimensions, save to object separately to be able to go back
integrated_regulons <- RunUMAP(integrated_regulons, umap.method = "umap-learn", dims = 1:10, min.dist = 0.5)
umap_dims <- integrated_regulons@reductions$umap@cell.embeddings
integrated_regulons[["umap_10dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP10_", assay = DefaultAssay(integrated_regulons))

integrated_regulons <- RunUMAP(integrated_regulons, umap.method = "umap-learn", dims = 1:20, min.dist = 0.5)
umap_dims <- integrated_regulons@reductions$umap@cell.embeddings
integrated_regulons[["umap_20dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP20_", assay = DefaultAssay(integrated_regulons))

integrated_regulons <- RunUMAP(integrated_regulons, umap.method = "umap-learn", dims = 1:30, min.dist = 0.5)
umap_dims <- integrated_regulons@reductions$umap@cell.embeddings
integrated_regulons[["umap_30dims"]] <- CreateDimReducObject(embeddings = umap_dims, key = "UMAP30_", assay = DefaultAssay(integrated_regulons))
```
```
options(repr.plot.width=7, repr.plot.height=5)
DimPlot(integrated_regulons, reduction = "umap_10dims", group.by = "Annotation", pt.size = 2, label = TRUE, repel = TRUE)
DimPlot(integrated_regulons, reduction = "umap_20dims", group.by = "Annotation", pt.size = 2, label = TRUE, repel = TRUE)
DimPlot(integrated_regulons, reduction = "umap_30dims", group.by = "Annotation", pt.size = 2, label = TRUE, repel = TRUE)
```
```
options(repr.plot.width=7, repr.plot.height=5)
DimPlot(integrated_regulons, reduction = "umap_10dims", group.by = "epithelial", pt.size = 2, label = TRUE, repel = TRUE)
DimPlot(integrated_regulons, reduction = "umap_20dims", group.by = "epithelial", pt.size = 2, label = TRUE, repel = TRUE)
DimPlot(integrated_regulons, reduction = "umap_30dims", group.by = "epithelial", pt.size = 2, label = TRUE, repel = TRUE)
```
Import pySCENIC results to original Seurat object of subclustered dental epithelium:
```
integrated[["Regulons"]] <- CreateAssayObject(counts = auc_mtx_T)
DefaultAssay(integrated) <- "Regulons"
```

## Diseased HTA setup 

```
GSM5509264 <- Read10X(data.dir = "./Tong_GEO_2021/SampleA") #Lin dataset
Tong_CAP_1 <- CreateSeuratObject(counts = GSM5509264, project = "GSM5509264", min.cells = 3, min.features = 200)

GSM5509265 <- Read10X(data.dir = "./Tong_GEO_2021/SampleB") #Lin dataset
Tong_CAP_2 <- CreateSeuratObject(counts = GSM5509265, project = "GSM5509265", min.cells = 3, min.features = 200)

GSM5509266 <- Read10X(data.dir = "./Tong_GEO_2021/SampleC") #Lin dataset
Tong_Periapical_Granuloma <- CreateSeuratObject(counts = GSM5509266, project = "GSM5509266", min.cells = 3, min.features = 200)
```

```
Tong_CAP_1[["percent.mito"]] <- PercentageFeatureSet(Tong_CAP_1, pattern = "^MT-")
Tong_CAP_2[["percent.mito"]] <- PercentageFeatureSet(Tong_CAP_2, pattern = "^MT-")
Tong_Periapical_Granuloma[["percent.mito"]] <- PercentageFeatureSet(Tong_Periapical_Granuloma, pattern = "^MT-")
```

```
# Add extra metadata to Seurat object
Tong_CAP_1@meta.data$Dataset <- "Tong_CAP_1"
Tong_CAP_1@meta.data$Tissue <- "Periapical Tissue"
Tong_CAP_1@meta.data$Extraction <- "Maxillary Posterior Teeth Right"
Tong_CAP_1@meta.data$Species <- "Human"
Tong_CAP_1@meta.data$Age <- "26"
Tong_CAP_1@meta.data$Conditon <- "Chronical Apical Peridontitis"
Tong_CAP_1@meta.data$Technology <- "10X"

Tong_CAP_2@meta.data$Dataset <- "Tong_CAP_2"
Tong_CAP_2@meta.data$Tissue <- "Periapical Tissue"
Tong_CAP_2@meta.data$Extraction <- "Maxillary Posterior Teeth Right"
Tong_CAP_2@meta.data$Species <- "Human"
Tong_CAP_2@meta.data$Age <- "44"
Tong_CAP_2@meta.data$Conditon <- "Chronical Apical Peridontitis"
Tong_CAP_2@meta.data$Technology <- "10X"

Tong_Periapical_Granuloma@meta.data$Dataset <- "Tong_Periapical_Granuloma"
Tong_Periapical_Granuloma@meta.data$Tissue <- "Periapical Tissue"
Tong_Periapical_Granuloma@meta.data$Extraction <- "Mandible Posterior Teeth Left"
Tong_Periapical_Granuloma@meta.data$Species <- "Human"
Tong_Periapical_Granuloma@meta.data$Age <- "27"
Tong_Periapical_Granuloma@meta.data$Conditon <- "Periapical Granuloma"
Tong_Periapical_Granuloma@meta.data$Technology <- "10X"
```

```
# Add Dataset of origin to cell name to avoid possible identical cell names between datasets
Tong_CAP_1 <- RenameCells(Tong_CAP_1, new.names = paste0(Tong_CAP_1$Dataset, "_", Cells(Tong_CAP_1)))
Tong_CAP_2 <- RenameCells(Tong_CAP_2, new.names = paste0(Tong_CAP_2$Dataset, "_", Cells(Tong_CAP_2)))
Tong_Periapical_Granuloma <- RenameCells(Tong_Periapical_Granuloma, new.names = paste0(Tong_Periapical_Granuloma$Dataset, "_", Cells(Tong_Periapical_Granuloma)))

```

```
Tong_CAP_1@meta.data$Dataset_merged <- "Tong_CAP"
Tong_CAP_2@meta.data$Dataset_merged <- "Tong_CAP"
Tong_Periapical_Granuloma@meta.data$Dataset_merged <- "Tong_Periapical_Granuloma"

Opasawatchai_Deep_Carries_1@meta.data$Dataset_merged <- "Opasawatchai_Deep_Carries"
Opasawatchai_Deep_Carries_2@meta.data$Dataset_merged <- "Opasawatchai_Deep_Carries"
Opasawatchai_Enamel_Carries@meta.data$Dataset_merged <- "Opasawatchai_Enamel_Carries"
```

## Diseased HTA quality control

```
merged <- merge(x = Tong_CAP_1, y = c(Tong_CAP_2, Tong_Periapical_Granuloma, Opasawatchai_Deep_Carries_1, Opasawatchai_Deep_Carries_2, Opasawatchai_Enamel_Carries))

tooth.list <- SplitObject(merged, split.by = "Dataset_merged")

Tong_CAP <- tooth.list$Tong_CAP
Tong_Periapical_Granuloma <- tooth.list$Tong_Periapical_Granuloma
Opasawatchai_Deep_Carries <- tooth.list$Opasawatchai_Deep_Carries
Opasawatchai_Enamel_Carries <- tooth.list$Opasawatchai_Enamel_Carries
```

```
Tong_CAP <- subset(Tong_CAP, subset = nFeature_RNA > 1000 & nFeature_RNA < 3500 & percent.mito < 6 & nCount_RNA < 40000)
Tong_Periapical_Granuloma <- subset(Tong_Periapical_Granuloma, subset = nFeature_RNA > 1100 & nFeature_RNA < 3500 & percent.mito < 5 & nCount_RNA < 40000)
Opasawatchai_Deep_Carries <- subset(Opasawatchai_Deep_Carries, subset = nFeature_RNA > 900 & nFeature_RNA < 3500 & percent.mito < 3 & nCount_RNA < 40000)
Opasawatchai_Enamel_Carries <- subset(Opasawatchai_Enamel_Carries, subset = nFeature_RNA > 800 & nFeature_RNA < 3000 & percent.mito < 5 & nCount_RNA < 40000)

```

## Diseased HTA initial integration

```
#Merge datasets
diseased_merged_final_qc <- merge(x = Tong_CAP, y = c(Tong_Periapical_Granuloma, Opasawatchai_Deep_Carries,
                                                      Opasawatchai_Enamel_Carries))
```
