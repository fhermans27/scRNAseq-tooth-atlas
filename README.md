# scRNAseq-tooth-atlas

This repository contains the R and Python code used perform the single-cell RNA-sequencing (scRNAseq) analysis used in :

*xxx add publication*

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

- [Human Tooth Atlas](#human-tooth-atlas)
  - [Setup](#hta-setup)
  - [Quality Control](#hta-quality-control)

# Mouse Tooth Atlas
## MTA Setup

```
suppressMessages({
#set working directory
setwd("/staging/leuven/stg_00073/SC_Florian/")
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

## MTA Quality Control
```
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
```# Merge datasets
merged_final_qc <- merge(x = Chiba, y = c(Takahashi, Sharir, Krivanek_incisors, Krivanek_molars, Chen, Wen, Chiba_2_epi, Chiba_2_mes, Nagata, Zhao_1, Zhao_2))
```

## Initial Integration on three groups 
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
xxx add mouse genes conversion
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
    plot.list[[i]] <- FeaturePlot(seu, reduction = "umap_50dims", features = DESC_features[i], pt.size = 1, cols = (c("lightgrey", "steelblue", "blue4")))
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
meta_data <- cbind(rownames(incisors@meta.data), incisors@meta.data[,"Annotation", drop=F])   # cluster is the user’s specific cluster column
#write.table(meta_data, "./incisor_cpdb_meta.txt", sep="\t", quote=F, row.names=F)
```

For the next part of the CellPhoneDB analysis of Ligand-Receptor interactions in the MTA we use our Python Conda environment with the CellPhoneDB tool installed, and perform the analysis using the command line. 

``` 
conda activate py_cpdb
```

```
# Example for incisor group, same done for molar/periodontal groups
cellphonedb method statistical_analysis --counts-data=gene_name incisor_cpdb_meta.txt incisor_cpdb_count.txt --iterations=1000 --threshold=0.2 --threads=8 --output-path=CellPhoneDB --project-name=MTA_incisor --quiet
```

Next, make a rows.txt and columns.txt file with the L/R pairs and cluster-cluster pairs for further analysis. Run the CellPhoneDB plotting function for each group.

```
# Example for incisor group, same done for molar/periodontal groups
cellphonedb plot dot_plot --means-path=./CellPhoneDB/MTA_incisor/means.txt --pvalues-path=./CellPhoneDB/MTA_incisor/pvalues.txt --rows ./CellPhoneDB/MTA_incisor/rows.txt --columns ./CellPhoneDB/MTA_incisor/columns.txt --output-path=./CellPhoneDB/MTA_incisor/ --output-name=cpdb_MTA_incisor_dotplot.pdf
```

# Human Tooth Atlas

## HTA Setup

## HTA Quality Control


