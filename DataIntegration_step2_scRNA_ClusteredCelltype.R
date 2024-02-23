#' @param dataset Input Seurat object for preclustered mono-cell type dataset 
#' @param path Output path 

ClusteredCelltype <- function(dataset, path) {
    library(dplyr)
    library(Seurat)
    library(ggplot2)
    library(doMC)
    library(harmony)
    library(RColorBrewer)
    registerDoMC(50)
    
    dir.create(paste0(path, gsub(".rds", "", x)))
    a <- NormalizeData(dataset, normalization.method = "LogNormalize", scale.factor = 10000)
    a <- FindVariableFeatures(a, selection.method = "vst", nfeatures = 2000)
    batch <- c("percent.mt", "percent.rp", "S.Score", "G2M.Score")
    if (length(unique(a@meta.data$Tech)) > 1) {
        batch <- c(batch, "Tech")
    } 
    if (length(unique(a@meta.data$SampleID)) > 1) {
        batch <- c(batch, "SampleID")
    } 

    a <- SCTransform(a, vars.to.regress = batch, verbose = T)
    a <- RunPCA(a, features = VariableFeatures(object = a))
    a <- RunUMAP(a, dims = 1:30, label = T)

    source("scRNA_Harmony.R")
    group.by.vars = c()
    if (length(unique(a@meta.data$Tech)) > 1) {
        group.by.vars <- c(group.by.vars, "Tech")
    } 
    if (length(unique(a@meta.data$SampleID)) > 1) {
        group.by.vars <- c(group.by.vars, "SampleID")
    } 
    if(length(group.by.vars)>1){
        a <- BatchRemove_Harmony(
            seurat_object = a, cores = 20, group.by.vars = group.by.vars,
            colors = c(brewer.pal(8, "Set1"), colorRampPalette(brewer.pal(8, "Accent"))(length(unique(a@meta.data$SampleID)))), max.dim = 30,
            out_dir = paste0(path, gsub(".rds", "", x))
        )        
    }

    # Calculation of marker genes
    source("Seurat_marker.R")
    if (length(unique(a@meta.data$Tech)) == 1) {
        logfc.threshold <- ifelse(unique(a@meta.data$Tech) == "smart-seq2", 1, 0.25)
    } else {
        logfc.threshold <- 0.25
    }
    a <- Seurat_marker(seurat_object = a, FindNeighbor_reduction = "harmony", resolution = 2, cores = 10, top_n = 50, out_dir = paste0(path, gsub(".rds", "", x)))
    return(a)

})

