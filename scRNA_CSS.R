#' @param data Input Seurat object

CSS_scaling <- function(data){
    library(Seurat)
    library(dplyr)
    library(simspec) 
    a <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
    a <- FindVariableFeatures(a, selection.method = "vst", nfeatures = 2000)

    a <- ScaleData(a, features = VariableFeatures(a), vars.to.regress = c("percent.mt", "percent.rp", "S.Score", "G2M.Score","dig.score", "DonorID", "Tech"))
    a <- RunPCA(a, features = VariableFeatures(object = a))

    a <- cluster_sim_spectrum(object = a, label_tag = "DonorID", spectrum_type = "corr_ztransform",corr_method = "spearman")
    a <- RunUMAP(a, reduction = "css", dims = 1:ncol(Embeddings(a, "css")))
    return(a)
}
