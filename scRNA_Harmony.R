#' @TODO Perform batch correction using harmony
#' @param seurat_object Input Seurat object
#' @param cores Number of threads for parallel computation
#' @param group.by.vars Which variables need to be removed (one or more, as a vector)
#' @param colors Colors used in Dimplot
#' @param assay.use used assay, default:RNA
#' @param max.dim The maximum dimension for dimensionality reduction; the larger the dimension, the finer the subsequent cell clustering. Use the inflection point in the Elbow plot
#' @param out_dir Path for the output image
#' @returnType Seurat object
#' @return 
#' 
#' author LX
#' 
#' 

BatchRemove_Harmony <- function(seurat_object = NULL, cores = NULL, group.by.vars = NULL, colors = NULL, assay.use = "RNA", max.dim = NULL, out_dir = NULL){
    library(dplyr)
    library(Seurat)
    library(harmony)
    library(ggplot2)
    library(ggthemes)
    library(doMC)
    registerDoMC(cores)

    # Dimplot before batch correction
    for(i in 1:length(group.by.vars)){
        p1 <- DimPlot(object = seurat_object, reduction = "umap", pt.size = .2, group.by = group.by.vars[i], raster = FALSE, cols = colors) + coord_fixed()
        ggsave(filename=paste0(out_dir,"/umap_before_Harmony_",group.by.vars[i],".pdf"))

    }

    seurat_object <- RunHarmony(seurat_object, assay.use = assay.use, group.by.vars)
    seurat_object <- RunUMAP(seurat_object, assay = assay.use, reduction = "harmony", dims = 1:max.dim)
  
    for(i in 1:length(group.by.vars)){
        p1 <- DimPlot(object = seurat_object, reduction = "umap", pt.size = .2, group.by = group.by.vars[i], raster = FALSE, cols = colors) + coord_fixed()
        ggsave(filename=paste0(out_dir,"/umap_after_Harmony_",group.by.vars[i],".pdf"),plot=p1)
    }
    
    return(seurat_object)
}
