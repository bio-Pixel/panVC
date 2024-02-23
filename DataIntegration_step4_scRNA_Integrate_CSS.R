#' @param dir.path Path to the folder containing the processed Seurat objects of single-cell type

Integrate_CSS <- function(dir.path){
  library(Seurat)
  library(dplyr)
  library(simspec) 
  seu = lapply(list.files(dir.path), function(file){
    a <- readRDS(file)
    return(a)
  })
  seu = merge(seu[[1]], seu[-1])
  source("scRNA_CSS.R")
  seu_css = CSS_scaling(data = seu)
  return(seu_css)
}
