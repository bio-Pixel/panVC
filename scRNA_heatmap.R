#' @param data Input Seurat object
#' @param label Genes of interest labeled in the heatmap

dittoHeatmap <- function(data = seu_obj, label = label){
  library(dittoSeq)
  library(viridis)
  library(ggplot2)
  library(Seurat)
  library(ComplexHeatmap)
  
  dittoHeatmap(data, scaled.to.max = T,
               cluster_row = FALSE, 
               genes = rownames(data), complex = T, 
               raster_quality = 6,
               heatmap.colors.max.scaled = colorRampPalette(c("#FCFDBFFF",rev(magma(10))))(50),
               show_colnames = FALSE)+
    rowAnnotation(label = anno_mark(at = match(label, rownames(data)), labels = label))

}

