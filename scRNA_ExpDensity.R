#' @param data Input Seurat object
#' @param features Genes of interest 

Edensity <- function(data = seu_obj, features = genes){
    
  library(Seurat)
  library(Nebulosa)
  library(ggplot2)
  library(patchwork)

  p = lapply(features,function(z){
              p = plot_density(data, reduction = "umap", z, method = "wkde", 
                          size = 0.2,  raster = T)+ 
              coord_fixed()+theme_void()+
              scale_color_gradientn(colours = c("lightblue","lightyellow","#EF3B36"))
              return(p)
          })
  pp = patchwork::wrap_plots(p, ncol = 4)

}
