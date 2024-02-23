#' @param seu_st_obj Input Seurat object
#' @param sc_count single cell count matrix 
#' @param sc_meta single cell meta file including column "cellType" 

RunCARD <- function(seu_st_obj, sc_count, sc_meta){

  library(Seurat)
  library(CARD)
  library(ggplot2)
  library(patchwork)
  
  spatial_count <- seu_st_obj@assays$Spatial@counts
  spatial_loca <- seu_st_obj@images$slice1@coordinates
  spatial_location <- spatial_loca[,2:3]
  colnames(spatial_location) <- c("x","y")
  CARD_obj = createCARDObject(
                sc_count = sc_count,
                sc_meta = sc_meta,
                spatial_count = spatial_count,
                spatial_location = spatial_location,
                ct.varname = "cellType",
                ct.select = unique(sc_meta$cellType),
                sample.varname = "sampleInfo",
                minCountGene = 100,
                minCountSpot = 5
  ) 
  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
  CARD_obj = CARD.imputation(CARD_obj,NumGrids = 2500,ineibor = 10,exclude = NULL)
  re = cbind(CARD_obj@spatial_location, CARD_obj@Proportion_CARD)
  return(list(re, CARD_obj))

}

