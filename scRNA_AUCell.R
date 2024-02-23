#' @param data Input Seurat object
#' @param geneSets Gene sets of interest 

aucell <- function(data, geneSets = NULL){

  cells_ranking = AUCell_buildRankings(data@assays$RNA@data)
  if(length(geneSets) > 0){
    if(class(geneSets) == "data.frame"){
      geneSets = apply(geneSets, 2, function(x){
        x = c(x,"")
        x[-which(x == "")]
      })
    }
  }
  cells_AUC <- AUCell_calcAUC(geneSets, cells_ranking, aucMaxRank=nrow(cells_ranking)*0.1)
  aucs = getAUC(cells_AUC)
  aucs = data.frame(t(aucs))
  cells_AUC <- cbind(data@meta.data, aucs)
  return(cells_AUC)
}

