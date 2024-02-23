#' @param seu_st_obj Input Seurat object
#' @param bk.dat bulk expression matrix 
#' @param features cell type signature genes

RunBayesPrism <- function(seu_obj, bk.dat, features){
  library(Seurat)
  library(dplyr)
  library(BayesPrism)
  suba <- subset(seu_obj, features = features)
  sc.dat.filtered.pc.sig <- t(as.matrix(suba@assays$RNA@counts))[,features]
  cell.type.labels <- as.character(Idents(suba))
  cell.state.labels <- cell.type.labels
  
  ###################### Construct a prism object ###################### 
  myPrism <- new.prism(reference=sc.dat.filtered.pc.sig,
                       mixture=bk.dat, 
                       input.type="count.matrix", 
                       cell.type.labels = cell.type.labels, 
                       cell.state.labels = cell.state.labels, 
                       key = NULL,
                       outlier.cut=0.01,
                       outlier.fraction=0.1)
  bp.res <- run.prism(prism = myPrism, n.cores=10)
  theta <- get.fraction(bp=bp.res, which.theta="final",
                        state.or.type="type")
  
  # extract coefficient of variation (CV) of cell type fraction
  theta.cv <- bp.res@posterior.theta_f@theta.cv
  re = list(theta = theta, theta.cv = theta.cv)
  return(re)
}

