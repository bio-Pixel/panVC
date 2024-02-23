############### slingshot ############### 
strace <- function(seu_obj = a){
  library(Seurat)
  library(RColorBrewer)
  library(slingshot) 
  library(SingleCellExperiment) 
  library(ggplot2)
  library(ggthemes)
  seu_obj$idents = as.character(Idents(seu_obj))
  a.sce <- as.SingleCellExperiment(seu_obj) 
  sim <- slingshot(a.sce, clusterLabels = 'celltype', reducedDim = 'UMAP') 
  pdf("trace.pdf",height = 10, width = 10)
  col = brewer.pal(8,'Set1')
  names(col)  = unique(sim$celltype)
  plot(reducedDims(sim)$UMAP, col = col[sim$celltype], pch=16, asp = 1) 
  lines(SlingshotDataSet(sim), lwd=2, type = 'lineages', col = 'black')
  dev.off()
  return(sim)
}

############### Diffusion map ############### 
DCmap <- function(seu_obj = a){
  library(Seurat)
  library(SeuratDisk)
  library(destiny)
  library(SingleCellExperiment)
  library(ggplot2)
  library(ggthemes)
  css = Embeddings(seu_obj,"css")
  dmm <- DiffusionMap(css)
  dpt = DPT(dmm)
  loc = data.frame(eigenvectors(dmm)[,1:3], time = dpt$dpt, type = rownames(css))
  return(list(dpt = dpt, loc = loc, dmm = dmm))
}

############### monocle3 ############### 
mono3 <- function(seu_obj = a){
  library(Seurat)
  library(monocle3)
  library(tidyverse)
  library(patchwork)
  data <- GetAssayData(seu_obj, assay = 'RNA', slot = 'counts')
  cell_metadata <- seu_obj@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  cds <- new_cell_data_set(data,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  cds <- preprocess_cds(cds, num_dim = 50)
  cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA
  
  cds.embed <- cds@int_colData$reducedDims$UMAP
  int.embed <- Embeddings(seu_obj, reduction = "umap")
  int.embed <- int.embed[rownames(cds.embed),]
  cds@int_colData$reducedDims$UMAP <- int.embed
  
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds)
  saveRDS(cds, "monocle3.rds")
  
  p = plot_cells(cds, color_cells_by = "cell_type", label_groups_by_cluster=FALSE,
       label_leaves=FALSE, label_branch_points=FALSE)
}

############### pyVIA data prepare ############### 
## Seurat to anndata
library(Seurat)
library(SeuratDisk)
Seu2Ann = function(seu_obj = a, filename){
    library(Seurat)
    library(SeuratDisk)
    SaveH5Seurat(seu_obj, filename = paste0(filename,".h5Seurat"))
    Convert(paste0(filename,".h5Seurat"), dest = "h5ad")
}

Seu2Ann(seu_obj = a, filename = "Anndata")
write.csv(data.frame(row.names = 1:ncol(a), cell_id = rownames(a@meta.data), group_id = a$celltype), "Anndata_id.csv",quote=F)




