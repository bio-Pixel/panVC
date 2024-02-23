#' @TODO Cluster the cells and return marker genes
#' @param seurat_object Input Seurat object
#' @param FindNeighbor_reduction Calculate the dimensions used for computing the nearest-neighbor graph among cells
#' @param resolution Resolution used for clustering
#' @param cores Number of threads for parallel computation
#' @param top_n Number of top genes in the output
#' @param out_dir Output marker genes and the path of images
#' @returnType Seurat object
#' @return 
#' 
#' author LX
#' 

Seurat_marker <- function(seurat_object = NULL, FindNeighbor_reduction = NULL, resolution = resolution, cores = NULL, top_n = NULL, out_dir = NULL){
    library(Seurat)
    library(ggplot2)
    library(ggthemes)
    library(dplyr)
    library(doMC)
    library(openxlsx)
    #library(xlsx)

    registerDoMC(cores)

    seurat_object <- FindNeighbors(seurat_object, reduction = FindNeighbor_reduction)
    seurat_object <- FindClusters(seurat_object, resolution = resolution)
    p1 <- DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = .2) + coord_fixed()
    ggsave(p1, file = paste0(out_dir, "/umap_resolution_",resolution,".pdf"))

    p1 <- DimPlot(seurat_object, reduction = "tsne", label = TRUE, pt.size = .2) + coord_fixed()
    ggsave(p1, file = paste0(out_dir, "/tsne_resolution_",resolution,".pdf"))

    p1 <- DimPlot(seurat_object, reduction = "umap", split.by="SampleID", label = TRUE, pt.size = .2) + coord_fixed()
    ggsave(p1, file = paste0(out_dir, "/umap_SampleID_resolution_",resolution,".pdf"))

    p1 <- DimPlot(seurat_object, reduction = "tsne", split.by="SampleID", label = TRUE, pt.size = .2) + coord_fixed()
    ggsave(p1, file = paste0(out_dir, "/tsne_SampleID_resolution_",resolution,".pdf"))
    #calculate marker genes
    markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    markers_sig <- markers[which(markers$avg_log2FC > 0.25 & markers$p_val_adj < 0.01),]
    openxlsx::write.xlsx(markers_sig, file= paste0(out_dir, "/MarkerGene_",resolution,"_sig.xlsx"), overwrite = TRUE)

    markers <- markers_sig %>% group_by(cluster) %>% top_n(n = top_n, wt = avg_log2FC)
    
    top_genes <- markers %>% select(gene,cluster)
    top_genes$cluster <- paste0('cluster',top_genes$cluster)
    #Select clusters with fewer than top_n.
    clusters <- names(table(top_genes$cluster)[which(table(top_genes$cluster) < top_n)])
    if(length(clusters) > 0){
    for(i in 1:length(clusters)){
        this.supple <- data.frame(gene = rep(NA, top_n - as.numeric(table(top_genes$cluster)[clusters[i]])),
        cluster = rep(clusters[i], top_n - as.numeric(table(top_genes$cluster)[clusters[i]])), stringsAsFactors = FALSE)
        top_genes <- rbind(top_genes, this.supple)
    }
    }

    #Select clusters with more than top_n
    clusters <- names(table(top_genes$cluster)[which(table(top_genes$cluster) > top_n)])
    if(length(clusters) > 0){
    for(i in 1:length(clusters)){
        this.cluster <- top_genes[which(top_genes$cluster == clusters[i]),]
        this.cluster <- this.cluster[1:top_n,]
        top_genes <- top_genes[which(top_genes$cluster != clusters[i]),]
        top_genes <- rbind(top_genes, this.cluster)
    }
    }

    unstack(top_genes) %>% write.csv(paste0(out_dir, "/MarkerGene_",resolution,".csv"),row.names=F)
    unstack(top_genes) %>% openxlsx::write.xlsx(file= paste0(out_dir, "/MarkerGene_",resolution,".xlsx"), overwrite = TRUE)
    return(seurat_object)
}
