#' @param seu_obj Input Seurat object for preclustered mono-cell type dataset 
#' @param Focusedcelltype The cell type of interest in "Epithelial_cell" "Endothelial" "Melanocyte" "T_NK" "B_Plasma" "mast" "Neutrophil" "Myeloid" "Fibroblast" "MC" "Red_blood" 
#' @param canonical.marker Canonical marker genes for the 'Focusedcelltype'

FindContamination <- function(seu_obj, Focusedcelltype, canonical.marker){
    library(Seurat)
    library(ggplot2)
    a = seu_obj
    ################################ filter cells according to scoring rules ########################################
    ############################ calculate the expression proportion of canonical markers ###########################
    canonical.marker <- intersect(canonical.marker, rownames(a))
    p <- DotPlot(a, features = canonical.marker)
    per.ec = p$data
    per.ec = per.ec[which(per.ec$pct.exp >50),]
    clust.high = paste0("cluster",unique(per.ec$id))
    clust.high = data.frame(id = clust.high, annotation = "Focused.pct", score = 1)
    print("Stat clust.high done!")

    ################# whether each cluster's specifically expressed genes include canonical marker genes #############
    mg <- read.delim(paste0(x, "/MarkerGene_2.txt"))
    hg <- apply(mg,2,function(g) intersect(g, canonical.marker))
    clust.marker = names(hg)[which(as.numeric(summary(hg)[,1])>0)]
    clust.marker = data.frame(id = clust.marker, annotation = "Fibro.marker", score = 2)
    print("Stat clust.marker done!")

    ############### Calculate the expression proportion of potentially contaminated cell markers ###############
    marker_gene <- list(
      Epithelial_cell = c("EPCAM"),
      Endothelial = ("PECAM1"),
      Melanocyte = c("PMEL", "MLANA", "TYRP1"),
      T_NK = c("CD2", "CD3D","CD3E","PTPRC","KLRD1","KLRD2","IL2RB", "TRBC2", "TRBC1", "TRAC",'KLRB1','NCR1'),
      B_Plasma = c("CD79A", "MS4A1"),
      mast = c('TPSAB1' , 'TPSB2'),
      Neutrophil = c("CSF3R", "FCGR3B"),
      Myeloid = c('FCGR3A', 'CD68', 'CD163', 'CD14', 'C1QA',  'C1QB','LAMP3', 'IDO1','IDO2', 'FCGR3A', 'CD1E','CD1C'),
      Fibroblast = c("COL1A1"),
      MC = c("RGS5", "PDGFRB", "ACTA2", "TAGLN"),
      Red_blood = c("HBB","HBA1","HBA2")
    )
    black.list.pct <- lapply(setdiff(names(marker_gene),Focusedcelltype), function(cell) {
        g = marker_gene[[cell]]
        g = intersect(g, rownames(a@assays$RNA@counts))
        index <- length(intersect(g, rownames(a@assays$RNA@counts)))
        if (index > 0) {
            p <- DotPlot(a, features = g)
            per <- p$data
            per <- per[which(per$avg.exp.scaled > 0), ]
            per <- per[which(per$pct.exp > 50), ]
            re = data.frame(id = paste0("cluster", unique(per$id)), annotation = paste0(cell,".pct"), score = -0.5)
            return(re)
        }else {
            re = data.frame(id = "cluster", annotation = paste0(cell,".pct"), score = -0.5)
            return(re)
        }
    })
    black.list.pct <- do.call(rbind,black.list.pct)
    black.list.pct <-  black.list.pct[-which(black.list.pct$id == "cluster"),]
    print("Stat black.list.pct done!")
    
    ############### Whether each cluster's specifically expressed genes include other cell markers ############### 
    black.list.hg <- lapply(names(marker_gene), function(cell){
        g = marker_gene[[cell]]
        re = apply(mg,2,function(y) intersect(y,g))
        if(length(re) >0 ){
            re = data.frame(id = names(re)[which(as.numeric(summary(re)[,1])>0)], annotation = paste0(cell,".marker"), score = -1)
            return(re)
        }
    })
    black.list.hg <- do.call(rbind,black.list.hg)
    print("Stat black.list.hg done!")

    black.list <- rbind(black.list.hg, black.list.pct)
    check.beta = paste0(setdiff(names(marker_gene),Focusedcelltype),".marker")
    
    cluster = rbind(clust.high,clust.marker,black.list)
    cluster$beta = 1
    cluster$beta[which(cluster$annotation %in% check.beta)] = 0
    res <- cbind(
        aggregate(cluster$annotation, by = list(id = cluster$id), function(x) {
            paste0(x, collapse = ",")
        }),
        score = aggregate(cluster$score, by = list(id = cluster$id), sum)[,2] * aggregate(cluster$beta, by = list(id = cluster$id), min)[,2]
    )

    ec.qc <- res$id[which(res$score>0)]
    Total = dim(a@meta.data)[1]
    ec.meta = a@meta.data[which(Idents(a) %in% gsub("cluster","",ec.qc)),]
    Remain = dim(ec.meta)[1]
    Delete = Total - Remain
    write.table(data.frame(x, Total, Delete, Remain),
        "Statistic.txt", row.names = F, col.names = F, quote = F, append = T, sep = "\t"
    )
    res$num = table(Idents(a))[gsub("cluster","",res$id)]
    res = data.frame(Set = x, res)
    write.table(res, "Cluster_table.txt", row.names = F, col.names = F, quote = F, append = T, sep = "\t") ## The output should also be checked manually!!!
    aa <- subset(a, idents = gsub("cluster","",ec.qc))
    saveRDS(aa, paste0(a$DataSet[1],".rds"))
    return(aa)
  }
  
