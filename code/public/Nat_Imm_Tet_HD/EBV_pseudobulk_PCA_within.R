### Making a Pseudobulk. Identifying the differential genes from the single cell experiments.
### Remove samples the samples which are an outlier and which you want to remove like DMSO
## Cluster 1 is first cluster Vs Cluster2 is the second cluster which you want to perform the differential
## You can have multiple cluster in both the cluster 1 and cluster 2

pseudobulk_within_cluster_AJ <- function(obj,savedir, group1, group2, grouping_by, cluster="all",
                                         cluster_group = "seurat_clusters", cell_freq=10,
                                         remove_samples=NULL, gene_min_counts = 5, sample_col = "donor",
                                         batch_col = "chip"){
  message("Loading Packages")
  set.seed(123)
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(DESeq2)
  library(edgeR)
  library(RColorBrewer)
  library(Glimma)
  library(ggplot2)
  library(ggrepel)
  library(session)
  library(DEFormats)
  library(EnhancedVolcano)
  library(limma)
  
  if(cluster=="all"){
    print("all")
    savedir2 <- paste(savedir,"pseudobulk/",sep = "")
    dir.create(savedir2, showWarnings = FALSE)
    remove_samples_name <- paste(remove_samples,collapse="_and_")
    savedir3 <- paste(savedir2,"clus_",cluster,"_removed_",remove_samples_name,"_",group1,"_vs_",group2,"/",sep = "")
    dir.create(savedir3,showWarnings = FALSE)
    obj@meta.data[,paste(sample_col,grouping_by,sep = "_")] <- paste(obj@meta.data[,sample_col], obj@meta.data[,grouping_by], sep="_")
    
    DefaultAssay(obj) <- "RNA"
    cell_samples <- table(obj@meta.data[,paste(sample_col,grouping_by,sep = "_")]) %>% as.data.frame()
    sample_low <- cell_samples[cell_samples$Freq < cell_freq,1]
    sample_low <- gsub("_","-",sample_low)
    sample_low_2 <- c(sample_low,remove_samples)
    message(paste("Removing this sample",sample_low_2,"\n",sep = " "))
    Treg_cts <- AggregateExpression(obj, group.by = paste(sample_col,grouping_by,sep = "_"),
                                    assays = 'RNA', slot = "counts", return.seurat = FALSE)
    Treg_cts_2 <- as.data.frame(Treg_cts$RNA)
    if(length(sample_low_2) != 0){
      index <- grep(paste(sample_low_2,collapse="|"),colnames(Treg_cts_2), invert=TRUE)
      Treg_cts_3 <- Treg_cts_2[,index]
    } else {
      Treg_cts_3 <- Treg_cts_2
    }
    
    Treg_metadata <- data.frame(samples = colnames(Treg_cts_3))
    # Function to split and convert to dataframe
    split_to_dataframe <- function(x) {
      split_elements <- strsplit(x, "-")[[1]]
      data.frame(t(split_elements))
    }
    # Convert split elements to dataframe
    # Treg_metadata_2 <- do.call(rbind, lapply(Treg_metadata$samples, split_to_dataframe))
    # colnames(Treg_metadata_2) <- c("Sample","Vaccine",grouping_by)
    # Treg_metadata_2$orig.ident <- Treg_metadata$samples
    # Treg_metadata_2$Run <- gsub("BB22028|BB22029|BB22030|BB22031","Run01",Treg_metadata_2$Sample)
    # Treg_metadata_2$Run <- gsub("BB22032|BB22033|BB22034|BB22035","Run02",Treg_metadata_2$Run)
    # Treg_metadata_2$Run <- gsub("BB22036|IMAG015|IMAG017|IMAG016","Run03",Treg_metadata_2$Run)
    # 
    # stopifnot(all(Treg_metadata_2$orig.ident == colnames(Treg_cts_3))) # if TRUE move forward
    # # Treg_metadata$samplename <- gsub("gE_|DMSO_|_DMSO|_gE|B22015","",Treg_metadata$samples)
    # ## Run will be treated as a covariate in the regression model,
    # 
    # design_AJ <- as.formula(paste("~",paste(colnames(Treg_metadata_2[,c(5,3)]),collapse = "+")))
    # Treg_dds <- DESeqDataSetFromMatrix(countData = Treg_cts_3,
    #                                    colData = Treg_metadata_2,
    #                                    design = design_AJ)
    
    Treg_metadata_2 <- do.call(rbind, lapply(Treg_metadata$samples, split_to_dataframe))
    colnames(Treg_metadata_2) <- c("Disease","samplenumber",grouping_by)
    Treg_metadata_2$orig.ident <- Treg_metadata$samples
    Treg_metadata_2$Run <- Treg_metadata_2$orig.ident
    
    obj@meta.data[,paste(sample_col,batch_col,sep = "_")] <- paste(obj@meta.data[,sample_col],obj@meta.data[,batch_col],sep = "_")
    batch_group <- unique(obj@meta.data[,batch_col])
    for (i in 1:length(batch_group)) {
      donor_change <- unique(obj@meta.data[,paste(sample_col,batch_col,sep = "_")]) %>% grep(batch_group[i],.,value=TRUE) %>% 
        gsub(paste(batch_group[i],sep = ""),"",.) %>% gsub("_","-",.) %>% paste(.,collapse=".*.|")
      donor_change <- paste(donor_change,".*.",sep="")
      Treg_metadata_2$Run <- gsub(donor_change,batch_group[i],Treg_metadata_2$Run)
    }
    Treg_metadata_2$Run <- gsub("-.*.","",Treg_metadata_2$Run)
    
    # Treg_metadata_2$Run <- gsub("T1D-5.*.|T1D-6.*.|T1D-4.*.|healthy-4.*.|healthy-3.*.|healthy-5.*.","Run02",Treg_metadata_2$orig.ident)
    # Treg_metadata_2$Run <- gsub("T1D-2.*.|T1D-1.*.|healthy-2.*.|T1D-3.*.|healthy-1.*.","Run03",Treg_metadata_2$orig.ident)
    # Treg_metadata_2$Run <- gsub("T1D-2.*.|T1D-1.*.|healthy-2.*.|T1D-3.*.|healthy-1.*.","Run04",Treg_metadata_2$orig.ident)
    
    stopifnot(all(Treg_metadata_2$orig.ident == colnames(Treg_cts_3))) # if TRUE move forward
    # Treg_metadata$samplename <- gsub("gE_|DMSO_|_DMSO|_gE|B22015","",Treg_metadata$samples)
    ## Run will be treated as a covariate in the regression model,
    
    design_AJ <- as.formula(paste("~",paste(colnames(Treg_metadata_2[,c("Run",grouping_by)]),collapse = "+")))
    Treg_dds <- DESeqDataSetFromMatrix(countData = Treg_cts_3,
                                       colData = Treg_metadata_2,
                                       design = design_AJ)
    
    keep = filterByExpr(Treg_dds, group=colData(Treg_dds)[,grouping_by], min.count=gene_min_counts)
    table(keep)
    Treg_dds_2 <- Treg_dds[keep,]
    
    rld<-rlog(Treg_dds_2)
    PCA <- DESeq2::plotPCA(object = rld,
                           intgroup = colnames(Treg_metadata_2),
                           returnData=TRUE,ntop = 5000)
    percentVar <- round(100 * attr(PCA, "percentVar"))
    
    savedir4 <- paste(savedir3,"PCA/",sep = "")
    dir.create(savedir4, showWarnings = FALSE)
    
    PC1 <-ggplot(PCA, aes_string("PC1", "PC2", label="name", color="Run", shape=grouping_by)) +
      geom_point(size=5, alpha=0.7) +
      geom_text_repel(size=2) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      ggtitle("PCA") +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    pdf(paste(savedir4,"PCA_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".pdf",sep = ""),
        width = 9, height = 8)
    print(PC1)
    dev.off()
    
    assay(rld) <- limma::removeBatchEffect(assay(rld),
                                           batch=rld$Run)
    
    PCA_batch <- DESeq2::plotPCA(object = rld,
                                 intgroup = colnames(Treg_metadata_2),
                                 returnData=TRUE,ntop = 5000)
    
    percentVar <- round(100 * attr(PCA_batch, "percentVar"))
    
    PC1 <-ggplot(PCA_batch, aes_string("PC1", "PC2", label="name", color=grouping_by, shape=grouping_by)) +
      geom_point(size=5, alpha=0.7) +
      geom_text_repel(size=2) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    pdf(paste(savedir4,"PCA_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".pdf",sep = ""),
        width = 9, height = 8)
    print(PC1)
    dev.off()
    
    message("Please find the PCA at this location ",savedir4)
    return(Treg_dds_2)
  } else {
    print("cluster")
    cluster_name <- paste(cluster,collapse = "_")
    savedir2 <- paste(savedir,"pseudobulk/",sep = "")
    dir.create(savedir2, showWarnings = FALSE)
    remove_samples_name <- paste(remove_samples,collapse="_and_")
    savedir3 <- paste(savedir2,"clus_",cluster_name,"_removed_",remove_samples_name,"_",group1,"_vs_",group2,"/",sep = "")
    dir.create(savedir3,showWarnings = FALSE)
    cluster_num <- paste("^",cluster,"$",sep = "")
    if (cluster_group == "seurat_clusters") {
      obj@meta.data$merge_cluster <- as.integer(gsub(paste(cluster_num, collapse="|"),100,obj@meta.data[,cluster_group]))
    } else{
      obj@meta.data$merge_cluster <- gsub(paste(cluster_num, collapse="|"),100,obj@meta.data[,cluster_group])
    }
    
    Idents(obj) <- obj@meta.data$merge_cluster
    obj_subset <- subset(obj,idents = c(100))
    obj_subset@meta.data[,paste(sample_col,grouping_by,sep = "_")] <- paste(obj_subset@meta.data[,sample_col], obj_subset@meta.data[,grouping_by], sep="_")
    
    DefaultAssay(obj_subset) <- "RNA"
    cell_samples <- table(obj_subset@meta.data[,paste(sample_col,grouping_by,sep = "_")]) %>% as.data.frame()
    sample_low <- cell_samples[cell_samples$Freq < cell_freq,1]
    sample_low <- gsub("_","-",sample_low)
    sample_low_2 <- c(sample_low,remove_samples)
    message(paste("Removing this sample",sample_low_2,"\n",sep = " "))
    Treg_cts <- AggregateExpression(obj_subset, group.by = paste(sample_col,grouping_by,sep = "_"),
                                    assays = 'RNA', slot = "counts", return.seurat = FALSE)
    Treg_cts_2 <- as.data.frame(Treg_cts$RNA)
    if(length(sample_low_2) != 0){
      index <- grep(paste(sample_low_2,collapse="|"),colnames(Treg_cts_2), invert=TRUE)
      Treg_cts_3 <- Treg_cts_2[,index]
    } else {
      Treg_cts_3 <- Treg_cts_2
    }
    
    Treg_metadata <- data.frame(samples = colnames(Treg_cts_3))
    # Function to split and convert to dataframe
    split_to_dataframe <- function(x) {
      split_elements <- strsplit(x, "-")[[1]]
      data.frame(t(split_elements))
    }
    
    # Convert split elements to dataframe
    Treg_metadata_2 <- do.call(rbind, lapply(Treg_metadata$samples, split_to_dataframe))
    colnames(Treg_metadata_2) <- c("Disease","samplenumber",grouping_by)
    Treg_metadata_2$orig.ident <- Treg_metadata$samples
    Treg_metadata_2$Run <- Treg_metadata_2$orig.ident
    
    obj@meta.data[,paste(sample_col,batch_col,sep = "_")] <- paste(obj@meta.data[,sample_col],obj@meta.data[,batch_col],sep = "_")
    batch_group <- unique(obj@meta.data[,batch_col])
    for (i in 1:length(batch_group)) {
      donor_change <- unique(obj@meta.data[,paste(sample_col,batch_col,sep = "_")]) %>% grep(batch_group[i],.,value=TRUE) %>% 
        gsub(paste(batch_group[i],sep = ""),"",.) %>% gsub("_","-",.) %>% paste(.,collapse=".*.|")
      donor_change <- paste(donor_change,".*.",sep="")
      Treg_metadata_2$Run <- gsub(donor_change,batch_group[i],Treg_metadata_2$Run)
    }
    Treg_metadata_2$Run <- gsub("-.*.","",Treg_metadata_2$Run)
    
    stopifnot(all(Treg_metadata_2$orig.ident == colnames(Treg_cts_3))) # if TRUE move forward
    # Treg_metadata$samplename <- gsub("gE_|DMSO_|_DMSO|_gE|B22015","",Treg_metadata$samples)
    ## Run will be treated as a covariate in the regression model,
    
    design_AJ <- as.formula(paste("~",paste(colnames(Treg_metadata_2[,c("Run",grouping_by)]),collapse = "+")))
    Treg_dds <- DESeqDataSetFromMatrix(countData = Treg_cts_3,
                                       colData = Treg_metadata_2,
                                       design = design_AJ)
    
    keep = filterByExpr(Treg_dds, group=colData(Treg_dds)[,grouping_by], min.count=gene_min_counts)
    table(keep)
    Treg_dds_2 <- Treg_dds[keep,]
    
    rld<-rlog(Treg_dds_2)
    PCA <- DESeq2::plotPCA(object = rld,
                           intgroup = colnames(Treg_metadata_2),
                           returnData=TRUE,ntop = 5000)
    percentVar <- round(100 * attr(PCA, "percentVar"))
    
    savedir4 <- paste(savedir3,"PCA/",sep = "")
    dir.create(savedir4, showWarnings = FALSE)
    
    PC1 <-ggplot(PCA, aes_string("PC1", "PC2", label="name", color="Run", shape=grouping_by)) +
      geom_point(size=5, alpha=0.7) +
      geom_text_repel(size=2) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      ggtitle("PCA") +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    pdf(paste(savedir4,"PCA_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".pdf",sep = ""),
        width = 9, height = 8)
    print(PC1)
    dev.off()
    
    assay(rld) <- limma::removeBatchEffect(assay(rld),
                                           batch=rld$Run)
    
    PCA_batch <- DESeq2::plotPCA(object = rld,
                                 intgroup = colnames(Treg_metadata_2),
                                 returnData=TRUE,ntop = 5000)
    
    percentVar <- round(100 * attr(PCA_batch, "percentVar"))
    
    PC1 <-ggplot(PCA_batch, aes_string("PC1", "PC2", label="name", color=grouping_by, shape="Run")) +
      geom_point(size=5, alpha=0.7) +
      geom_text_repel(size=2) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    pdf(paste(savedir4,"PCA_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".pdf",sep = ""),
        width = 9, height = 8)
    print(PC1)
    dev.off()
    
    message("Please find the PCA at this location ",savedir4)
    return(Treg_dds_2)
  }
}


