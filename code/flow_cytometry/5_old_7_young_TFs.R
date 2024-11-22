#### 50k Cells ######
### CD8 lytic clustering with Seurat and Phenograph 
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.

rawData_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/5_old_7_young/TFs/raw_data/"
bulk_csv_files = list.files(rawData_path,pattern = ".csv")
dir.create("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/5_old_7_young/TFs/analysis",showWarnings = FALSE)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/5_old_7_young/TFs/analysis/"

column = c("APC.eFluor.780.A....Eomes","Alexa.Fluor.488.A....BLIMP","Alexa.Fluor.700.A....RUNX3","BUV496.A....CCR7","BUV661.A....BCL6",
           "BUV737.A....CD28","BUV805.A....CD27","BV421.A....TCF1","BV605.A....CD45RA","BV650.A....CD45RO","BV785.A....Tbet",
           "Pacific.Blue.A....Helios","alexa.fluor.594.A....Tox.Tox2")

sample_id = gsub(".*.TF_|_2022_.*","",bulk_csv_files)
Run_id = c("Exp01","Exp01","Exp02","Exp02","Exp02","Exp02","Exp02","Exp02","Exp03","Exp03","Exp03","Exp03")
age = c("young","old","young","young","young","old","old","old","young","old","young","young")
name=paste(age,sample_id,Run_id,sep = "_")

for (i in 1:length(bulk_csv_files)) {
  bulk_sample = read.csv(paste(rawData_path,bulk_csv_files[i],sep = ""))
  bulk_sample = bulk_sample[,which(colnames(bulk_sample) %in% column)]
  colnames(bulk_sample) = gsub("FJComp.","",colnames(bulk_sample))
  rownames(bulk_sample) = paste("cell",1:nrow(bulk_sample),sep = "")
  bulk_sample_obj = CreateSeuratObject(t(bulk_sample), project = name[i], assay = "ADT")
  obj_name= paste(name[i],"obj",sep = "_")
  assign(obj_name,bulk_sample_obj)
}

### merging the 7 young datasets
combine <- merge(young_40187_Exp01_obj, y = c(old_40188_Exp01_obj,young_40189_Exp02_obj,young_40191_Exp02_obj,young_40193_Exp02_obj,old_40197_Exp02_obj,
                                              old_40198_Exp02_obj,old_40199_Exp02_obj,young_40222_Exp03_obj,old_40223_Exp03_obj,young_40224_Exp03_obj,young_40225_Exp03_obj), 
                 add.cell.ids = name, project = "Combined")

combine@meta.data$age = gsub("_.*","",combine@meta.data$orig.ident)
combine@meta.data$Run = gsub(".*_","",combine@meta.data$orig.ident)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")
objname = "CD8"
process = "integrated"
Assay = "ADT"

CD8_ADT_integrated <- ADT_merging(combine, savedir, dims = 10, numfeatures=13,
                                  Assay=Assay, process=process, objname=objname, 
                                  sample_tree = NULL, split_by = "Run",
                                  reference=NULL)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/scRNA_integration.R")
objname = "CD8_TF"
process = "integration"
Assay = "ADT"
CD8_ADT_integrated <- RNA_integration(obj_path = CD8_ADT_integrated, saveDir = savedir, dims = 10,
                                      RNA_features = c("Alexa.Fluor.488.A....BLIMP","BUV805.A....CD27"),Assay = Assay,
                                      process = process, objname = objname, ncol = 3)


library(ggplot2)
library(ArchR)
rm(plot_list)
DefaultAssay(CD8_ADT_integrated) = "ADT"
plot_list = list()
for (i in 1:length(rownames(CD8_ADT_integrated))) {
  p = FeaturePlot(CD8_ADT_integrated,rownames(CD8_ADT_integrated)[i], reduction = "umap") + scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] = p
}


dir.create(paste(savedir,"featureplot",sep = ""), showWarnings = FALSE)
require(gridExtra)
pdf(paste(savedir,"featureplot/featureplot_ADT_SolarExtra_defaultAssay_ADT.pdf",sep = ""), width = 20, height = 20)
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
             plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
             plot_list[[13]], ncol = 4, nrow = 4)
dev.off()

DefaultAssay(CD8_ADT_integrated) = "MAGIC_ADT"
rm(plot_list)
plot_list = list()
for (i in 1:length(rownames(CD8_ADT_integrated))) {
  p = FeaturePlot(CD8_ADT_integrated,rownames(CD8_ADT_integrated)[i], reduction = "umap") + 
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] = p
}

require(gridExtra)
pdf(paste(savedir,"featureplot/featureplot_MAGIC_ADT.pdf",sep = ""), width = 20, height = 20)
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
             plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
             plot_list[[13]], ncol = 4, nrow = 4)
dev.off()

### Nearest Neighbour and Clustering
DefaultAssay(CD8_ADT_integrated) <- "integrated"
CD8_ADT_integrated_NN <- FindNeighbors(CD8_ADT_integrated, dims = 1:10)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/scRNA_cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8)
for (i in 1:length(res)) {
  process = paste("ADT_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  CD8_ADT_integrated_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD8_ADT_integrated_NN, dims = 10, res = res[i], saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_ADT","nFeature_ADT"), objname = "CD8_ADT_integrated_NN_cluster",
                                                       process = process)
}

dir.create(paste(savedir,"saveRDS_obj/",sep = ""),showWarnings = FALSE)
saveRDS(CD8_ADT_integrated_NN_cluster,paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster_0.8.RDS",sep = ""))

## Making a heatmap
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/Heatmap.R")
AJ_heatmap(obj = CD8_ADT_integrated_NN_cluster, savedir = savedir)

### Separating out all the clusters ###
reference=CD8_ADT_integrated_NN_cluster
for (i in 1:length(unique(reference@meta.data$seurat_clusters))) {
  j=i-1
  reference_subset=subset(reference, subset = seurat_clusters==j)
  p=DimPlot(reference_subset)
  
  pdf(paste(savedir,"UMAP/Bulk_cluster",j,".pdf",sep = ""), width = 6.5, height = 6)
  print(p)
  dev.off()
}

### Easier and better way then down
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/5_old_7_young/TFs/analysis/"

p <- DimPlot(reference, reduction = "umap", label = TRUE, label.size = 5,
             cols = c( '0' = "khaki1", '1' = "yellow1", '2' = "goldenrod2", '3' = "forestgreen", '4' = "tomato", 
                       "5"="steelblue1", "6"="firebrick3", "7"="seagreen3", "8" ="gold4","9"="yellow3", 
                       "10"="steelblue1", "11"="red4", "12"="blue1","13"="firebrick1","14"="darkgoldenrod3",
                       "15"="gold","16"="royalblue","17"="limegreen","18"="gold")) + NoLegend()

pdf(paste(savedir,"UMAP/bulk_mapping.pdf",sep = ""), width = 6, height = 6)
p
dev.off()

#### 50k Query and Reference Mapping #####
reference = CD8_ADT_integrated_NN_cluster
reference = RunUMAP(reference, dims = 1:10, reduction = "pca", return.model = TRUE)
saveRDS(reference,paste(savedir,"saveRDS_obj/reference.RDS",sep = ""))
# reference=readRDS(paste(savedir,"saveRDS_obj/reference.RDS",sep = ""))

tetramer_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_panel_tetramer/"
teteramer_files = list.files(tetramer_path, pattern = "csv",full.names = TRUE)
tetramer_basename = basename(teteramer_files)
tetramer_name = gsub(".*TF_|_2022_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]|.csv","",tetramer_basename)

for (i in 1:length(teteramer_files)) {
  query_csv <- read.csv(teteramer_files[i])
  query_csv = query_csv[,which(colnames(query_csv) %in% column)]
  # colnames(query_csv) = gsub("FJComp.","",colnames(query_csv))
  rownames(query_csv) = paste("query",tetramer_name[i],1:nrow(query_csv),sep = "_")
  query_obj = CreateSeuratObject(t(query_csv), project = tetramer_name[i], assay = "ADT")
  query_obj <- NormalizeData(query_obj, normalization.method = 'CLR', margin = 2)
  query_obj <- FindVariableFeatures(query_obj, selection.method = "vst", nfeatures = 13)
  query.anchors <- FindTransferAnchors(reference = reference, query = query_obj,dims = 1:10, reference.reduction = "pca")
  predictions <- TransferData(anchorset = query.anchors, refdata = reference$seurat_clusters,dims = 1:10)
  query_obj <- AddMetaData(query_obj, metadata = predictions)
  query_obj <- MapQuery(anchorset = query.anchors, reference = reference, query = query_obj,
                        refdata = list(seurat_clusters = "seurat_clusters"), reference.reduction = "pca", reduction.model = "umap")
  
  dir.create(paste(savedir,"saveRDS_obj",sep = ""),showWarnings = FALSE)
  saveRDS(query_obj, paste(savedir,"saveRDS_obj/",tetramer_name[i],"_50k_bulk_reference_mapped_obj.RDS",sep = ""))
  
  p1 <- DimPlot(reference, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3,
                repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
  p2 <- DimPlot(query_obj, reduction = "ref.umap", group.by = "predicted.seurat_clusters", label = TRUE,
                label.size = 3, repel = TRUE) + NoLegend() + ggtitle(paste(tetramer_name[i],"reference mapped",sep = " "))
  
  pdf(paste(savedir,"UMAP/",tetramer_name[i],"_reference_mapping_bulk_50k.pdf",sep = ""),width = 12, height = 6)
  print(p1 + p2)
  dev.off()
}

### Downstream files #####
# here we are recoloring the UMAP with the desired color as in feature plot it will provide the different colors
# Also we are saving the number of cells mapped to the reference 
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/5_old_7_young/TFs/analysis/"
mapped_files=list.files(paste(savedir,"saveRDS_obj/",sep = ""),pattern = "_mapped_obj.RDS",full.names = TRUE)
basename=gsub("_50k_bulk_reference_mapped_obj.RDS","",basename(mapped_files))
dir.create(paste(savedir,"reference_mapping",sep = ""),showWarnings = FALSE)
dir.create(paste(savedir,"reference_mapping/table",sep = ""),showWarnings = FALSE)
dir.create(paste(savedir,"reference_mapping/UMAP",sep = ""),showWarnings = FALSE)

for (i in 1:length(basename)) {
  obj=readRDS(mapped_files[i])
  mapped_cluster=as.data.frame(table(obj@meta.data$predicted.seurat_clusters))
  mapped_cluster_order=mapped_cluster[order(-mapped_cluster$Freq),]
  colnames(mapped_cluster_order) = c("Cluster","Number_of_Cells_mapped")
  write.table(mapped_cluster_order,paste(savedir,"reference_mapping/table/",basename[i],".txt",sep = ""), 
              quote = FALSE,row.names = FALSE, col.names = TRUE, sep = "\t")
  
  if (all(rownames(obj@reductions$ref.umap@cell.embeddings) == rownames(obj@meta.data))) {
    UMAP=as.data.frame(obj@reductions$ref.umap@cell.embeddings)
    UMAP$predicted_seurat_cluster=obj@meta.data$predicted.seurat_clusters
    UMAP$colors <- gsub("^0$","khaki1",
                        gsub("^1$","yellow1",
                             gsub("^2$","goldenrod2",
                                  gsub("^3$","forestgreen",
                                       gsub("^4$","tomato",
                                            gsub("^5$","steelblue1",
                                                 gsub("^6$","firebrick3",
                                                      gsub("^7$","seagreen3",
                                                           gsub("^8$","gold4",
                                                                gsub("^9$","yellow3",
                                                                     gsub("^10$","steelblue1",
                                                                          gsub("^11$","red4",
                                                                               gsub("^12$","blue1",
                                                                                    gsub("^13$","firebrick1",
                                                                                         gsub("^14$","darkgoldenrod3",
                                                                                              gsub("^15$","gold",
                                                                                                   gsub("^16$","royalblue",
                                                                                                        gsub("^17$","limegreen",
                                                                                                             gsub("^18$","gold",
                                                                                                        UMAP$predicted_seurat_cluster)))))))))))))))))))
    
    
    col <- as.character(UMAP$colors)
    names(col) <- as.character(UMAP$predicted_seurat_cluster)
    
    library(dplyr)
    label.df_2 <- UMAP %>% 
      group_by(predicted_seurat_cluster) %>% 
      summarize(refUMAP_1 = mean(refUMAP_1), refUMAP_2 = mean(refUMAP_2))
    
    p1 <- ggplot(UMAP,aes(x=refUMAP_1,y=refUMAP_2,label=predicted_seurat_cluster)) + 
      geom_point(aes(colour = predicted_seurat_cluster), show.legend = TRUE) +
      scale_color_manual(values=col) + NoLegend() + 
      ggrepel::geom_label_repel(data = label.df_2, aes(label = predicted_seurat_cluster), size=5) +
      ggtitle(paste(basename[i]," reference mapped")) +
      theme(panel.background = element_rect(fill = 'white',color='black'),
            plot.title = element_text(hjust = 0.5))
    
    pdf(paste(savedir,"reference_mapping/UMAP/",basename[i],"_reference_mapping_bulk_50k_specific_colors.pdf",sep = ""),width = 6, height = 6)
    print(p1)
    dev.off()
  }
}


#### Making a bubble plot
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/5_old_7_young/TFs/analysis/"
CD8_C_S <- read.table(paste(maindir,"Table/CD8_ADT_integrated_NN_cluster_integrated_ADT_UMAP_QC_0.8_cluster_and_samples.txt",
                            sep = ""), header = T)
colnames(CD8_C_S) <- c(paste("C",0:(length(CD8_C_S)-1),sep = ""))
# colnames(CD8_C_S) <- c(paste("C",0:(length(unique(CD8_modality_integrate_cluster@meta.data$seurat_clusters))-1),sep = ""))
CD8_C_S[is.na(CD8_C_S)] <- 0
CD8_C_S_sum <- as.vector(rowSums(CD8_C_S))

### Making an empty dataframe 
df <- data.frame(matrix(nrow = nrow(CD8_C_S), ncol = ncol(CD8_C_S)))
rownames(df) <- rownames(CD8_C_S)
colnames(df) <- colnames(CD8_C_S)

for (i in 1:nrow(CD8_C_S)) {
  df[i,] <- (CD8_C_S[i,]/CD8_C_S_sum[i])*100
}

df$Samples <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("Samples","Cluster","Cell_percentage")

df_melted$Sampleid = gsub("old_|young_","",df_melted$Samples)
df_melted$Age = gsub("_.*","",df_melted$Samples)
df_melted$Sampleid = factor(df_melted$Sampleid,levels = unique(df_melted$Sampleid))
df_melted$Age = factor(df_melted$Age,levels = unique(df_melted$Age))

p <- ggplot(data=df_melted, aes(x=Cluster, y=Sampleid, size=Cell_percentage, color = Age))

pdf(paste(maindir,"Table/CD8_bubble_plot_sample.pdf",sep = ""),width = 9, height = 8)
p + geom_point(alpha=.75) + 
  scale_size(range = c(1,15), breaks=seq(0,40,by=15)) +  
  theme_bw() + ggtitle("CD8_Cluster_Sample")
dev.off()

#### Colwise sum
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/5_old_7_young/TFs/analysis/"
CD8_C_S <- read.table(paste(maindir,"Table/CD8_ADT_integrated_NN_cluster_integrated_ADT_UMAP_QC_0.8_cluster_and_samples.txt",
                            sep = ""), header = T)
colnames(CD8_C_S) <- c(paste("C",0:(length(CD8_C_S)-1),sep = ""))
# colnames(CD8_C_S) <- c(paste("C",0:(length(unique(CD8_modality_integrate_cluster@meta.data$seurat_clusters))-1),sep = ""))
CD8_C_S[is.na(CD8_C_S)] <- 0
CD8_C_S_sum <- as.vector(colSums(CD8_C_S))

### Making an empty dataframe 
df <- data.frame(matrix(nrow = nrow(CD8_C_S), ncol = ncol(CD8_C_S)))
rownames(df) <- rownames(CD8_C_S)
colnames(df) <- colnames(CD8_C_S)

for (i in 1:ncol(CD8_C_S)) {
  df[,i] <- (CD8_C_S[,i]/CD8_C_S_sum[i])*100
}

df$Samples <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("Samples","Cluster","Cell_percentage")

df_melted$Sampleid = gsub("old_|young_","",df_melted$Samples)
df_melted$Age = gsub("_.*","",df_melted$Samples)
df_melted$Sampleid = factor(df_melted$Sampleid,levels = unique(df_melted$Sampleid))
df_melted$Age = factor(df_melted$Age,levels = unique(df_melted$Age))


p <- ggplot(data=df_melted, aes(x=Cluster, y=Samples, size=Cell_percentage, color = Age))

pdf(paste(maindir,"Table/CD8_bubble_plot_ColSums.pdf",sep = ""),width = 9, height = 8)
p + geom_point(alpha=.75) + 
  scale_size(range = c(1,15), breaks=seq(0,100,by=15)) +  
  theme_bw() + ggtitle("CD8 modality Integrated ColSums")
dev.off()

### Antigenic Specific #####
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/5_old_7_young/TFs/analysis/"
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/5_old_7_young/TFs/analysis/reference_mapping/"
saveRDS = list.files(paste(maindir,"saveRDS_obj/",sep = ""),pattern = "_50k_bulk_reference_mapped_obj.RDS",full.names = TRUE)
filename = gsub("_50k_bulk_reference_mapped_obj.RDS","",basename(saveRDS))

for (i in 1:length(saveRDS)) {
  A40225_lyt2 = readRDS(saveRDS[i])
  df = as.data.frame(table(A40225_lyt2@meta.data$predicted.seurat_clusters))
  colnames(df) = c("Cluster","Cellnumber")
  cellsum=sum(df$Cellnumber)
  df$percentage_cell = (df$Cellnumber/cellsum)*100
  write.table(df,paste(savedir,"table/",filename[i],"cells_cluster.txt",sep = ""),col.names = T, row.names = F, sep = "\t",quote = FALSE)
  
  p <- ggplot(df, aes(x=Cluster, y=percentage_cell)) +
    geom_bar(stat = "identity",color="black", fill="grey", position = "dodge") +
    theme_bw() + ggtitle(filename[i])
  
  dir.create(paste(savedir,"graphs",sep =""),showWarnings = FALSE)
  pdf(paste(savedir,"graphs/",filename[i],"_bargraph.pdf",sep = ""), width = 8, height = 6)
  print(p)
  dev.off()
}

library(ggplot2)
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/5_old_7_young/TFs/analysis/"
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/5_old_7_young/TFs/analysis/reference_mapping/UMAP/500_cells/"

dir.create(savedir,showWarnings = FALSE)
reference=readRDS(paste(maindir,"saveRDS_obj/reference.RDS",sep = ""))
tetramer_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_panel_tetramer/500_cells/"
teteramer_files = list.files(tetramer_path, pattern = "500_cell",full.names = TRUE)

## Since we are only interested in the young and upsampled so we will extract those out
teteramer_files_upsample=teteramer_files[grep("upsample",basename(teteramer_files))]

tetramer_basename = basename(teteramer_files_upsample)
tetramer_name = gsub(".csv","",tetramer_basename)

column = c("APC.eFluor.780.A....Eomes","Alexa.Fluor.488.A....BLIMP","Alexa.Fluor.700.A....RUNX3","BUV496.A....CCR7","BUV661.A....BCL6",
           "BUV737.A....CD28","BUV805.A....CD27","BV421.A....TCF1","BV605.A....CD45RA","BV650.A....CD45RO","BV785.A....Tbet",
           "Pacific.Blue.A....Helios","alexa.fluor.594.A....Tox.Tox2")

for (i in 1:length(tetramer_name)) {
  query_csv <- read.csv(teteramer_files_upsample[i])
  # colnames(query_csv) = gsub("FJComp.","",colnames(query_csv))
  query_csv = query_csv[,which(colnames(query_csv) %in% column)]
  rownames(query_csv) = paste("query",tetramer_name[i],1:nrow(query_csv),sep = "_")
  query_obj = CreateSeuratObject(t(query_csv), project = tetramer_name[i], assay = "ADT")
  query_obj <- NormalizeData(query_obj, normalization.method = 'CLR', margin = 2)
  query_obj <- FindVariableFeatures(query_obj, selection.method = "vst", nfeatures = 18)
  query.anchors <- FindTransferAnchors(reference = reference, query = query_obj,dims = 1:10, reference.reduction = "pca")
  predictions <- TransferData(anchorset = query.anchors, refdata = reference$seurat_clusters,dims = 1:10)
  query_obj <- AddMetaData(query_obj, metadata = predictions)
  query_obj <- MapQuery(anchorset = query.anchors, reference = reference, query = query_obj,
                        refdata = list(seurat_clusters = "seurat_clusters"), reference.reduction = "pca", reduction.model = "umap")
  
  dir.create(paste(savedir,"saveRDS_obj",sep = ""),showWarnings = FALSE)
  saveRDS(query_obj, paste(savedir,"saveRDS_obj/",tetramer_name[i],"_50k_bulk_reference_mapped_obj.RDS",sep = ""))
  
  p1 <- DimPlot(reference, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3,
                repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
  p2 <- DimPlot(query_obj, reduction = "ref.umap", group.by = "predicted.seurat_clusters", label = TRUE,
                label.size = 3, repel = TRUE) + NoLegend() + ggtitle(paste(tetramer_name[i],"reference mapped",sep = " "))
  
  pdf(paste(savedir,tetramer_name[i],"_reference_mapping_bulk_50k.pdf",sep = ""),width = 12, height = 6)
  print(p1 + p2)
  dev.off()
}

### Downsampling ####
ref_map=list.files(paste(maindir,"saveRDS_obj",sep = ""),pattern = "_50k_bulk_reference_mapped_obj.RDS",full.names = TRUE)
tetramer_name=gsub("_50k_bulk_reference_mapped_obj.RDS","",basename(ref_map))
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/5_old_7_young/TFs/analysis/reference_mapping/UMAP/500_cells/"

for (i in 1:length(tetramer_name)) {
  obj=readRDS(ref_map[i])
  if(ncol(obj) >500){
    obj.subsampled <- obj[, sample(colnames(obj), size =500, replace=F)]
    p1 <- DimPlot(obj, reduction = "ref.umap", group.by = "predicted.seurat_clusters", label = TRUE,
                  label.size = 3, repel = TRUE) + NoLegend() + ggtitle(paste(tetramer_name[i],"reference mapped",sep = " "))
    p2 <- DimPlot(obj.subsampled, reduction = "ref.umap", group.by = "predicted.seurat_clusters", label = TRUE,
                  label.size = 3, repel = TRUE) + NoLegend() + ggtitle(paste(tetramer_name[i],"reference downsampled",sep = " "))
    
    dir.create(paste(savedir,"saveRDS_obj",sep = ""),showWarnings = FALSE)
    saveRDS(obj.subsampled, paste(savedir,"saveRDS_obj/",tetramer_name[i],"_50k_bulk_reference_mapped_obj_downsampled.RDS",sep = ""))
    
    pdf(paste(savedir,tetramer_name[i],"_reference_mapping_bulk_50k_downsampled.pdf",sep = ""),width = 12, height = 6)
    print(p1 + p2)
    dev.off()
  }
}


### Combining ####
maindir <- savedir 
library(tidyverse)
library(ggrepel)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

antigen = c("Lyt1","Lyt2","Lat2","Lat3")
for (j in 1:length(antigen)) {
  rm(file_list)
  file_list <- list()
  lyt1_files = list.files(paste(maindir,"saveRDS_obj/",sep = ""), full.names = TRUE, pattern = antigen[j])
  young=c(40187,40189,40193,40191,40225,40222,40224)
  young_lyt1=lyt1_files[grep("40187|40189|40193|40191|40225|40222|40224",lyt1_files)]
  
  for (i in 1:length(young_lyt1)) {  
    obj = readRDS(young_lyt1[i])
    if (all(rownames(obj@reductions$ref.umap@cell.embeddings) == rownames(obj@meta.data))==TRUE) {
      UMAP=obj@reductions$ref.umap@cell.embeddings
      UMAP_df=as.data.frame(UMAP)
      UMAP_df$seurat_cluster=obj@meta.data$predicted.seurat_clusters
      UMAP_df$cell_name= rownames(UMAP_df)
      file_list[[i]]=UMAP_df
    }
  }
  
  filelist_df=as.data.frame(rbindlist(file_list))
  set.seed(1)
  filelist_df$seurat_cluster <- factor(filelist_df$seurat_cluster, levels = 0:18)
  
  filelist_df$colors <- gsub("^0$","khaki1",
                      gsub("^1$","yellow1",
                           gsub("^2$","goldenrod2",
                                gsub("^3$","forestgreen",
                                     gsub("^4$","tomato",
                                          gsub("^5$","steelblue1",
                                               gsub("^6$","firebrick3",
                                                    gsub("^7$","seagreen3",
                                                         gsub("^8$","gold4",
                                                              gsub("^9$","yellow3",
                                                                   gsub("^10$","steelblue1",
                                                                        gsub("^11$","red4",
                                                                             gsub("^12$","blue1",
                                                                                  gsub("^13$","firebrick1",
                                                                                       gsub("^14$","darkgoldenrod3",
                                                                                            gsub("^15$","gold",
                                                                                                 gsub("^16$","royalblue",
                                                                                                      gsub("^17$","limegreen",
                                                                                                           gsub("^18$","gold",
                                                                                                                filelist_df$seurat_cluster)))))))))))))))))))
  
  
  col <- as.character(filelist_df$colors)
  names(col) <- as.character(filelist_df$seurat_cluster)
  
  label.df_2 <- filelist_df %>% 
    group_by(seurat_cluster) %>% 
    summarize(refUMAP_1 = mean(refUMAP_1), refUMAP_2 = mean(refUMAP_2))
  
  set.seed(1)
  p1 <- ggplot(filelist_df,aes(x=refUMAP_1,y=refUMAP_2,label=seurat_cluster))+
    geom_point(aes(colour = seurat_cluster), size=0.5, show.legend = TRUE) +
    ggrepel::geom_label_repel(data = label.df_2, aes(label = seurat_cluster), size=5) + 
    ggtitle(paste(antigen[j]," Young Combined")) + NoLegend() +
    scale_color_manual(values=col) + 
    theme(panel.background = element_rect(fill = 'white',color='black'),
          plot.title = element_text(hjust = 0.5))
  
  p2 <- ggplot(filelist_df,aes(x=refUMAP_1,y=refUMAP_2,label=seurat_cluster))+
    geom_point(aes(colour = seurat_cluster), size=0.5, show.legend = TRUE) +
    #ggrepel::geom_label_repel(data = label.df_2, aes(label = seurat_cluster), size=5) + 
    ggtitle(paste(antigen[j]," Young Combined")) + NoLegend() +
    scale_color_manual(values=col) + 
    theme(panel.background = element_rect(fill = 'white',color='black'),
          plot.title = element_text(hjust = 0.5))
  
  dir.create(paste(maindir,"combine/",sep = ""),showWarnings = FALSE)
  dir.create(paste(maindir,"combine/UMAP",sep = ""),showWarnings = FALSE)
  pdf(paste(maindir,"combine/UMAP/",antigen[j],"_combine_young.pdf",sep = ""),width = 12, height = 6)
  print(p1+p2)
  dev.off()
  
  clus_cells=as.data.frame(table(filelist_df$seurat_cluster))
  colnames(clus_cells) = c("Cluster","Cell_number")
  clus_cells_order=clus_cells[order(-clus_cells$Cell_number),]
  
  dir.create(paste(maindir,"combine/Table/",sep = ""),showWarnings = FALSE)
  write.table(clus_cells_order,paste(maindir,"combine/Table/",antigen[j],"_Young_order.txt",sep = ""), 
              quote = FALSE, sep = "\t",col.names = T, row.names = F)
}


## Old

for (j in 1:length(antigen)) {
  rm(file_list)
  file_list <- list()
  lyt1_files = list.files(paste(maindir,"saveRDS_obj/",sep = ""), full.names = TRUE, pattern = antigen[j])
  old = c(40188,40197,40198,40199,40223)
  Old_lyt1=lyt1_files[grep("40188|40197|40198|40199|40223",lyt1_files)]
  
  for (i in 1:length(Old_lyt1)) {  
    main_obj = readRDS(Old_lyt1[i])
    obj <- main_obj[, sample(colnames(main_obj), size =400, replace=F)]
    if (all(rownames(obj@reductions$ref.umap@cell.embeddings) == rownames(obj@meta.data))==TRUE) {
      UMAP=obj@reductions$ref.umap@cell.embeddings
      UMAP_df=as.data.frame(UMAP)
      UMAP_df$seurat_cluster=obj@meta.data$predicted.seurat_clusters
      UMAP_df$cell_name= rownames(UMAP_df)
      file_list[[i]]=UMAP_df
    }
  }
  
  filelist_df=as.data.frame(rbindlist(file_list))
  set.seed(1)
  filelist_df$seurat_cluster <- factor(filelist_df$seurat_cluster, levels = 0:18)
  
  filelist_df$colors <- gsub("^0$","khaki1",
                             gsub("^1$","yellow1",
                                  gsub("^2$","goldenrod2",
                                       gsub("^3$","forestgreen",
                                            gsub("^4$","tomato",
                                                 gsub("^5$","steelblue1",
                                                      gsub("^6$","firebrick3",
                                                           gsub("^7$","seagreen3",
                                                                gsub("^8$","gold4",
                                                                     gsub("^9$","yellow3",
                                                                          gsub("^10$","steelblue1",
                                                                               gsub("^11$","red4",
                                                                                    gsub("^12$","blue1",
                                                                                         gsub("^13$","firebrick1",
                                                                                              gsub("^14$","darkgoldenrod3",
                                                                                                   gsub("^15$","gold",
                                                                                                        gsub("^16$","royalblue",
                                                                                                             gsub("^17$","limegreen",
                                                                                                                  gsub("^18$","gold",
                                                                                                                       filelist_df$seurat_cluster)))))))))))))))))))
  
  col <- as.character(filelist_df$colors)
  names(col) <- as.character(filelist_df$seurat_cluster)
  
  label.df_2 <- filelist_df %>% 
    group_by(seurat_cluster) %>% 
    summarize(refUMAP_1 = mean(refUMAP_1), refUMAP_2 = mean(refUMAP_2))
  
  set.seed(1)
  p1 <- ggplot(filelist_df,aes(x=refUMAP_1,y=refUMAP_2,label=seurat_cluster))+
    geom_point(aes(colour = seurat_cluster), size=0.5, show.legend = TRUE) +
    ggrepel::geom_label_repel(data = label.df_2, aes(label = seurat_cluster), size=5) + 
    ggtitle(paste(antigen[j]," old Combined")) + NoLegend() +
    scale_color_manual(values=col) + 
    theme(panel.background = element_rect(fill = 'white',color='black'),
          plot.title = element_text(hjust = 0.5))
  
  p2 <- ggplot(filelist_df,aes(x=refUMAP_1,y=refUMAP_2,label=seurat_cluster))+
    geom_point(aes(colour = seurat_cluster), size=0.5, show.legend = TRUE) +
    #ggrepel::geom_label_repel(data = label.df_2, aes(label = seurat_cluster), size=5) + 
    ggtitle(paste(antigen[j]," old Combined")) + NoLegend() +
    scale_color_manual(values=col) + 
    theme(panel.background = element_rect(fill = 'white',color='black'),
          plot.title = element_text(hjust = 0.5))
  
  dir.create(paste(maindir,"combine/UMAP",sep = ""),showWarnings = FALSE)
  pdf(paste(maindir,"combine/UMAP/",antigen[j],"_combined_old_400.pdf",sep = ""),width = 12, height = 6)
  print(p1+p2)
  dev.off()
  
  clus_cells=as.data.frame(table(filelist_df$seurat_cluster))
  colnames(clus_cells) = c("Cluster","Cell_number")
  clus_cells_order=clus_cells[order(-clus_cells$Cell_number),]
  
  write.table(clus_cells_order,paste(maindir,"combine/Table/",antigen[j],"_old_order_400.txt",sep = ""), 
              quote = FALSE, sep = "\t",col.names = T, row.names = F)
  
  
}






