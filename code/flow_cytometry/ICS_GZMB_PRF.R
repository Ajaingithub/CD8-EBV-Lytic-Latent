#### 50k Cells ######
### CD8 lytic clustering with Seurat and Phenograph 
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.

rawData_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/rawdata/bulk_50k_cells/"
bulk_csv_files = list.files(rawData_path,pattern = ".csv")
# dir.create("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/5_old_7_young/TFs/analysis",showWarnings = FALSE)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream_include_PRF_GZMB/"

column = c("FJComp.APC.eFluor.780.A....Eomes","FJComp.Alexa.Fluor.700.A....RUNX3","FJComp.BUV496.A....CCR7",
           "FJComp.BUV737.A....CD28","FJComp.BUV805.A....CD27","FJComp.BV421.A....TCF1","FJComp.BV605.A....CD45RA",
           "FJComp.BV650.A....CD45RO","FJComp.BV785.A....Tbet","FJComp.Pacific.Blue.A....Helios",
           "FJComp.alexa.fluor.594.A....Tox.Tox2","FJComp.BV510.A....Granzyme.B","FJComp.BV711.A....Perforin")

sample_id = gsub(".*.TF_|_2022_.*|_2023_.*","",bulk_csv_files)
Run_id = c("Exp05","Exp05","Exp06","Exp06","Exp06","Exp06","Exp06","Exp07",
           "Exp07","Exp11","Exp11","Exp12","Exp12","Exp12","Exp13","Exp13","Exp13")
age = c("Y","O","Y","Y","O","O","O","O","Y","O","O","Y","O","O","Y","Y","Y")
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

objects <- ls(pattern="_obj") %>% grep("sample",.,invert=TRUE, value=TRUE)
name <- gsub("_obj","",objects)

### merging the 7 young datasets
combine <- merge(get(objects[1]), y = c(get(objects[2]),get(objects[3]),get(objects[4]),get(objects[5]),get(objects[6]),
                                        get(objects[7]),get(objects[8]),get(objects[9]),get(objects[10]),get(objects[11]),
                                        get(objects[12]),get(objects[13]),get(objects[14]),get(objects[15]),get(objects[16]),
                                        get(objects[17])), 
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
                                      RNA_features = c("BV510.A....Granzyme.B","BV711.A....Perforin"),Assay = Assay,
                                      process = process, objname = objname, ncol = 3)


### Imputing outside the conda environment
RNA_integrated <- CD8_ADT_integrated
Assay <- "ADT"
dims = 10
DefaultAssay(RNA_integrated) <- Assay
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/miniconda3/envs/seurat_new/bin/python3.8")
library(Rmagic)
py_discover_config("magic")
RNA_integrated_impute <- magic(RNA_integrated, npca=dims) ## imputing the RNA data as for RNA PCs are 20

savedir <- paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream_include_PRF_GZMB/")
saveRDS(RNA_integrated_impute, paste(savedir,"saveRDS_obj/CD8_ADT_integrated_impute.RDS",sep = ""))

CD8_ADT_integrated <- RNA_integrated_impute

library(ggplot2)
library(ArchR)
rm(plot_list)
DefaultAssay(CD8_ADT_integrated) = "ADT"
plot_list = list()
for (i in 1:length(rownames(CD8_ADT_integrated))) {
  p = FeaturePlot(CD8_ADT_integrated,rownames(CD8_ADT_integrated)[i], reduction = "umap", raster = TRUE) + 
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
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
  p = FeaturePlot(CD8_ADT_integrated,rownames(CD8_ADT_integrated)[i], reduction = "umap", raster = TRUE) + 
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

CD8_ADT_integrated_NN <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream_include_PRF_GZMB/saveRDS_obj/CD8_ADT_integrated_NN_cluster_1.RDS")
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/scRNA_cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8,1)
res = 0.8
for (i in 1:length(res)) {
  process = paste("ADT_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  fact = list()
  fact[[1]] <- unique(CD8_ADT_integrated@meta.data$orig.ident)[order(unique(CD8_ADT_integrated@meta.data$orig.ident))]
  fact[[2]] <- c("O","Y")
  CD8_ADT_integrated_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD8_ADT_integrated_NN, dims = 10, res = res[i], saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_ADT","nFeature_ADT"), objname = "CD8_ADT_integrated_NN_cluster",
                                                       process = process, col_sel = c("orig.ident","age"), factorized = fact)
}

dir.create(paste(savedir,"saveRDS_obj/",sep = ""),showWarnings = FALSE)
saveRDS(CD8_ADT_integrated_NN_cluster,paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster_0.8.RDS",sep = ""))

## Making a heatmap
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/Heatmap.R")
AJ_heatmap(obj = CD8_ADT_integrated_NN_cluster, savedir = savedir, objname = "CD8_ADT_integrated_NN", assay = "integrated")

### Separating out all the clusters ###
reference=CD8_ADT_integrated_NN_cluster
for (i in 1:length(unique(reference@meta.data$seurat_clusters))) {
  j=i-1
  reference_subset=subset(reference, subset = seurat_clusters==j)
  p=DimPlot(reference_subset, raster = TRUE, pt.size = 5)
  
  pdf(paste(savedir,"UMAP/Bulk_cluster",j,".pdf",sep = ""), width = 6.5, height = 6)
  print(p)
  dev.off()
}

### Easier and better way then down
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream/"


p <- DimPlot(reference, reduction = "umap", label = TRUE, label.size = 5,
             cols = c( '0' = "khaki1", '1' = "yellow1", '2' = "goldenrod2", '3' = "forestgreen", '4' ='red4',
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
reference=readRDS(paste(savedir,"saveRDS_obj/reference.RDS",sep = ""))

tetramer_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/rawdata/EBV_specific/"
teteramer_files = list.files(tetramer_path, pattern = "csv",full.names = TRUE)
### Since it could be completed in one go and restarted running it again on remaining
teteramer_files <- teteramer_files[1:length(teteramer_files)]
tetramer_basename = basename(teteramer_files)
tetramer_name = gsub(".*TF_|_2022_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]|_2023_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]|.csv","",tetramer_basename)

column = c("FJComp.APC.eFluor.780.A....Eomes","FJComp.Alexa.Fluor.700.A....RUNX3","FJComp.BUV496.A....CCR7",
           "FJComp.BUV737.A....CD28","FJComp.BUV805.A....CD27","FJComp.BV421.A....TCF1","FJComp.BV605.A....CD45RA",
           "FJComp.BV650.A....CD45RO","FJComp.BV785.A....Tbet","FJComp.Pacific.Blue.A....Helios",
           "FJComp.alexa.fluor.594.A....Tox.Tox2","FJComp.BV510.A....Granzyme.B","FJComp.BV711.A....Perforin")

library(ggplot2)
for (i in 1:length(teteramer_files)) {
  query_csv <- read.csv(teteramer_files[i])
  query_csv = query_csv[,which(colnames(query_csv) %in% column)]
  colnames(query_csv) = gsub("FJComp.","",colnames(query_csv))
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
                repel = TRUE, raster = TRUE) + NoLegend() + ggtitle("Reference annotations")
  p2 <- DimPlot(query_obj, reduction = "ref.umap", group.by = "predicted.seurat_clusters", label = TRUE,
                label.size = 3, repel = TRUE) + NoLegend() + ggtitle(paste(tetramer_name[i],"reference mapped",sep = " "))
  
  pdf(paste(savedir,"UMAP/",tetramer_name[i],"_reference_mapping_bulk_50k.pdf",sep = ""),width = 12, height = 6)
  print(p1 + p2)
  dev.off()
}

### Downstream files #####
# here we are recoloring the UMAP with the desired color as in feature plot it will provide the different colors
# Also we are saving the number of cells mapped to the reference 
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream/"
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
  write.table(mapped_cluster_order,paste(savedir,"reference_mapping/table/",basename[i],"_cell_mapped.txt",sep = ""), 
              quote = FALSE,row.names = FALSE, col.names = TRUE, sep = "\t")
}

total_cluster =16
df <- data.frame(matrix(ncol=1,nrow = (total_cluster+1)))
colnames(df) <- c("Cluster")
df$Cluster <- c(0:total_cluster)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream/reference_mapping/table/"
table_files <- list.files(savedir,pattern = "cell_mapped.txt",
                          full.names = TRUE)
table_basename <- gsub("_cell_mapped.txt","",basename(table_files))
for (i in 1:length(table_basename)) {
  tab <- read.table(table_files[i], header = TRUE)
  colnames(tab) <- c("Cluster",table_basename[i])
  df=join(df,tab)
}

df[is.na(df)] <- 0

write.table(df,paste(savedir,"combined_table_matrix.txt",sep = ""), quote = F, row.names = F,col.names = T, sep = "\t")

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

#### upsampling ####
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/resampling.R")
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream/500_cells/"
dir.create(savedir,showWarnings = FALSE)

checkdir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream/original/"
dir.create(checkdir,showWarnings = FALSE)

column = c("FJComp.APC.eFluor.780.A....Eomes","FJComp.Alexa.Fluor.488.A....BLIMP","FJComp.Alexa.Fluor.700.A....RUNX3",
           "FJComp.BUV496.A....CCR7","FJComp.BUV661.A....BCL6","FJComp.BUV737.A....CD28","FJComp.BUV805.A....CD27",
           "FJComp.BV421.A....TCF1","FJComp.BV605.A....CD45RA","FJComp.BV650.A....CD45RO","FJComp.BV785.A....Tbet",
           "FJComp.Pacific.Blue.A....Helios","FJComp.alexa.fluor.594.A....Tox.Tox2")

tetramer_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/rawdata/EBV_specific/"

teteramer_files = list.files(tetramer_path, pattern = "csv",full.names = TRUE)
### Since it could be completed in one go and restarted running it again on remaining
teteramer_files <- teteramer_files[1:length(teteramer_files)]
tetramer_basename = basename(teteramer_files)
tetramer_name = gsub(".*TF_|_2022_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]|_2023_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]|.csv","",tetramer_basename)

for (i in 1:length(teteramer_files)) {
  tetramer = read.csv(teteramer_files[i])
  tetramer_sel = tetramer[,which(colnames(tetramer) %in% column)]
  write.csv(tetramer_sel, paste(checkdir,tetramer_name[i],"_original.csv",sep = ""),
            col.names = TRUE, row.names = FALSE, quote = FALSE)
  # rownames(tetramer_sel) = paste("cell",1:nrow(tetramer_sel),sep = "")
  if(nrow(tetramer_sel)<500){
    tetramer_sel_rem=resample_data(tetramer_sel, N = 500-nrow(tetramer_sel))
    tetramer_sel_500=rbind(tetramer_sel,tetramer_sel_rem)
    write.csv(tetramer_sel_500, paste(savedir,tetramer_name[i],"_500_cells_upsample.csv",sep = ""),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  } else if(nrow(tetramer_sel)>500){
    tetramer_sel_rem=resample_data(tetramer_sel, N = 500)
    write.csv(tetramer_sel_rem, paste(savedir,tetramer_name[i],"_500_cells_downsample.csv",sep = ""),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}

## Reference mapping for upsampled files #####
# As reference mapping takes a lot of time so we will do reference mapping of the upsampled files only, while downsampling we can use the same cell number as we have earlier
library(ggplot2)
library(Seurat)
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream/"
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream/reference_mapping/UMAP/500_cells_resample/"

dir.create(savedir,showWarnings = FALSE)
reference=readRDS(paste(maindir,"saveRDS_obj/reference.RDS",sep = ""))
tetramer_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream/500_cells/"
teteramer_files = list.files(tetramer_path, pattern = "csv",full.names = TRUE)

## Since we are only interested in the young and upsampled so we will extract those out
teteramer_files_upsample=teteramer_files[grep("upsample",basename(teteramer_files))]
# young
# teteramer_files_upsample_young=teteramer_files_upsample[grep("40187|40189|40193|40191|40225|40222|40224",teteramer_files_upsample)]
teteramer_files_upsample_young <- teteramer_files_upsample

tetramer_basename = basename(teteramer_files_upsample_young)
tetramer_name = gsub(".csv","",tetramer_basename)
column = c("FJComp.APC.eFluor.780.A....Eomes","FJComp.Alexa.Fluor.488.A....BLIMP","FJComp.Alexa.Fluor.700.A....RUNX3",
           "FJComp.BUV496.A....CCR7","FJComp.BUV661.A....BCL6","FJComp.BUV737.A....CD28","FJComp.BUV805.A....CD27",
           "FJComp.BV421.A....TCF1","FJComp.BV605.A....CD45RA","FJComp.BV650.A....CD45RO","FJComp.BV785.A....Tbet",
           "FJComp.Pacific.Blue.A....Helios","FJComp.alexa.fluor.594.A....Tox.Tox2")
column <- gsub("FJComp.","",column)

for (i in 1:length(tetramer_name)) {
  query_csv <- read.csv(teteramer_files_upsample_young[i])
  colnames(query_csv) = gsub("FJComp.","",colnames(query_csv))
  query_csv = query_csv[,which(colnames(query_csv) %in% column)]
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
  
  pdf(paste(savedir,tetramer_name[i],"_reference_mapping_bulk_50k.pdf",sep = ""),width = 12, height = 6)
  print(p1 + p2)
  dev.off()
}

### downsampling the RDS object #####
# As we have already performed the reference mapping and save the RDS object
ref_map=list.files(paste(maindir,"saveRDS_obj",sep = ""),pattern = "_50k_bulk_reference_mapped_obj.RDS",full.names = TRUE)
tetramer_name=gsub("_50k_bulk_reference_mapped_obj.RDS","",basename(ref_map))
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream/reference_mapping/UMAP/500_cells_resample/"

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

### Combining the Antigen specific UMAPs ####
# Since we are combining only the 500 cells from both the upsampled and downsampled files that we have wrote down code to this one
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream/reference_mapping/UMAP/500_cells_resample/"
saveRDS = list.files(paste(maindir,"saveRDS_obj/",sep = ""),pattern = ".RDS",full.names = TRUE)
# saveRDS_req=saveRDS[grep("40188",saveRDS,invert = TRUE)]
filename = gsub("_50k_bulk_reference_mapped_.*.RDS|_500_cells_upsample","",basename(saveRDS))

dir.create(paste(maindir,"combine/",sep = ""),showWarnings = FALSE)
savedir = paste(maindir,"combine/",sep = "")
dir.create(paste(savedir,"Table/",sep = ""),showWarnings = FALSE)
dir.create(paste(savedir,"graphs/",sep = ""),showWarnings = FALSE)

for (i in 1:length(saveRDS)) {
  A40225_lyt2 = readRDS(saveRDS[i])
  df = as.data.frame(table(A40225_lyt2@meta.data$predicted.seurat_clusters))
  colnames(df) = c("Cluster","Cellnumber")
  cellsum=sum(df$Cellnumber)
  df$percentage_cell = (df$Cellnumber/cellsum)*100
  write.table(df,paste(savedir,"Table/",filename[i],".txt",sep = ""),col.names = T, row.names = F, sep = "\t",quote = FALSE)
  
  p <- ggplot(df, aes(x=Cluster, y=percentage_cell)) +
    geom_bar(stat = "identity",color="black", fill="grey", position = "dodge") +
    theme_bw() + ggtitle(filename[i])
  
  pdf(paste(savedir,"graphs/",filename[i],"_bargraph.pdf",sep = ""), width = 8, height = 6)
  print(p)
  dev.off()
}

## Lyt1
antigen = "Lat3"
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream/reference_mapping/UMAP/500_cells_resample/"
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream/reference_mapping/UMAP/500_cells_resample/combine/"
library(tidyverse)
library(ggrepel)
rm(file_list)
file_list <- list()
lyt1_files = list.files(paste(maindir,"saveRDS_obj/",sep = ""), full.names = TRUE, pattern = antigen)

for (i in 1:length(lyt1_files)) {
  obj = readRDS(lyt1_files[i])
  if (all(rownames(obj@reductions$ref.umap@cell.embeddings) == rownames(obj@meta.data))==TRUE) {
    UMAP=obj@reductions$ref.umap@cell.embeddings
    UMAP_df=as.data.frame(UMAP)
    UMAP_df$seurat_cluster=obj@meta.data$predicted.seurat_clusters
    UMAP_df$cell_name= rownames(UMAP_df)
    file_list[[i]]=UMAP_df
  }
}
library(data.table)
filelist_df=as.data.frame(rbindlist(file_list))

library(dplyr)
library(ggplot2)
library(ggrepel)

set.seed(1)
filelist_df$seurat_cluster <- factor(filelist_df$seurat_cluster, levels = 0:16)

filelist_df$colors <- gsub("^0$","goldenrod2",
                           gsub("^1$","khaki1",
                                gsub("^2$","dodgerblue",
                                     gsub("^3$","tomato",
                                          gsub("^4$","seagreen3",
                                               gsub("^5$","blue1",
                                                    gsub("^6$","yellow1",
                                                         gsub("^7$","orangered",
                                                              gsub("^8$","steelblue1",
                                                                   gsub("^9$","red3",
                                                                        gsub("^10$","firebrick",
                                                                             gsub("^11$","steelblue1",
                                                                                  gsub("^12$","gold4",
                                                                                       gsub("^13$","darkblue",
                                                                                            gsub("^14$","forestgreen",
                                                                                                 gsub("^15$","tomato",
                                                                                                      gsub("^16$","yellow1",
                                                                                                           filelist_df$seurat_cluster)))))))))))))))))

col <- as.character(filelist_df$colors)
names(col) <- as.character(filelist_df$seurat_cluster)

label.df_2 <- filelist_df %>% 
  group_by(seurat_cluster) %>% 
  summarize(refUMAP_1 = mean(refUMAP_1), refUMAP_2 = mean(refUMAP_2))

young=c(40187,40189,40193,40515,40225,40025,40102,40104)
young_lyt1=lyt1_files[grep(paste(young,collapse = "|"),lyt1_files)]

filelist_df$Age <- "O"
filelist_df[grep(paste(young,collapse = "|"),filelist_df$cell_name),"Age"]  <- "Y"

set.seed(1)
p1 <- ggplot(filelist_df,aes(x=refUMAP_1,y=refUMAP_2,label=seurat_cluster)) + facet_wrap(~Age) +
  geom_point(aes(colour = seurat_cluster), size=0.5, show.legend = TRUE) +
  ggrepel::geom_label_repel(data = label.df_2, aes(label = seurat_cluster), size=5) + 
  ggtitle(paste(antigen,"Combined")) + NoLegend() +
  scale_color_manual(values=col) + 
  theme(panel.background = element_rect(fill = 'white',color='black'),
        plot.title = element_text(hjust = 0.5))

p2 <- ggplot(filelist_df,aes(x=refUMAP_1,y=refUMAP_2,label=seurat_cluster)) + facet_wrap(~Age) +
  geom_point(aes(colour = seurat_cluster), size=0.5, show.legend = TRUE) +
  #ggrepel::geom_label_repel(data = label.df_2, aes(label = seurat_cluster), size=5) + 
  ggtitle(paste(antigen,"Combined")) + NoLegend() +
  scale_color_manual(values=col) + 
  theme(panel.background = element_rect(fill = 'white',color='black'),
        plot.title = element_text(hjust = 0.5))

dir.create(paste(maindir,"combine/UMAP",sep = ""),showWarnings = FALSE)
pdf(paste(maindir,"combine/UMAP/",antigen,"_combined.pdf",sep = ""),width = 12, height = 12)
print(p1/p2)
dev.off()

### UMAP
library("Seurat")
library("ggplot2")
library("ArchR")
ICS = readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream/saveRDS_obj/reference.RDS")
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")

metadata <- ICS@meta.data
metadata$cellnames <- rownames(metadata)
metadata_equal <- metadata %>% group_by(orig.ident) %>% slice_sample(n=25000)

ICS_subset <- subset(ICS, cells = metadata_equal$cellnames)

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/29_Aug_2023/ICS/Downstream/"
pdf(paste(savedir,"UMAP/Age_split_lowcells.pdf",sep = ""), width = 11, height = 6)
DimPlot(ICS_subset, group.by = "age", split.by = "age", cols = c("darkturquoise","hotpink1"), raster = TRUE)
dev.off()
