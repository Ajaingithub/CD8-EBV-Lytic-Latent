#### 50k Cells ######
### CD8 lytic clustering with Seurat and Phenograph Transcription factor
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure

rawData_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/7_young/TFs/raw_data/"
bulk_csv_files = list.files(rawData_path,pattern = ".csv")
dir.create("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/7_young/TFs/analysis/",showWarnings = FALSE)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/7_young/TFs/analysis/"

young = c("40187","40189","40193","40191","40225","40222","40224")
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/R_functions_Abhinav/partial_matching_Abhinav.R")
bulk_csv_files_young=bulk_csv_files[partial_match_Abhinav(young,bulk_csv_files)]

young_id = gsub(".*.cytoTF_|.*._TF_|_2022_.*","",bulk_csv_files_young)
Run_id = c("Exp01","Exp02","Exp02","Exp02","Exp03","Exp03","Exp03")
name=paste("young",young_id,Run_id,sep = "_")

column = c("APC.eFluor.780.A....Eomes","Alexa.Fluor.488.A....BLIMP","Alexa.Fluor.700.A....RUNX3","BUV496.A....CCR7","BUV661.A....BCL6",
           "BUV737.A....CD28","BUV805.A....CD27","BV421.A....TCF1","BV605.A....CD45RA","BV650.A....CD45RO","BV785.A....Tbet",
           "Pacific.Blue.A....Helios","alexa.fluor.594.A....Tox.Tox2")

for (i in 1:length(bulk_csv_files_young)) {
  bulk_sample = read.csv(paste(rawData_path,bulk_csv_files_young[i],sep = ""))
  bulk_sample = bulk_sample[,which(colnames(bulk_sample) %in% column)]
  rownames(bulk_sample) = paste("cell",1:nrow(bulk_sample),sep = "")
  bulk_sample_obj = CreateSeuratObject(t(bulk_sample), project = name[i], assay = "ADT")
  obj_name= paste(name[i],"obj",sep = "_")
  assign(obj_name,bulk_sample_obj)
}

### merging the 9 dataset 
combine <- merge(young_40187_Exp01_obj, y = c(young_40189_Exp02_obj,young_40191_Exp02_obj,young_40193_Exp02_obj,
                                              young_40222_Exp03_obj,young_40224_Exp03_obj,young_40225_Exp03_obj), 
                 add.cell.ids = c("young_40187_Exp01","young_40189_Exp02","young_40191_Exp02","young_40193_Exp02",
                                  "young_40222_Exp03","young_40224_Exp03","young_40225_Exp03"), project = "Combined")

combine@meta.data$age <- combine@meta.data$orig.ident
combine@meta.data$age = gsub("_.*","",combine@meta.data$age)

combine@meta.data$Run = combine@meta.data$orig.ident
combine@meta.data$Run = gsub(".*_","",combine@meta.data$Run)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")
objname = "CD8_TF"
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

process = paste("ADT_UMAP_QC_0.6",sep = "_")
Assay = "integrated"
CD8_ADT_integrated_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD8_ADT_integrated_NN, dims = 10, res = 0.6, saveDir=savedir, Assay = Assay,
                                                     QC_features = c("nCount_ADT","nFeature_ADT"), objname = "CD8_ADT_integrated_NN_cluster",
                                                     process = process)

dir.create(paste(savedir,"saveRDS_obj/",sep = ""),showWarnings = FALSE)
saveRDS(CD8_ADT_integrated_NN_cluster,paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster_0.6.RDS",sep = ""))

reference=readRDS(paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster_0.6.RDS",sep = ""))
cluster_cells=as.data.frame(table(reference@meta.data$seurat_clusters))
colnames(cluster_cells) = c("Cluster","Number_of_cells")
write.table(cluster_cells,paste(savedir,"Table/Bulk_cluster_cell_number.txt",sep = ""), quote = FALSE,
            row.names = FALSE, col.names = TRUE)

n_cells <- FetchData(reference, 
                     vars = c("seurat_clusters", "orig.ident")) %>%
  dplyr::count(seurat_clusters, orig.ident) %>%
  tidyr::spread(seurat_clusters, n)

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/7_young/TFs/analysis/"
write.table(n_cells, paste(savedir, "Table/Bulk_reference_sample_and_cluster_cm.txt",sep = ""),
            quote = FALSE, row.names = FALSE, col.names = TRUE)


### Separating out all the clusters ###
for (i in 1:length(unique(reference@meta.data$seurat_clusters))) {
  j=i-1
  reference_subset=subset(reference, subset = seurat_clusters==j)
  p=DimPlot(reference_subset)
  
  pdf(paste(savedir,"UMAP/Bulk_cluster",j,".pdf",sep = ""), width = 6.5, height = 6)
  print(p)
  dev.off()
}

### Easier and better way to make bulk UMAP ###

p <- DimPlot(reference, reduction = "umap", label = TRUE, label.size = 5,
             cols = c( '0' = "gold4", '1' = "khaki1", '2' = "seagreen3", '3' = "yellow3", '4' = "red4", 
                       "5"="blue1", "6"="tomato", "7"="goldenrod2", "8" ="royalblue","9"="gold", 
                       "10"="yellow1", "11"="firebrick1", "12"="steelblue1", "13"="forestgreen", "14"="darkgoldenrod3")) + NoLegend()

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/7_young/surface_panel/analysis/"
pdf(paste(savedir,"UMAP/bulk_mapping_3.pdf",sep = ""), width = 6, height = 6)
p
dev.off()

### Making a bulk UMAP specific colours given by Ines
UMAP=as.data.frame(reference@reductions$umap@cell.embeddings)
UMAP$cluster=reference@meta.data$seurat_clusters
UMAP$colors <- gsub("^0$","gold4",
                    gsub("^1$","khaki1",
                         gsub("^2$","seagreen3",
                              gsub("^3$","yellow3",
                                   gsub("^4$","red4",
                                        gsub("^5$","blue1",
                                             gsub("^6$","tomato",
                                                  gsub("^7$","goldenrod2",
                                                       gsub("^8$","royalblue",
                                                            gsub("^9$","gold",
                                                                 gsub("^10$","yellow1",
                                                                      gsub("^11$","firebrick1",
                                                                           gsub("^12$","steelblue1",
                                                                                gsub("^13$","forestgreen",
                                                                                     gsub("^14$","darkgoldenrod3",
                                                                                          UMAP$cluster)))))))))))))))

col <- as.character(UMAP$colors)
names(col) <- as.character(UMAP$cluster)

library(dplyr)
label.df_2 <- UMAP %>% 
  group_by(cluster) %>% 
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

p1 <- ggplot(UMAP,aes(x=UMAP_1,y=UMAP_2,label=cluster)) + 
  geom_point(aes(colour = cluster), size=0.01, show.legend = TRUE) +
  scale_color_manual(values=col) + NoLegend() + 
  ggrepel::geom_label_repel(data = label.df_2, aes(label = cluster), size=3) +
  ggtitle(paste("Reference Annotation")) +
  theme(panel.background = element_rect(fill = 'white',color='black'),
        plot.title = element_text(hjust = 0.5))

pdf(paste(savedir,"UMAP/Bulk_mapping_2.pdf",sep = ""),width = 6, height = 6)
print(p1)
dev.off()


# CD8_ADT_integrated_NN_cluster=readRDS(paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster_0.6.RDS",sep = ""))
CD8_counts = CD8_ADT_integrated_NN_cluster@assays$ADT@counts
all(colnames(CD8_counts) == rownames(CD8_ADT_integrated_NN_cluster@meta.data)) ## if TRUE move ahead

## Making an empty dataframe for heatmap
CD8_df = data.frame(matrix(ncol=length(unique(CD8_ADT_integrated_NN_cluster@meta.data$seurat_clusters)),nrow = nrow(CD8_counts)))
rownames(CD8_df) = rownames(CD8_counts) 
colnames(CD8_df) = paste("C",0:(length(unique(CD8_ADT_integrated_NN_cluster@meta.data$seurat_clusters))-1),sep = "")

for (i in 1:nrow(CD8_df)) {
  for (j in 1:ncol(CD8_df)) {
    CD8_df[i,j] = median(CD8_counts[i,grep(paste("^",j-1,"$",sep = ""),CD8_ADT_integrated_NN_cluster@meta.data$seurat_clusters)]) 
  }
}

write.table(CD8_df,paste(savedir,"heatmap/cluster_protein_median.txt",sep = ""),quote = FALSE,col.names = T, row.names = T, sep = "\t")

CD8_df_normalize = data.frame(matrix(ncol=length(unique(CD8_ADT_integrated_NN_cluster@meta.data$seurat_clusters)),nrow = nrow(CD8_counts)))
rownames(CD8_df_normalize) = rownames(CD8_counts) 
colnames(CD8_df_normalize) = paste("C",0:(length(unique(CD8_ADT_integrated_NN_cluster@meta.data$seurat_clusters))-1),sep = "")

i=1
for (i in 1:nrow(CD8_df)) {
  CD8_df_normalize[i,] = CD8_df[i,]/rowSums(CD8_df)[i]
}

write.table(CD8_df_normalize,paste(savedir,"heatmap/cluster_protein_median_normalize.txt",sep = ""),quote = FALSE,col.names = T, row.names = T, sep = "\t")

CD8_df_scaled = scale(CD8_df_normalize)
write.table(CD8_df_scaled,paste(savedir,"heatmap/cluster_protein_median_normalize_scaled.txt",sep = ""),quote = FALSE,col.names = T, row.names = T, sep = "\t")

library(ComplexHeatmap)
h <- Heatmap(CD8_df_scaled, 
             name = "z-score",
             column_title = "Seurat_Clusters",
             row_title = "Surface_protein",
             row_names_gp = gpar(fontsize = 12),
             column_names_gp = gpar(fontsize = 12),
             cluster_columns = TRUE,
             cluster_rows = TRUE
)

dir.create(paste(savedir,"heatmap",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"heatmap/cluster_and_protein_normalize_2.pdf",sep = ""), width = 10, height = 8)
h
dev.off()

#### 50k Query and Reference Mapping #####
reference = CD8_ADT_integrated_NN_cluster
reference = RunUMAP(reference, dims = 1:10, reduction = "pca", return.model = TRUE)
saveRDS(reference,paste(savedir,"saveRDS_obj/reference.RDS",sep = ""))
reference=readRDS(paste(savedir,"saveRDS_obj/reference.RDS",sep = ""))

tetramer_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_panel_tetramer/"
teteramer_files = list.files(tetramer_path, pattern = "csv",full.names = TRUE)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/R_functions_Abhinav/Bash_grep.R")
teteramer_files_young = teteramer_files[bash_grep_Abhinav(young,as.data.frame(teteramer_files))]
tetramer_basename = basename(teteramer_files_young)
tetramer_name = gsub(".*TF_|_2022_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]|.csv","",tetramer_basename)

for (i in 1:length(teteramer_files_young)) {
  query_csv <- read.csv(teteramer_files_young[i])
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
  UMAP$colors <- gsub("^0$","gold4",
                      gsub("^1$","khaki1",
                           gsub("^2$","seagreen3",
                                gsub("^3$","yellow3",
                                     gsub("^4$","red4",
                                          gsub("^5$","blue1",
                                               gsub("^6$","tomato",
                                                    gsub("^7$","goldenrod2",
                                                         gsub("^8$","royalblue",
                                                              gsub("^9$","gold",
                                                                   gsub("^10$","yellow1",
                                                                        gsub("^11$","firebrick1",
                                                                             gsub("^12$","steelblue1",
                                                                                  gsub("^13$","forestgreen",
                                                                                       gsub("^14$","darkgoldenrod3",
                                                                                            UMAP$predicted_seurat_cluster)))))))))))))))
  
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


#### Extra 3 TFs ######
### Other 3 TF that has not been considered for clustering and UMAP but we will use it for the feature plots
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure

rawData_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/7_young/TFs/raw_data/"
bulk_csv_files = list.files(rawData_path,pattern = ".csv")
dir.create("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/7_young/TFs/analysis/",showWarnings = FALSE)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/7_young/TFs/analysis/"

young = c("40187","40189","40193","40191","40225","40222","40224")
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/R_functions_Abhinav/partial_matching_Abhinav.R")
bulk_csv_files_young=bulk_csv_files[partial_match_Abhinav(young,bulk_csv_files)]

young_id = gsub(".*.cytoTF_|.*._TF_|_2022_.*","",bulk_csv_files_young)
Run_id = c("Exp01","Exp02","Exp02","Exp02","Exp03","Exp03","Exp03")
name=paste("young",young_id,Run_id,sep = "_")

column = c("APC.eFluor.780.A....Eomes","Alexa.Fluor.488.A....BLIMP","Alexa.Fluor.700.A....RUNX3","BUV496.A....CCR7","BUV661.A....BCL6",
           "BUV737.A....CD28","BUV805.A....CD27","BV421.A....TCF1","BV605.A....CD45RA","BV650.A....CD45RO","BV785.A....Tbet",
           "Pacific.Blue.A....Helios","alexa.fluor.594.A....Tox.Tox2","BV711.A....Perforin","PE.Cy7.A....Granzyme.K",
           "BV510.A....Granzyme.B")

for (i in 1:length(bulk_csv_files)) {
  bulk_sample = read.csv(paste(rawData_path,bulk_csv_files[i],sep = ""))
  bulk_sample = bulk_sample[,which(colnames(bulk_sample) %in% column)]
  rownames(bulk_sample) = paste("cell",1:nrow(bulk_sample),sep = "")
  bulk_sample_obj = CreateSeuratObject(t(bulk_sample), project = name[i], assay = "ADT")
  obj_name= paste(name[i],"obj",sep = "_")
  assign(obj_name,bulk_sample_obj)
}

### merging the 9 dataset 
combine <- merge(young_40187_Exp01_obj, y = c(young_40189_Exp02_obj,young_40191_Exp02_obj,young_40193_Exp02_obj,
                                              young_40222_Exp03_obj,young_40224_Exp03_obj,young_40225_Exp03_obj), 
                 add.cell.ids = c("young_40187_Exp01","young_40189_Exp02","young_40191_Exp02","young_40193_Exp02",
                                  "young_40222_Exp03","young_40224_Exp03","young_40225_Exp03"), project = "Combined")

combine <- NormalizeData(combine, normalization.method = 'CLR', margin = 2)
combine <- FindVariableFeatures(combine, selection.method = "vst", nfeatures = 16)
combine <- ScaleData(combine, verbose = FALSE)
combine <- RunPCA(combine, npcs = 16, approx=FALSE)
reference=readRDS(paste(savedir,"saveRDS_obj/reference.RDS",sep = ""))

combine@reductions$umap=reference@reductions$umap
saveRDS(combine,paste(savedir,"saveRDS_obj/CD8_3_more_TFs.RDS",sep = ""))
                                
rm(plot_list)
library("ggplot2")
library("ArchR")
DefaultAssay(combine) = "ADT"
plot_list = list()
for (i in 1:length(rownames(combine))) {
  p = FeaturePlot(combine,rownames(combine)[i], reduction = "umap") + scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] = p
}


dir.create(paste(savedir,"featureplot",sep = ""), showWarnings = FALSE)
require(gridExtra)
pdf(paste(savedir,"featureplot/featureplot_ADT_SolarExtra_defaultAssay_ADT_combine_2.pdf",sep = ""), width = 20, height = 20)
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
             plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
             plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]], ncol = 4, nrow = 4)
dev.off()

# combine=readRDS(paste(savedir,"saveRDS_obj/CD8_3_more_TFs.RDS",sep = ""))
CD8_counts = combine@assays$ADT@counts
all(colnames(CD8_counts) == rownames(combine@meta.data)) ## if TRUE move ahead
all(rownames(combine@meta.data) == rownames(reference@meta.data))
combine@meta.data$seurat_clusters = reference@meta.data$seurat_clusters

## Making an empty dataframe for heatmap
CD8_df = data.frame(matrix(ncol=length(unique(combine@meta.data$seurat_clusters)),nrow = nrow(CD8_counts)))
rownames(CD8_df) = rownames(CD8_counts) 
colnames(CD8_df) = paste("C",0:(length(unique(combine@meta.data$seurat_clusters))-1),sep = "")

for (i in 1:nrow(CD8_df)) {
  for (j in 1:ncol(CD8_df)) {
    CD8_df[i,j] = median(CD8_counts[i,grep(paste("^",j-1,"$",sep = ""),combine@meta.data$seurat_clusters)]) 
  }
}

write.table(CD8_df,paste(savedir,"heatmap/cluster_protein_median_16_TFs.txt",sep = ""),
            quote = FALSE,col.names = T, row.names = T, sep = "\t")

CD8_df_normalize = data.frame(matrix(ncol=length(unique(combine@meta.data$seurat_clusters)),nrow = nrow(CD8_counts)))
rownames(CD8_df_normalize) = rownames(CD8_counts) 
colnames(CD8_df_normalize) = paste("C",0:(length(unique(combine@meta.data$seurat_clusters))-1),sep = "")


i=1
for (i in 1:nrow(CD8_df)) {
  CD8_df_normalize[i,] = CD8_df[i,]/rowSums(CD8_df)[i]
}
write.table(CD8_df_normalize,paste(savedir,"heatmap/cluster_protein_median_normalize_16_TFs.txt",sep = ""),
            quote = FALSE,col.names = T, row.names = T, sep = "\t")

CD8_df_scaled = scale(CD8_df_normalize)
write.table(CD8_df_scaled,paste(savedir,"heatmap/cluster_protein_median_normalize_scaled_16_TFs.txt",sep = ""),
            quote = FALSE,col.names = T, row.names = T, sep = "\t")

library(ComplexHeatmap)
h <- Heatmap(CD8_df_scaled, 
             name = "z-score",
             column_title = "Seurat_Clusters",
             row_title = "Surface_protein",
             row_names_gp = gpar(fontsize = 12),
             column_names_gp = gpar(fontsize = 12),
             cluster_columns = TRUE,
             cluster_rows = TRUE
)

dir.create(paste(savedir,"heatmap",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"heatmap/cluster_and_protein_normalize_with_3_TFs.pdf",sep = ""), width = 10, height = 8)
h
dev.off()

saveRDS(combine,paste(savedir,"saveRDS_obj/combined_with_16_features.RDS",sep = ""))
# 
# reference_UMAP = as.data.frame(reference@reductions$umap@cell.embeddings)
# norm_counts = t(combine@assays$ADT@data)
# all(rownames(reference_UMAP)==rownames(norm_counts))
# nrow_count_UMAP=merge(reference_UMAP,norm_counts, by="row.names", all=TRUE)
# 
# nrow_count_UMAP = read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_bulk/Analysis/normalized_counts_UMAP_16_features.txt", header = TRUE)
# p1 = ggplot(nrow_count_UMAP,aes(x=UMAP_1,y=UMAP_2,color=PE.Cy7.A....Granzyme.K)) + geom_point(size=0.0001) + scale_color_gradientn(colours = ArchRPalettes$solarExtra) + theme_bw()
# p2 = ggplot(nrow_count_UMAP,aes(x=UMAP_1,y=UMAP_2,color=BV510.A....Granzyme.B)) + geom_point(size=0.0001) + scale_color_gradientn(colours = ArchRPalettes$solarExtra) + theme_bw()
# p3 = ggplot(nrow_count_UMAP,aes(x=UMAP_1,y=UMAP_2,color=BV711.A....Perforin)) + geom_point(size=0.0001) + scale_color_gradientn(colours = ArchRPalettes$solarExtra) + theme_bw()
# 
# require(gridExtra)
# pdf(paste(savedir,"featureplot/extra_3_TF.pdf",sep = ""),width = 15, height = 6)
# grid.arrange(p1,p2,p3,ncol=3)
# dev.off()
#
# obj_name= paste(name[i],"obj",sep = "_")
# assign(obj_name,bulk_sample_obj)


#### Making a bubble plot
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/7_young/TFs/analysis/"
CD4_C_S <- read.table(paste(maindir,"Table/CD8_ADT_integrated_NN_cluster_integrated_ADT_UMAP_QC_0.6_cluster_and_samples.txt",
                            sep = ""), header = T)
colnames(CD4_C_S) <- c(paste("C",0:8,sep = ""))
# colnames(CD4_C_S) <- c(paste("C",0:(length(unique(CD4_modality_integrate_cluster@meta.data$seurat_clusters))-1),sep = ""))
CD4_C_S[is.na(CD4_C_S)] <- 0
CD4_C_S_sum <- as.vector(rowSums(CD4_C_S))

### Making an empty dataframe 
df <- data.frame(matrix(nrow = nrow(CD4_C_S), ncol = ncol(CD4_C_S)))
rownames(df) <- rownames(CD4_C_S)
colnames(df) <- colnames(CD4_C_S)

for (i in 1:nrow(CD4_C_S)) {
  df[i,] <- (CD4_C_S[i,]/CD4_C_S_sum[i])*100
}

df$Samples <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("Samples","Cluster","Cell_percentage")

df_melted$Sampleid = gsub("old_|young_","",df_melted$Samples)
df_melted$Age = gsub("_.*","",df_melted$Samples)
df_melted$Sampleid = factor(df_melted$Sampleid,levels = unique(df_melted$Sampleid))
df_melted$Age = factor(df_melted$Age,levels = unique(df_melted$Age))

p <- ggplot(data=df_melted, aes(x=Cluster, y=Sampleid, size=Cell_percentage, color = Age))

pdf(paste(maindir,"Table/CD8_bubble_plot_sample.pdf",sep = ""),width = 8, height = 7)
p + geom_point(alpha=.75) + 
  scale_size(range = c(1,15), breaks=seq(0,40,by=15)) +  
  theme_bw() + ggtitle("CD8_Cluster_Sample")
dev.off()

#### Colwise sum
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_bulk/Analysis/"
CD4_C_S <- read.table(paste(maindir,"Table/CD8_ADT_integrated_NN_cluster_integrated_ADT_UMAP_QC_0.4_cluster_and_samples.txt",
                            sep = ""), header = T)
colnames(CD4_C_S) <- c(paste("C",0:8,sep = ""))
# colnames(CD4_C_S) <- c(paste("C",0:(length(unique(CD4_modality_integrate_cluster@meta.data$seurat_clusters))-1),sep = ""))
CD4_C_S[is.na(CD4_C_S)] <- 0
CD4_C_S_sum <- as.vector(colSums(CD4_C_S))

### Making an empty dataframe 
df <- data.frame(matrix(nrow = nrow(CD4_C_S), ncol = ncol(CD4_C_S)))
rownames(df) <- rownames(CD4_C_S)
colnames(df) <- colnames(CD4_C_S)

for (i in 1:ncol(CD4_C_S)) {
  df[,i] <- (CD4_C_S[,i]/CD4_C_S_sum[i])*100
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

maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_bulk/Analysis/"
saveRDS = list.files(paste(maindir,"saveRDS_obj/",sep = ""),pattern = "_bulk_reference_mapped_obj.RDS",full.names = TRUE)
filename = gsub("_50k_bulk_reference_mapped_obj.RDS","",basename(saveRDS))

dir.create(paste(maindir,"Antigen_specific/",sep = ""),showWarnings = FALSE)
savedir = paste(maindir,"Antigen_specific/",sep = "")
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


### Combining 
## Lyt1
antigen = "Lat3"
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_bulk/Analysis/"
library(tidyverse)
library(ggrepel)
rm(file_list)
file_list <- list()
lyt1_files = list.files(paste(maindir,"saveRDS_obj/",sep = ""), full.names = TRUE, pattern = antigen)
young=c(40187,40189,40193,40225)
young_lyt1=lyt1_files[grep("40187|40189|40193|40225",lyt1_files)]

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

library(dplyr)
library(ggplot2)
library(ggrepel)

set.seed(1)
filelist_df$seurat_cluster <- factor(filelist_df$seurat_cluster)

label.df_2 <- filelist_df %>% 
  group_by(seurat_cluster) %>% 
  summarize(refUMAP_1 = mean(refUMAP_1), refUMAP_2 = mean(refUMAP_2))

set.seed(1)
p <- ggplot(filelist_df,aes(x=refUMAP_1,y=refUMAP_2,color=seurat_cluster))+geom_point(size=0.5) + theme_bw() +
  ggrepel::geom_label_repel(data = label.df_2, aes(label = seurat_cluster), size=5) + 
  ggtitle(paste(antigen," Young Combined")) + NoLegend() +
  scale_color_manual(values = c("aquamarine3","brown2","black","dodgerblue1","lightcoral",
                                "lightgreen","cyan4","chartreuse3",
                                "bisque3"))
# "aquamarine3","brown2","black","dodgerblue1","lightcoral",
# "lightgreen","burlywood","firebrick3","cyan4","chartreuse3","bisque3","darksalmon",
# "darkseagreen4","mediumorchid3","lightslategray","red"
#                       # "C17" = "brown1",
# "C18" = "chocolate2",
# "C19" = "dodgerblue2",
# "C20" = "coral",
# "C21" = "green",
# "C22" = "burlywood3",
# "C23" = "firebrick2")
savedir=maindir
dir.create(paste(savedir,"Antigen_specific/combined",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"Antigen_specific/combined/",antigen,"_Young.pdf",sep = ""),width = 6, height = 5)
p
dev.off()

clus_cells=as.data.frame(table(filelist_df$seurat_cluster))
colnames(clus_cells) = c("Cluster","Cell_number")
clus_cells_order=clus_cells[order(-clus_cells$Cell_number),]

write.table(clus_cells_order,paste(savedir,"Antigen_specific/combined/",antigen,"_Young_order.txt",sep = ""), 
            quote = FALSE, sep = "\t",col.names = T, row.names = F)

## Old
old = c(40188,40197,40198,40199,40223)
Old_lyt1=lyt1_files[grep("40188|40197|40198|40199|40223",lyt1_files)]

rm(file_list)
file_list = list()
for (i in 1:length(Old_lyt1)) {
  obj = readRDS(Old_lyt1[i])
  if (all(rownames(obj@reductions$ref.umap@cell.embeddings) == rownames(obj@meta.data))==TRUE) {
    UMAP=obj@reductions$ref.umap@cell.embeddings
    UMAP_df=as.data.frame(UMAP)
    UMAP_df$seurat_cluster=obj@meta.data$predicted.seurat_clusters
    UMAP_df$cell_name= rownames(UMAP_df)
    file_list[[i]]=UMAP_df
  }
}
filelist_df=as.data.frame(rbindlist(file_list))

library(dplyr)
library(ggplot2)
library(ggrepel)

set.seed(1)
filelist_df$seurat_cluster <- factor(filelist_df$seurat_cluster)

label.df_2 <- filelist_df %>% 
  group_by(seurat_cluster) %>% 
  summarize(refUMAP_1 = mean(refUMAP_1), refUMAP_2 = mean(refUMAP_2))

set.seed(1)
p <- ggplot(filelist_df,aes(x=refUMAP_1,y=refUMAP_2,color=seurat_cluster))+geom_point(size=0.5) + theme_bw() +
  ggrepel::geom_label_repel(data = label.df_2, aes(label = seurat_cluster), size=5) + 
  ggtitle(paste(antigen,"Old Combined")) + NoLegend() +
  scale_color_manual(values = c("aquamarine3","brown2","black","dodgerblue1","lightcoral",
                                "lightgreen","firebrick3","cyan4","chartreuse3"))
# "aquamarine3","brown2","black","dodgerblue1","lightcoral",
# "lightgreen","burlywood","firebrick3","cyan4","chartreuse3","bisque3","darksalmon",
# "darkseagreen4","mediumorchid3","lightslategray","red"
#                       # "C17" = "brown1",
# "C18" = "chocolate2",
# "C19" = "dodgerblue2",
# "C20" = "coral",
# "C21" = "green",
# "C22" = "burlywood3",
# "C23" = "firebrick2")

dir.create(paste(savedir,"Antigen_specific/combined",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"Antigen_specific/combined/",antigen,"_Old.pdf",sep = ""),width = 6, height = 5)
p
dev.off()

clus_cells=as.data.frame(table(filelist_df$seurat_cluster))
colnames(clus_cells) = c("Cluster","Cell_number")
clus_cells_order=clus_cells[order(-clus_cells$Cell_number),]

write.table(clus_cells_order,paste(savedir,"Antigen_specific/combined/",antigen,"_Old_order.txt",sep = ""), 
            quote = FALSE, sep = "\t",col.names = T, row.names = F)


#### upsampling ####
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/resampling.R")
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_panel_tetramer/500_cells/"
dir.create(savedir,showWarnings = FALSE)

checkdir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_panel_tetramer/original/"
dir.create(checkdir,showWarnings = FALSE)

tetramer_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_panel_tetramer/"
files=list.files(tetramer_path,pattern = ".csv",full.names = TRUE)
# file_less=files[grep("40187_.*._Lat3|40189_.*._Specific_Lat3|40225_.*_Lat2|40225_.*._Lat3|40225_.*._Lyt2",files,ignore.case = TRUE)]
file_name=basename(files)
file_name_2=gsub("export_EBV_.*.TF_|2022_.*.[0-9][0-9]_|.csv","",file_name)

column = c("APC.eFluor.780.A....Eomes","Alexa.Fluor.488.A....BLIMP","Alexa.Fluor.700.A....RUNX3","BUV496.A....CCR7","BUV661.A....BCL6",
           "BUV737.A....CD28","BUV805.A....CD27","BV421.A....TCF1","BV605.A....CD45RA","BV650.A....CD45RO","BV785.A....Tbet",
           "Pacific.Blue.A....Helios","alexa.fluor.594.A....Tox.Tox2","BV711.A....Perforin","PE.Cy7.A....Granzyme.K",
           "BV510.A....Granzyme.B")

for (i in 1:length(file_name_2)) {
  tetramer = read.csv(files[i])
  tetramer_sel = tetramer[,which(colnames(tetramer) %in% column)]
  write.csv(tetramer_sel, paste(checkdir,file_name_2[i],"original.csv",sep = ""),
            col.names = TRUE, row.names = FALSE, quote = FALSE)
  # rownames(tetramer_sel) = paste("cell",1:nrow(tetramer_sel),sep = "")
  if(nrow(tetramer_sel)<500){
    tetramer_sel_rem=resample_data(tetramer_sel, N = 500-nrow(tetramer_sel))
    tetramer_sel_500=rbind(tetramer_sel,tetramer_sel_rem)
    write.csv(tetramer_sel_500, paste(savedir,file_name_2[i],"_500_cells_upsample.csv",sep = ""),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  } else if(nrow(tetramer_sel)>500){
    tetramer_sel_rem=resample_data(tetramer_sel, N = 500)
    write.csv(tetramer_sel_rem, paste(savedir,file_name_2[i],"_500_cells_downsample.csv",sep = ""),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}

## Reference mapping for upsampled files #####
# As reference mapping takes a lot of time so we will do reference mapping of the upsampled files only, while downsampling we can use the same cell number as we have earlier
library(ggplot2)
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/7_young/TFs/analysis/"
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/7_young/TFs/analysis/reference_mapping/UMAP/500_cells_resample/"

dir.create(savedir,showWarnings = FALSE)
reference=readRDS(paste(maindir,"saveRDS_obj/reference.RDS",sep = ""))
tetramer_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_panel_tetramer/500_cells/"
teteramer_files = list.files(tetramer_path, pattern = "csv",full.names = TRUE)

## Since we are only interested in the young and upsampled so we will extract those out
teteramer_files_upsample=teteramer_files[grep("upsample",basename(teteramer_files))]
# young
teteramer_files_upsample_young=teteramer_files_upsample[grep("40187|40189|40193|40191|40225|40222|40224",teteramer_files_upsample)]

tetramer_basename = basename(teteramer_files_upsample_young)
tetramer_name = gsub(".csv","",tetramer_basename)
column = c("APC.eFluor.780.A....Eomes","Alexa.Fluor.488.A....BLIMP","Alexa.Fluor.700.A....RUNX3","BUV496.A....CCR7","BUV661.A....BCL6",
           "BUV737.A....CD28","BUV805.A....CD27","BV421.A....TCF1","BV605.A....CD45RA","BV650.A....CD45RO","BV785.A....Tbet",
           "Pacific.Blue.A....Helios","alexa.fluor.594.A....Tox.Tox2")

for (i in 1:length(tetramer_name)) {
  query_csv <- read.csv(teteramer_files_upsample_young[i])
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

### downsampling the RDS object #####
# As we have already performed the reference mapping and save the RDS object
ref_map=list.files(paste(maindir,"saveRDS_obj",sep = ""),pattern = "_50k_bulk_reference_mapped_obj.RDS",full.names = TRUE)
tetramer_name=gsub("_50k_bulk_reference_mapped_obj.RDS","",basename(ref_map))
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/7_young/TFs/analysis/reference_mapping/UMAP/500_cells_resample/"

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
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/7_young/TFs/analysis/reference_mapping/UMAP/500_cells_resample/"
saveRDS = list.files(paste(maindir,"saveRDS_obj/",sep = ""),pattern = ".RDS",full.names = TRUE)
saveRDS_req=saveRDS[grep("40188",saveRDS,invert = TRUE)]
filename = gsub("_50k_bulk_reference_mapped_.*.RDS|_500_cells_upsample","",basename(saveRDS_req))

dir.create(paste(maindir,"combine/",sep = ""),showWarnings = FALSE)
savedir = paste(maindir,"combine/",sep = "")
dir.create(paste(savedir,"Table/",sep = ""),showWarnings = FALSE)
dir.create(paste(savedir,"graphs/",sep = ""),showWarnings = FALSE)

for (i in 1:length(saveRDS_req)) {
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
antigen = "Lat2"
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/7_young/TFs/analysis/reference_mapping/UMAP/500_cells_resample/"
library(tidyverse)
library(ggrepel)
rm(file_list)
file_list <- list()
lyt1_files = list.files(paste(maindir,"saveRDS_obj/",sep = ""), full.names = TRUE, pattern = antigen)
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
library(data.table)
filelist_df=as.data.frame(rbindlist(file_list))

library(dplyr)
library(ggplot2)
library(ggrepel)

set.seed(1)
filelist_df$seurat_cluster <- factor(filelist_df$seurat_cluster, levels = 0:14)

filelist_df$colors <- gsub("^0$","gold4",
                    gsub("^1$","khaki1",
                         gsub("^2$","seagreen3",
                              gsub("^3$","yellow3",
                                   gsub("^4$","red4",
                                        gsub("^5$","blue1",
                                             gsub("^6$","tomato",
                                                  gsub("^7$","goldenrod2",
                                                       gsub("^8$","royalblue",
                                                            gsub("^9$","gold",
                                                                 gsub("^10$","yellow1",
                                                                      gsub("^11$","firebrick1",
                                                                           gsub("^12$","steelblue1",
                                                                                gsub("^13$","forestgreen",
                                                                                     gsub("^14$","darkgoldenrod3",
                                                                                          filelist_df$seurat_cluster)))))))))))))))


col <- as.character(filelist_df$colors)
names(col) <- as.character(filelist_df$seurat_cluster)

label.df_2 <- filelist_df %>% 
  group_by(seurat_cluster) %>% 
  summarize(refUMAP_1 = mean(refUMAP_1), refUMAP_2 = mean(refUMAP_2))

set.seed(1)
p1 <- ggplot(filelist_df,aes(x=refUMAP_1,y=refUMAP_2,label=seurat_cluster))+
  geom_point(aes(colour = seurat_cluster), size=0.5, show.legend = TRUE) +
  ggrepel::geom_label_repel(data = label.df_2, aes(label = seurat_cluster), size=5) + 
  ggtitle(paste(antigen," Young Combined")) + NoLegend() +
  scale_color_manual(values=col) + 
  theme(panel.background = element_rect(fill = 'white',color='black'),
        plot.title = element_text(hjust = 0.5))

p2 <- ggplot(filelist_df,aes(x=refUMAP_1,y=refUMAP_2,label=seurat_cluster))+
  geom_point(aes(colour = seurat_cluster), size=0.5, show.legend = TRUE) +
  #ggrepel::geom_label_repel(data = label.df_2, aes(label = seurat_cluster), size=5) + 
  ggtitle(paste(antigen," Young Combined")) + NoLegend() +
  scale_color_manual(values=col) + 
  theme(panel.background = element_rect(fill = 'white',color='black'),
        plot.title = element_text(hjust = 0.5))

dir.create(paste(maindir,"combine/UMAP",sep = ""),showWarnings = FALSE)
pdf(paste(maindir,"combine/UMAP/",antigen,"_combined.pdf",sep = ""),width = 12, height = 6)
print(p1+p2)
dev.off()

#### Synthetic Data Generation #####
## Random oversampling is overfitting the data and does not look good. We will synthetically generate the data now.
library(synthpop)
synsavedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_panel_tetramer/syn_sample/"
dir.create(synsavedir,showWarnings = FALSE)

tetramer_dir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_panel_tetramer/"
files=list.files(tetramer_dir,pattern = ".csv",full.names = TRUE)
file_less=files[grep("40187_.*._Lat3|40189_.*._Specific_Lat3|40225_.*_Lat2|40225_.*._Lat3|40225_.*._Lyt2",files,ignore.case = TRUE)]
file_name=basename(file_less)
file_name_2=gsub("export_EBV_.*.TF_|2022_.*.[0-9][0-9]_|.csv","",file_name)

column = c("APC.eFluor.780.A....Eomes","Alexa.Fluor.488.A....BLIMP","Alexa.Fluor.700.A....RUNX3","BUV496.A....CCR7","BUV661.A....BCL6",
           "BUV737.A....CD28","BUV805.A....CD27","BV421.A....TCF1","BV605.A....CD45RA","BV650.A....CD45RO","BV785.A....Tbet",
           "Pacific.Blue.A....Helios","alexa.fluor.594.A....Tox.Tox2","BV711.A....Perforin","PE.Cy7.A....Granzyme.K",
           "BV510.A....Granzyme.B")

for (i in 1:length(file_less)) {
  tetramer = read.csv(file_less[i])
  tetramer_sel = tetramer[,which(colnames(tetramer) %in% column)]
  # colnames(tetramer_sel) = gsub("FJComp.","",colnames(tetramer_sel))
  tetramer_syn=syn(tetramer_sel, k=1000-nrow(tetramer_sel))
  tetramer_1000=rbind(tetramer_sel,tetramer_syn$syn)
  rownames(tetramer_1000) = paste("cell",1:nrow(tetramer_1000),sep = "")
  write.csv(tetramer_1000, paste(synsavedir,file_name_2[i],"_syn_sampled.csv",sep = ""), 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
}

## Reference mapping for synthetic sampled #####
library(ggplot2)
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_bulk/Analysis/"
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_bulk/Analysis/UMAP/syn_sample/"

dir.create(savedir,showWarnings = FALSE)
reference=readRDS(paste(maindir,"saveRDS_obj/reference.RDS",sep = ""))
tetramer_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Transcription_factor/TF_panel_tetramer/syn_sample/"
teteramer_files = list.files(tetramer_path, pattern = "csv",full.names = TRUE)

tetramer_basename = basename(teteramer_files)
tetramer_name = gsub(".csv","",tetramer_basename)
column = c("APC.eFluor.780.A....Eomes","Alexa.Fluor.488.A....BLIMP","Alexa.Fluor.700.A....RUNX3","BUV496.A....CCR7","BUV661.A....BCL6",
           "BUV737.A....CD28","BUV805.A....CD27","BV421.A....TCF1","BV605.A....CD45RA","BV650.A....CD45RO","BV785.A....Tbet",
           "Pacific.Blue.A....Helios","alexa.fluor.594.A....Tox.Tox2")

for (i in 1:length(tetramer_name)) {
  query_csv <- read.csv(teteramer_files[i])
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


