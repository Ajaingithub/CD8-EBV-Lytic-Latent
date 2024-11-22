#### 50k Cells ######
### CD8 lytic clustering with Seurat and Phenograph 
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.

rawData_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/clustering/raw_data/50000_downsampled_bulk/"
bulk_csv_files = list.files(rawData_path,pattern = ".csv")
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/clustering/raw_data/50000_downsampled_bulk/analysis/"

column = c("FJComp.APC.Fire.750.A....TIGIT","FJComp.Alexa.Fluor.488.A....CD69","FJComp.Alexa.Fluor.700.A....CD95","FJComp.BUV395.A....CD27","FJComp.BUV496.A....CCR7",
           "FJComp.BUV563.A....CD25","FJComp.BUV661.A....CX3CR1","FJComp.BUV737.A....CD28","FJComp.BUV805.A....CD85J","FJComp.BV421.A....CD137","FJComp.BV510.A....CD57",
           "FJComp.BV605.A....CD45RA","FJComp.BV650.A....CD45RO","FJComp.BV711.A....PD1","FJComp.BV785.A....KLRG1","FJComp.PE.Cy7.A....CD16","FJComp.PE.Dazzle594.A....TIM3",
           "FJComp.PE.Fire.700.A....CD127.IL7RA")

name = c("young_40187_Exp01","old_40188_Exp01","young_40189_Exp02","young_40193_Exp02","old_40197_Exp02","old_40198_Exp02","old_40199_Exp02","old_40223_Exp03","young_40225_Exp03")

for (i in 1:length(bulk_csv_files)) {
  bulk_sample = read.csv(paste(rawData_path,bulk_csv_files[i],sep = ""))
  bulk_sample = bulk_sample[,which(colnames(bulk_sample) %in% column)]
  colnames(bulk_sample) = gsub("FJComp.","",colnames(bulk_sample))
  rownames(bulk_sample) = paste("cell",1:nrow(bulk_sample),sep = "")
  bulk_sample_obj = CreateSeuratObject(t(bulk_sample), project = name[i], assay = "ADT")
  obj_name= paste(name[i],"obj",sep = "_")
  assign(obj_name,bulk_sample_obj)
}

### merging the 9 dataset 
combine <- merge(young_40187_Exp01_obj, y = c(old_40188_Exp01_obj,young_40189_Exp02_obj,young_40193_Exp02_obj,old_40197_Exp02_obj,
                                                                    old_40198_Exp02_obj,old_40199_Exp02_obj,old_40223_Exp03_obj,young_40225_Exp03_obj), 
                 add.cell.ids = c("young_40187_Exp01","old_40188_Exp01","young_40189_Exp02","young_40193_Exp02",
                                  "old_40197_Exp02","old_40198_Exp02","old_40199_Exp02","old_40223_Exp03",
                                  "young_40225_Exp03"), project = "Combined")

combine@meta.data$age <- combine@meta.data$orig.ident
combine@meta.data$age = gsub("_.*","",combine@meta.data$age)

combine@meta.data$Run = combine@meta.data$orig.ident
combine@meta.data$Run = gsub(".*_","",combine@meta.data$Run)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")
objname = "CD8"
process = "integrated"
Assay = "ADT"

CD4_ADT_integrated <- ADT_merging(combine, savedir, dims = 10, numfeatures=18,
                                  Assay=Assay, process=process, objname=objname, 
                                  sample_tree = NULL, split_by = "Run",
                                  reference=NULL)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "CD4"
process = "integration"
Assay = "ADT"
CD4_ADT_integrated <- RNA_integration(obj_path = CD4_ADT_integrated, saveDir = savedir, dims = 10,
                                      RNA_features = c("CD4-protein","CD8a-protein"),Assay = Assay,
                                      process = process, objname = objname, ncol = 3)


rm(plot_list)
plot_list = list()
for (i in 1:length(rownames(CD4_ADT_integrated))) {
  p = FeaturePlot(CD4_ADT_integrated,rownames(CD4_ADT_integrated)[i], reduction = "umap") + scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] = p
}


dir.create(paste(savedir,"featureplot",sep = ""), showWarnings = FALSE)
require(gridExtra)
pdf(paste(savedir,"featureplot/featureplot_ADT_SolarExtra.pdf",sep = ""), width = 20, height = 25)
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
             plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
             plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],plot_list[[17]],plot_list[[18]], ncol = 4, nrow = 5)
dev.off()


DefaultAssay(CD4_ADT_integrated) = "MAGIC_ADT"
rm(plot_list)
plot_list = list()
for (i in 1:length(rownames(CD4_ADT_integrated))) {
  p = FeaturePlot(CD4_ADT_integrated,rownames(CD4_ADT_integrated)[i], reduction = "umap")
  plot_list[[i]] = p
}

require(gridExtra)
pdf(paste(savedir,"featureplot/featureplot_MAGIC_ADT.pdf",sep = ""), width = 20, height = 25)
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
             plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
             plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],plot_list[[17]],plot_list[[18]], ncol = 4, nrow = 5)
dev.off()


### Nearest Neighbour and Clustering
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/NN_clustertree.R")
process = "ADT_NN"
Assay = "integrated"
CD4_ADT_integrated_NN <- nearest_neigbour(Obj_path = CD4_ADT_integrated, Dims = 10, saveDir = savedir,
                                          samplename = "CD4_ADT_integrated_NN", process = process, Assay = Assay)
# donot run the cluster tree takes a lot of time 
CD4_ADT_integrated_NN = GEX_obj_3

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8)
for (i in 1:length(res)) {
  process = paste("ADT_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  CD4_ADT_integrated_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD4_ADT_integrated_NN, dims = 10, res = res[i], saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_ADT","nFeature_ADT"), objname = "CD4_ADT_integrated_NN_cluster", 
                                                       process = process)
}

saveRDS(CD4_ADT_integrated_NN_cluster,paste(savedir,"saveRDS_obj/CD4_ADT_integrated_NN_cluster.RDS",sep = ""))

CD4_counts = CD4_ADT_integrated_NN_cluster@assays$ADT@counts

all(colnames(CD4_counts) == rownames(CD4_ADT_integrated_NN_cluster@meta.data)) ## if TRUE move ahead

## Making an empty dataframe for heatmap
CD4_df = data.frame(matrix(ncol=length(unique(CD4_ADT_integrated_NN_cluster@meta.data$seurat_clusters)),nrow = nrow(CD4_counts)))
rownames(CD4_df) = rownames(CD4_counts) 
colnames(CD4_df) = paste("C",0:16,sep = "")

for (i in 1:nrow(CD4_df)) {
  for (j in 1:17) {
    CD4_df[i,j] = median(CD4_counts[i,grep(j-1,CD4_ADT_integrated_NN_cluster@meta.data$seurat_clusters)]) 
  }
}

CD4_df_normalize = data.frame(matrix(ncol=length(unique(CD4_ADT_integrated_NN_cluster@meta.data$seurat_clusters)),nrow = nrow(CD4_counts)))
rownames(CD4_df_normalize) = rownames(CD4_counts) 
colnames(CD4_df_normalize) = paste("C",0:16,sep = "")

i=1
for (i in 1:nrow(CD4_df)) {
  CD4_df_normalize[i,] = CD4_df[i,]/rowSums(CD4_df)[i]
}

CD4_df_scaled = scale(CD4_df_normalize)
library(ComplexHeatmap)
h <- Heatmap(CD4_df_scaled, 
             name = "z-score",
             column_title = "Seurat_Clusters",
             row_title = "Surface_protein",
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             cluster_columns = TRUE,
             cluster_rows = TRUE
)

dir.create(paste(savedir,"heatmap",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"heatmap/cluster_and_protein_normalize.pdf",sep = ""), width = 10, height = 8)
h
dev.off()

#### 50k Query and Reference Mapping #####
reference = CD4_ADT_integrated_NN_cluster
reference = RunUMAP(reference, dims = 1:10, reduction = "pca", return.model = TRUE)
tetramer_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/clustering/raw_data/all_teteramer_pos/"
teteramer_files = list.files(tetramer_path, pattern = "csv",full.names = TRUE)
tetramer_basename = basename(teteramer_files)
tetramer_name = gsub(".*_surface_|_2022_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]|.csv","",tetramer_basename)

for (i in 2:length(teteramer_files)) {
  query_csv <- read.csv(teteramer_files[i])
  query_csv = query_csv[,which(colnames(query_csv) %in% column)]
  colnames(query_csv) = gsub("FJComp.","",colnames(query_csv))
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
  
  pdf(paste(savedir,"UMAP/",tetramer_name[i],"_reference_mapping_bulk_50k.pdf",sep = ""),width = 12, height = 6)
  print(p1 + p2)
  dev.off()
}

# obj_name= paste(name[i],"obj",sep = "_")
# assign(obj_name,bulk_sample_obj)


### 5 millions ####
### Making the reference with the 5 million cells also
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.


rawData_path = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/clustering/raw_data/5million_cells/"
savedir = paste(rawData_path,"analysis/",sep = "")
dir.create(savedir, showWarnings = FALSE)

bulk_csv_files = list.files(rawData_path,pattern = ".csv")
bulk_csv_files = bulk_csv_files[-1] # since we donot need the scaled matrix
column = c("FJComp.APC.Fire.750.A....TIGIT","FJComp.Alexa.Fluor.488.A....CD69","FJComp.Alexa.Fluor.700.A....CD95","FJComp.BUV395.A....CD27","FJComp.BUV496.A....CCR7",
           "FJComp.BUV563.A....CD25","FJComp.BUV661.A....CX3CR1","FJComp.BUV737.A....CD28","FJComp.BUV805.A....CD85J","FJComp.BV421.A....CD137","FJComp.BV510.A....CD57",
           "FJComp.BV605.A....CD45RA","FJComp.BV650.A....CD45RO","FJComp.BV711.A....PD1","FJComp.BV785.A....KLRG1","FJComp.PE.Cy7.A....CD16","FJComp.PE.Dazzle594.A....TIM3",
           "FJComp.PE.Fire.700.A....CD127.IL7RA")
column=gsub("FJ","",column)


name = c("young_40187_Exp01","old_40188_Exp01","young_40189_Exp02","young_40193_Exp02",
         "old_40197_Exp02","old_40198_Exp02","old_40199_Exp02","old_40223_Exp03","young_40225_Exp03")

for (i in 1:length(bulk_csv_files)) {
  bulk_sample = read.csv(paste(rawData_path,bulk_csv_files[i],sep = ""))
  bulk_sample = bulk_sample[,which(colnames(bulk_sample) %in% column)]
  colnames(bulk_sample) = gsub("Comp.","",colnames(bulk_sample))
  rownames(bulk_sample) = paste("cell",1:nrow(bulk_sample),sep = "")
  bulk_sample_obj = CreateSeuratObject(t(bulk_sample), project = name[i], assay = "ADT")
  obj_name= paste(name[i],"obj",sep = "_")
  assign(obj_name,bulk_sample_obj)
}

### merging the 9 dataset 
combine <- merge(young_40187_Exp01_obj, y = c(old_40188_Exp01_obj,young_40189_Exp02_obj,young_40193_Exp02_obj,old_40197_Exp02_obj,
                                              old_40198_Exp02_obj,old_40199_Exp02_obj,old_40223_Exp03_obj,young_40225_Exp03_obj), 
                 add.cell.ids = c("young_40187_Exp01","old_40188_Exp01","young_40189_Exp02","young_40193_Exp02",
                                  "old_40197_Exp02","old_40198_Exp02","old_40199_Exp02","old_40223_Exp03",
                                  "young_40225_Exp03"), project = "Combined")

combine@meta.data$age <- combine@meta.data$orig.ident
combine@meta.data$age = gsub("_.*","",combine@meta.data$age)

combine@meta.data$Run = combine@meta.data$orig.ident
combine@meta.data$Run = gsub(".*_","",combine@meta.data$Run)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")
objname = "CD8"
process = "integrated"
Assay = "ADT"

CD4_ADT_integrated <- ADT_merging(combine, savedir, dims = 10, numfeatures=18,
                                  Assay=Assay, process=process, objname=objname, 
                                  sample_tree = NULL, split_by = "Run",
                                  reference=NULL)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "CD4"
process = "integration"
Assay = "ADT"
CD4_ADT_integrated <- RNA_integration(obj_path = CD4_ADT_integrated, saveDir = savedir, dims = 10,
                                      RNA_features = c("CD4-protein","CD8a-protein"),Assay = Assay,
                                      process = process, objname = objname, ncol = 6)



### Nearest Neighbour and Clustering
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/NN_clustertree.R")
process = "ADT_NN"
Assay = "integrated"
CD4_ADT_integrated_NN <- nearest_neigbour(Obj_path = CD4_ADT_integrated, Dims = 10, saveDir = savedir,
                                          samplename = "CD4_ADT_integrated_NN", process = process, Assay = Assay)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8)
for (i in 1:length(res)) {
  process = paste("ADT_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  CD4_ADT_integrated_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD4_ADT_integrated_NN, dims = 10, res = res[i], saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_ADT","nFeature_ADT"), objname = "CD4_ADT_integrated_NN_cluster", 
                                                       process = process)
  
}

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/scRNA_cluster_UMAP_QC.R")
process = paste("ADT_UMAP_QC_0.8",sep = "_")
Assay = "integrated"
CD4_ADT_integrated_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD4_ADT_integrated_NN_cluster, dims = 10, res = 0.8, saveDir=savedir, Assay = Assay,
                                                     QC_features = c("nCount_ADT","nFeature_ADT"), objname = "CD4_ADT_integrated_NN_cluster",
                                                     process = process)


saveRDS(CD4_ADT_integrated_NN_cluster,paste(savedir,"saveRDS_obj/CD4_ADT_integrated_NN_cluster.RDS",sep = ""))


##### Performing Analysis of single sample #######
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.

column = c("FJComp.APC.Fire.750.A....TIGIT","FJComp.Alexa.Fluor.488.A....CD69","FJComp.Alexa.Fluor.700.A....CD95","FJComp.BUV395.A....CD27","FJComp.BUV496.A....CCR7",
           "FJComp.BUV563.A....CD25","FJComp.BUV661.A....CX3CR1","FJComp.BUV737.A....CD28","FJComp.BUV805.A....CD85J","FJComp.BV421.A....CD137","FJComp.BV510.A....CD57",
           "FJComp.BV605.A....CD45RA","FJComp.BV650.A....CD45RO","FJComp.BV711.A....PD1","FJComp.BV785.A....KLRG1","FJComp.PE.Cy7.A....CD16","FJComp.PE.Dazzle594.A....TIM3",
           "FJComp.PE.Fire.700.A....CD127.IL7RA")

bulk_sample = read.csv("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/clustering/raw_data/50000_downsampled_bulk/Matrix_BulkSample_EBV_phenotyping06_surface_40189_2022_04_28_12_23_33_down.csv")
bulk_sample = bulk_sample[,which(colnames(bulk_sample) %in% column)]
colnames(bulk_sample) = gsub("FJComp.","",colnames(bulk_sample))
rownames(bulk_sample) = paste("cell",1:nrow(bulk_sample),sep = "")
bulk_sample_obj = CreateSeuratObject(t(bulk_sample), project ="40189", assay = "ADT")

bulk_sample_obj <- NormalizeData(bulk_sample_obj, normalization.method = 'CLR', margin = 2)
bulk_sample_obj <- FindVariableFeatures(bulk_sample_obj, selection.method = "vst", nfeatures = 18)

bulk_sample_obj <- ScaleData(bulk_sample_obj, verbose = FALSE)
bulk_sample_obj <- RunPCA(bulk_sample_obj, npcs = 18, approx=FALSE)

bulk_sample_obj <- RunUMAP(bulk_sample_obj, dims = 1:10, reduction = "pca")
DimPlot(bulk_sample_obj)

rm(plot_list)
plot_list = list()
for (i in 1:length(rownames(bulk_sample_obj))) {
  p = FeaturePlot(bulk_sample_obj,rownames(bulk_sample_obj)[i], reduction = "umap")
  plot_list[[i]] = p
}

bulk_sample_obj <- FindNeighbors(bulk_sample_obj, dims = 1:10)
bulk_sample_obj <- FindClusters(bulk_sample_obj, resolution = 0.6, verbose = FALSE)
DimPlot(bulk_sample_obj)

grep("0",bulk_sample_obj@meta.data$seurat_clusters)


savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/clustering/raw_data/50000_downsampled_bulk/analysis/single_bulk/"
dir.create(savedir, showWarnings = F)
require(gridExtra)
pdf(paste(savedir,"featureplot_ADT.pdf",sep = ""), width = 20, height = 25)
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
             plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
             plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],plot_list[[17]],plot_list[[18]], ncol = 4, nrow = 5)
dev.off()

# matrix_scaled = read.csv("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/clustering/raw_data/50000_downsampled_bulk/Matrix_scaled_BulkSample_EBV_phenotyping06_surface_40189_2022_04_28_12_23_33_down.csv")
# rownames(matrix_scaled) = paste("cell",1:nrow(matrix_scaled),sep = "")
# matrix_scaled_sel = matrix_scaled[,which(colnames(matrix_scaled) %in% column)]
# cell_embedding = bulk_sample_obj@reductions$umap@cell.embeddings
# 
# cell_embedding_matrix_sclaed = merge(cell_embedding,matrix_scaled_sel, by="row.names")
 
pdf(paste(savedir,"UMAP/CD4_ADT_integrated_NN_cluster_Age_splitted.pdf",sep = ""),width = 10, height = 6)
DimPlot(CD4_ADT_integrated_NN_cluster, split.by = "age", group.by="age", ncol = 2)
dev.off()


#### Making a bubble plot
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/clustering/raw_data/50000_downsampled_bulk/analysis/"
CD4_C_S <- read.table(paste(maindir,"Table/CD4_ADT_integrated_NN_cluster_integrated_ADT_UMAP_QC_0.8_cluster_and_samples.txt",
                            sep = ""), header = T)
colnames(CD4_C_S) <- c(paste("C",0:16,sep = ""))
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

pdf(paste(maindir,"Table/CD4_bubble_plot_sample.pdf",sep = ""),width = 9, height = 8)
p + geom_point(alpha=.75) + 
  scale_size(range = c(1,15), breaks=seq(0,40,by=15)) +  
  theme_bw() + ggtitle("CD4_Cluster_Sample")
dev.off()

#### Colwise sum
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/clustering/raw_data/50000_downsampled_bulk/analysis/"
CD4_C_S <- read.table(paste(maindir,"Table/CD4_ADT_integrated_NN_cluster_integrated_ADT_UMAP_QC_0.8_cluster_and_samples.txt",
                            sep = ""), header = T)
colnames(CD4_C_S) <- c(paste("C",0:16,sep = ""))
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

pdf(paste(maindir,"Table/CD4_bubble_plot_ColSums.pdf",sep = ""),width = 9, height = 8)
p + geom_point(alpha=.75) + 
  scale_size(range = c(1,15), breaks=seq(0,100,by=15)) +  
  theme_bw() + ggtitle("CD4 modality Integrated ColSums")
dev.off()


maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/clustering/raw_data/50000_downsampled_bulk/analysis/"
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

## Lyt1
antigen = "Lat3"
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/clustering/raw_data/50000_downsampled_bulk/analysis/"
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
   "lightgreen","cyan4","chartreuse3","bisque3","darksalmon",
   "darkseagreen4","mediumorchid3","lightslategray","red"
  ))
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
     "lightgreen","firebrick3","cyan4","chartreuse3","bisque3","darksalmon",
     "darkseagreen4","mediumorchid3","lightslategray","red"))
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









