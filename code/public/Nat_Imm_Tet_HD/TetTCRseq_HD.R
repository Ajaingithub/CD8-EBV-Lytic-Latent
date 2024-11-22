### Performing analysis for BD Rhapsody using this paper https://www.nature.com/articles/s41590-021-01073-2
library(Seurat)
library(dplyr)
library("sctransform") 
library("glmGamPoi")

raw_data <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/processed_data_download/raw_data/"
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/Analysis/"

GEX_protein <- read.table(paste(raw_data,"Targeted_Gex_protein_expression.txt",sep = ""), header = TRUE)
rownames(GEX_protein) <- GEX_protein$rowname
GEX_protein <- GEX_protein[,-1]

### This is targeted gene expression 
### We start with the GEX
#### GEX ######
GEX <- t(GEX_protein[60:ncol(GEX_protein)])
rownames(GEX) <- gsub("\\.NM.*.|\\.ENS.*.","",rownames(GEX))
rownames(GEX) <-  make.unique(as.character(rownames(GEX)))

Tet_metadata <- read.table(paste(raw_data,"TCR_and_metadata.txt",sep = ""),sep = "\t",  header = TRUE)
rownames(Tet_metadata) <- Tet_metadata$CellID
GEX_obj <- CreateSeuratObject(GEX, meta.data = Tet_metadata)
GEX_obj <- NormalizeData(GEX_obj)
# GEX_obj <- GEX_obj %>%
#   NormalizeData() %>%
#   FindVariableFeatures(nfeatures = 300) %>%
#   ScaleData() %>%
#   RunPCA() %>%
#   RunUMAP(dims = 1:30)
# GEX_obj <- FindNeighbors(GEX_obj, dims = 1:30, verbose = FALSE)
# GEX_obj <- FindClusters(GEX_obj, resolution = 0.6, verbose = FALSE) # resolution 0.4 to 1.2 
# dir.create(paste(savedir,"UMAP/",sep = ""),showWarnings = FALSE)
# pdf(paste(savedir,"UMAP/GEX_umap.pdf",sep = ""), width = 6, height = 5)
# DimPlot(GEX_obj, reduction = "umap")
# dev.off()
# pdf(paste(savedir,"UMAP/GEX_umap_chip.pdf",sep = ""))
# DimPlot(GEX_obj, reduction = "umap", group.by = "chip", split.by = "chip", ncol = 2)
# dev.off()
# pdf(paste(savedir,"UMAP/GEX_umap_donor.pdf",sep = ""))
# DimPlot(GEX_obj, reduction = "umap", group.by = "donor", split.by = "donor", ncol = 5)
# dev.off()

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_sctransform_V2.R")
objname = "GEX"
Assay = "RNA"
process = "sctransform"
GEX_obj_sctransformed <- sctransform_V2_integration(obj = GEX_obj, saveDir = savedir, ngenes = 300,
                                                    regress = c("nCount_RNA"),
                                                    dims = 30,
                                                    Assay = Assay, process = process, objname = objname,
                                                    split_by = "chip",
                                                    reference = NULL,
                                                    sample_tree = NULL)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "GEX"
process = "integration"
Assay = "RNA"
GEX_integrated_RNA <- RNA_integration(GEX_obj_sctransformed, savedir, dims = 20, RNA_features = c("CD3E","CD8A"),
                                      Assay=Assay, process=process, objname=objname, ncol = 4, ndims = 50)

GEX_integrated_RNA_NN <- FindNeighbors(GEX_integrated_RNA, dims = 1:30)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8,1)
res = 0.6
for (i in 1:length(res)) {
  process = paste("RNA_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  GEX_integrated_RNA_NN_cluster <- cluster_UMAP_and_QC(obj_path = GEX_integrated_RNA_NN, dims = 30, res = res[i], saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_RNA","nFeature_RNA"), 
                                                       objname = "GEX_integrated_RNA_NN_cluster",
                                                       process = process, col_sel = c("chip","orig.ident"))
}

DefaultAssay(GEX_integrated_RNA_NN_cluster) <- "RNA"
GEX_integrated_RNA_NN_cluster <- NormalizeData(GEX_integrated_RNA_NN_cluster)
GEX_integrated_RNA_NN_cluster <- ScaleData(GEX_integrated_RNA_NN_cluster, features = rownames(GEX_integrated_RNA_NN_cluster))

dir.create(paste(savedir,"saveRDS_obj",sep = ""))
saveRDS(GEX_integrated_RNA_NN_cluster,paste(savedir,"saveRDS_obj/GEX_integrated_RNA_NN_cluster.RDS",sep = ""))

#### ADT ####
ADT <- t(GEX_protein[1:59])
rownames(ADT) <- gsub(".AH.*.|.Ab.*.","",rownames(ADT))
ADT_obj <- CreateSeuratObject(ADT, assay = "ADT", meta.data = Tet_metadata)
ADT_obj <- NormalizeData(ADT_obj, normalization.method = 'CLR', margin = 2)
ADT_obj <- ScaleData(ADT_obj,features = rownames(ADT_obj))

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")
objname = "obj"
process = "integrated"
Assay = "ADT"

obj_ADT_integrated <- ADT_merging(ADT_obj, savedir, dims = 10, numfeatures=40,
                                  Assay=Assay, process=process, objname=objname,
                                  sample_tree = NULL, split_by = "chip",
                                  reference=NULL, margin = 2)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "obj"
process = "integration"
Assay = "ADT"
obj_ADT_integrated <- RNA_integration(obj_path = obj_ADT_integrated, saveDir = savedir, dims = 10,
                                      RNA_features = c("CD45RO.PTPRC","CD45RA.PTPRC"),Assay = Assay,
                                      process = process, objname = objname, ncol = 2)

### Nearest Neighbour and Clustering
obj_ADT_integrated_NN <- FindNeighbors(obj_ADT_integrated, dims = 1:10)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8)
res = 0.8
for (i in 1:length(res)) {
  process = paste("ADT_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  obj_ADT_integrated_NN_cluster <- cluster_UMAP_and_QC(obj_path = obj_ADT_integrated_NN, dims = 10, res = res[i], saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_ADT","nFeature_ADT"), objname = "obj_ADT_integrated_NN_cluster",
                                                       process = process,col_sel = c("chip","orig.ident"))
}

saveRDS(obj_ADT_integrated_NN_cluster,paste(savedir,"saveRDS_obj/obj_ADT_integrated_NN_cluster.RDS",sep = ""))

#### All Modality Integration ######
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITEseq_modality_integration.R")
# Find the last save object for RNA and ADT
# we have RNA already in the path
library(Seurat)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/Analysis/"
# CD4_ADT_integrated_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD4_ADT_integrated_NN_cluster.RDS",sep = ""))
# CD4_integrated_RNA_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD4_integrated_RNA_NN_cluster.RDS",sep = ""))
RNA_obj_path <- GEX_integrated_RNA_NN_cluster
ADT_obj_path <- obj_ADT_integrated_NN_cluster

objname <- "CD4"
process <- "modality_integrate"
Assay <- "integrated"

obj_modality_integrate <- modality_integration(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path,
                                               RNA_dims = 30, ADT_dims = 10,
                                               saveDir = savedir, Assay = Assay, process = process, objname = objname)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/Modal_cluster_and_UMAP.R")
obj_path <- obj_modality_integrate_cluster
RNA_features = c("CD3E","CD8A")
ADT_features = c("CD45RO.PTPRC","CD45RA.PTPRC")

res = c(0.2,0.4,0.6,0.8,1)
res=1
for (i in 1:length(res)) {
  objname = "CD4"
  Assay = "integrated"
  process = paste("modality_cluster_and_QC",res[i],sep = "_")
  obj_modality_integrate_cluster <- Modal_Cluster_and_UMAP_QC(modal_integrate_obj_path = obj_path, saveDir = savedir,
                                                              res = res[i], RNA_features = RNA_features,
                                                              ADT_features = ADT_features,Assay = Assay, protein_assay = "ADT",
                                                              process = process, objname=objname, 
                                                              col_sel=c("chip","orig.ident","condition"))
}

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/ADT_featureplot.R")
ADT_featureplot(obj = obj_modality_integrate_cluster, savedir = savedir,objname = "obj_modality_integrate_cluster_2",
                process = "featureplots",reduction = "wnn.umap",x = "wnnUMAP_1",y = "wnnUMAP_2", slot="data")

obj_samplewise_counts <- table(obj_modality_integrate_cluster@meta.data$orig.ident) %>% as.data.frame()
write.table(obj_samplewise_counts, paste(savedir,"Table/obj_samplewise_counts.txt",sep = ""), sep = "\t",
            row.names = F, col.names = T, quote = F)

obj_modality_integrate_cluster@meta.data$condition <- gsub("_.*.","",obj_modality_integrate_cluster@meta.data$donor)

pdf(paste(savedir,"UMAP/obj_modality_integrate_cluster_condition.pdf",sep = ""))
DimPlot(obj_modality_integrate_cluster, group.by = "condition", reduction = "wnn.umap")
dev.off()

pdf(paste(savedir,"UMAP/obj_modality_integrate_cluster_condition_splitted.pdf",sep = ""), width = 8, height = 10)
DimPlot(obj_modality_integrate_cluster, group.by = "condition", reduction = "wnn.umap", split.by = "condition", ncol = 2)
dev.off()

obj_modality_integrate_cluster@meta.data$gene.final_condition <- paste(obj_modality_integrate_cluster@meta.data$gene.final,obj_modality_integrate_cluster@meta.data$condition,sep="_")

saveRDS(obj_modality_integrate_cluster,paste(savedir,"saveRDS_obj/obj_modality_integrate_cluster_1.RDS",sep = ""))

### featureplot for each gene
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
DefaultAssay(obj_modality_integrate_cluster) <- "RNA"
genes <- rownames(obj_modality_integrate_cluster)
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p <- featureplot_front(obj_modality_integrate_cluster, genes[i],
                         reduction = "wnn.umap", x="wnnUMAP_1", y="wnnUMAP_2",size=0.001) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- p
}

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/obj_All_genes_unimpute_1.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list[1:100])
dev.off()

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/obj_All_genes_unimpute_2.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list[101:205])
dev.off()

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/obj_All_genes_unimpute_3.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list[206:310])
dev.off()

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/obj_All_genes_unimpute_4.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list[311:418])
dev.off()


source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
DefaultAssay(obj_modality_integrate_cluster) <- "MAGIC_RNA"
genes <- rownames(obj_modality_integrate_cluster)
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p <- featureplot_front(obj_modality_integrate_cluster, genes[i],
                         reduction = "wnn.umap", x="wnnUMAP_1", y="wnnUMAP_2",size=0.001) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- p
}

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/obj_All_genes_impute_1.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list[1:100])
dev.off()

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/obj_All_genes_impute_2.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list[101:205])
dev.off()

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/obj_All_genes_impute_3.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list[206:310])
dev.off()

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/obj_All_genes_impute_4.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list[311:418])
dev.off()

obj_modality_integrate_cluster@meta.data$condition_gene_specificity <- paste(obj_modality_integrate_cluster@meta.data$condition,
                                                                         obj_modality_integrate_cluster@meta.data$gene.final,
                                                                         obj_modality_integrate_cluster@meta.data$specificity.final, sep = "_")


antigens <- obj_modality_integrate_cluster@meta.data[grep("^EBV$",obj_modality_integrate_cluster@meta.data$gene.final),"specificity.final"]  %>% unique() %>% grep("undetermined",.,value=TRUE,invert=TRUE)

for (i in 1:length(antigens)) {
  healthy_cellnames <- rownames(obj_modality_integrate_cluster@meta.data[grep(paste("healthy_EBV_",antigens[i],sep = ""),
                                                                      obj_modality_integrate_cluster@meta.data$condition_gene_specificity),])
  p1 <- DimPlot(obj_modality_integrate_cluster,cells.highlight = healthy_cellnames,pt.size = 0.05, sizes.highlight = 0.5, reduction = "wnn.umap") + ggtitle(paste("healthy_EBV_",antigens[i],sep = "")) +
    NoLegend() + 
    theme(plot.title = element_text(hjust = 0.5))
  
  T1D_cellnames <- rownames(obj_modality_integrate_cluster@meta.data[grep(paste("T1D_EBV_",antigens[i],sep = ""),
                                                                          obj_modality_integrate_cluster@meta.data$condition_gene_specificity),])
  p2 <- DimPlot(obj_modality_integrate_cluster,cells.highlight = T1D_cellnames,sizes.highlight = 0.5, pt.size = 0.05, reduction = "wnn.umap") + 
    ggtitle(paste("T1D_EBV_",antigens[i],sep = "")) +
    NoLegend() + 
    theme(plot.title = element_text(hjust = 0.5))
  pdf(paste(savedir,"featureplot/",antigens[i],"_healthy_T1D.pdf",sep = ""), width = 10, height = 5)
  print(p1+p2)
  dev.off()
}


### Extracting the EBV healthy individuals #####
EBV_cells <- obj_modality_integrate_cluster@meta.data[grep("^EBV_healthy$",obj_modality_integrate_cluster@meta.data$gene.final_condition),] %>% rownames()
GEX_protein <- read.table(paste(raw_data,"Targeted_Gex_protein_expression.txt",sep = ""), header = TRUE)
rownames(GEX_protein) <- GEX_protein$rowname
GEX_protein <- GEX_protein[,-1]
GEX_protein_EBV <- GEX_protein[match(EBV_cells, rownames(GEX_protein)),]

Tet_metadata <- read.table(paste(raw_data,"TCR_and_metadata.txt",sep = ""),sep = "\t",  header = TRUE)
rownames(Tet_metadata) <- Tet_metadata$CellID
Tet_metadata_EBV <- Tet_metadata[match(EBV_cells, rownames(Tet_metadata)),]

#### GEX EBV ######
GEX_EBV <- t(GEX_protein_EBV[60:ncol(GEX_protein_EBV)])
rownames(GEX_EBV) <- gsub("\\.NM.*.|\\.ENS.*.","",rownames(GEX_EBV))
rownames(GEX_EBV) <-  make.unique(as.character(rownames(GEX_EBV)))
all(rownames(Tet_metadata_EBV)== colnames(GEX_EBV))

GEX_EBV_obj <- CreateSeuratObject(GEX_EBV, meta.data = Tet_metadata_EBV)
GEX_EBV_obj <- NormalizeData(GEX_EBV_obj)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_sctransform_V2.R")
objname = "GEX"
Assay = "RNA"
process = "sctransform"
GEX_EBV_obj_sctransformed <- sctransform_V2_integration(obj = GEX_EBV_obj, saveDir = savedir, ngenes = 300,
                                                    regress = c("nCount_RNA"),
                                                    dims = 30,
                                                    Assay = Assay, process = process, objname = objname,
                                                    split_by = "chip",
                                                    reference = NULL,
                                                    sample_tree = NULL,
                                                    k_weight=55, 
                                                    k_filter = 55)

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_analysis/"
dir.create(savedir, showWarnings = FALSE)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "GEX_EBV"
process = "integration"
Assay = "RNA"
GEX_EBV_integrated_RNA <- RNA_integration(GEX_EBV_obj_sctransformed, savedir, dims = 20, RNA_features = c("CD3E","CD8A"),
                                      Assay=Assay, process=process, objname=objname, ncol = 3, ndims = 50)

GEX_EBV_integrated_RNA_NN <- FindNeighbors(GEX_EBV_integrated_RNA, dims = 1:30)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8,1)
res = 0.6
for (i in 1:length(res)) {
  process = paste("RNA_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  GEX_EBV_integrated_RNA_NN_cluster <- cluster_UMAP_and_QC(obj_path = GEX_EBV_integrated_RNA_NN, dims = 30, res = res[i], 
                                                           saveDir=savedir, Assay = Assay,
                                                           QC_features = c("nCount_RNA","nFeature_RNA"),
                                                           objname = "GEX_EBV_integrated_RNA_NN_cluster",
                                                           process = process, col_sel = c("chip","orig.ident"))
}

DefaultAssay(GEX_EBV_integrated_RNA_NN_cluster) <- "RNA"
GEX_EBV_integrated_RNA_NN_cluster <- NormalizeData(GEX_EBV_integrated_RNA_NN_cluster)
GEX_EBV_integrated_RNA_NN_cluster <- ScaleData(GEX_EBV_integrated_RNA_NN_cluster, features = rownames(GEX_EBV_integrated_RNA_NN_cluster))

dir.create(paste(savedir,"saveRDS_obj",sep = ""))
saveRDS(GEX_EBV_integrated_RNA_NN_cluster,paste(savedir,"saveRDS_obj/GEX_EBV_integrated_RNA_NN_cluster.RDS",sep = ""))

#### ADT EBV ####
ADT_EBV <- t(GEX_protein_EBV[1:59])
rownames(ADT_EBV) <- gsub(".AH.*.|.Ab.*.","",rownames(ADT_EBV))
all(colnames(ADT_EBV) == rownames(Tet_metadata_EBV))

ADT_EBV_obj <- CreateSeuratObject(ADT_EBV, assay = "ADT", meta.data = Tet_metadata_EBV)
ADT_EBV_obj <- NormalizeData(ADT_EBV_obj, normalization.method = 'CLR', margin = 2)
ADT_EBV_obj <- ScaleData(ADT_EBV_obj,features = rownames(ADT_EBV_obj))

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")
objname = "obj_EBV"
process = "integrated"
Assay = "ADT"
obj_ADT_EBV_integrated <- ADT_merging(ADT_EBV_obj, savedir, dims = 10, numfeatures=40,
                                  Assay=Assay, process=process, objname=objname,
                                  sample_tree = NULL, split_by = "chip",
                                  reference=NULL, margin = 2, k_filter = 55, k_weight = 40)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "obj_EBV"
process = "integration"
Assay = "ADT"
obj_ADT_EBV_integrated <- RNA_integration(obj_path = obj_ADT_EBV_integrated, saveDir = savedir, dims = 10,
                                      RNA_features = c("CD45RO.PTPRC","CD45RA.PTPRC"),Assay = Assay,
                                      process = process, objname = objname, ncol = 2)

### Nearest Neighbour and Clustering
obj_ADT_EBV_integrated_NN <- FindNeighbors(obj_ADT_EBV_integrated, dims = 1:10)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8)
res = 0.8
for (i in 1:length(res)) {
  process = paste("ADT_EBV_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  obj_ADT_EBV_integrated_NN_cluster <- cluster_UMAP_and_QC(obj_path = obj_ADT_EBV_integrated_NN, dims = 10, res = res[i], saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_ADT","nFeature_ADT"), objname = "obj_ADT_EBV_integrated_NN_cluster",
                                                       process = process,col_sel = c("chip","orig.ident"))
}

saveRDS(obj_ADT_EBV_integrated_NN_cluster,paste(savedir,"saveRDS_obj/obj_ADT_EBV_integrated_NN_cluster.RDS",sep = ""))

#### EBV Modality Integration ######
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITEseq_modality_integration.R")
# Find the last save object for RNA and ADT
# we have RNA already in the path
library(Seurat)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_analysis/"
# CD4_ADT_integrated_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD4_ADT_integrated_NN_cluster.RDS",sep = ""))
# CD4_integrated_RNA_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD4_integrated_RNA_NN_cluster.RDS",sep = ""))
RNA_obj_path <- GEX_EBV_integrated_RNA_NN_cluster
ADT_obj_path <- obj_ADT_EBV_integrated_NN_cluster

objname <- "obj_EBV"
process <- "modality_integrate"
Assay <- "integrated"

obj_EBV_modality_integrate <- modality_integration(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path,
                                               RNA_dims = 20, ADT_dims = 10,
                                               saveDir = savedir, Assay = Assay, process = process, objname = objname)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/Modal_cluster_and_UMAP.R")
obj_path <- obj_EBV_modality_integrate
RNA_features = c("CD3E","CD8A")
ADT_features = c("CD45RO.PTPRC","CD45RA.PTPRC")

res = c(0.2,0.4,0.6,0.8,1)
res=0.8
for (i in 1:length(res)) {
  objname = "obj_EBV"
  Assay = "integrated"
  process = paste("modality_cluster_and_QC",res[i],sep = "_")
  obj_EBV_modality_integrate_cluster <- Modal_Cluster_and_UMAP_QC(modal_integrate_obj_path = obj_path, saveDir = savedir,
                                                              res = res[i], RNA_features = RNA_features,
                                                              ADT_features = ADT_features,Assay = Assay, protein_assay = "ADT",
                                                              process = process, objname=objname,
                                                              col_sel=c("chip","orig.ident","specificity.final","donor"))
}

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/ADT_featureplot.R")
ADT_featureplot(obj = obj_EBV_modality_integrate_cluster, savedir = savedir,objname = "obj_EBV_modality_integrate_cluster_2",
                process = "featureplots",reduction = "wnn.umap",x = "wnnUMAP_1",y = "wnnUMAP_2", slot="data")

obj_EBV_samplewise_counts <- table(obj_EBV_modality_integrate_cluster@meta.data$orig.ident) %>% as.data.frame()
write.table(obj_EBV_samplewise_counts, paste(savedir,"Table/obj_EBV_samplewise_counts.txt",sep = ""), sep = "\t",
            row.names = F, col.names = T, quote = F)


pdf(paste(savedir,"UMAP/obj_EBV_modality_integrate_cluster_antigen.pdf",sep = ""))
DimPlot(obj_EBV_modality_integrate_cluster, group.by = "specificity.final", reduction = "wnn.umap")
dev.off()

pdf(paste(savedir,"UMAP/obj_EBV_modality_integrate_cluster_antigen_splitted.pdf",sep = ""), width = 14, height = 10)
DimPlot(obj_EBV_modality_integrate_cluster, group.by = "specificity.final", reduction = "wnn.umap", split.by = "specificity.final", ncol = 4)
dev.off()

saveRDS(obj_EBV_modality_integrate_cluster,paste(savedir,"saveRDS_obj/obj_EBV_modality_integrate_cluster_0.8.RDS",sep = ""))

### featureplot for each gene
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
DefaultAssay(obj_EBV_modality_integrate_cluster) <- "RNA"
genes <- rownames(obj_EBV_modality_integrate_cluster)
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p <- featureplot_front(obj_EBV_modality_integrate_cluster, genes[i],
                         reduction = "wnn.umap", x="wnnUMAP_1", y="wnnUMAP_2",size=0.5) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- p
}

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/All_genes_unimpute.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list)
dev.off()


source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
DefaultAssay(obj_EBV_modality_integrate_cluster) <- "MAGIC_RNA"
genes <- rownames(obj_EBV_modality_integrate_cluster)
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p <- featureplot_front(obj_EBV_modality_integrate_cluster, genes[i],
                         reduction = "wnn.umap", x="wnnUMAP_1", y="wnnUMAP_2",size=0.5) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- p
}

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/All_genes_imputed.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list)
dev.off()

#### EBV T1D and Healthy #####
library(Seurat)
library(dplyr)
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/"
dir.create(savedir, showWarnings = FALSE)
obj_modality_integrate_cluster <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/Analysis/saveRDS_obj/obj_modality_integrate_cluster_1.RDS")
raw_data <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/processed_data_download/raw_data/"
obj_modality_integrate_cluster@meta.data$condition <- gsub("_.*.","",obj_modality_integrate_cluster@meta.data$donor)
obj_modality_integrate_cluster@meta.data$gene.final_condition <- paste(obj_modality_integrate_cluster@meta.data$gene.final , obj_modality_integrate_cluster@meta.data$condition, sep = "_")
EBV_cells <- obj_modality_integrate_cluster@meta.data[grep("^EBV_T1D$|^EBV_healthy",obj_modality_integrate_cluster@meta.data$gene.final_condition),] %>% rownames()

GEX_protein <- read.table(paste(raw_data,"Targeted_Gex_protein_expression.txt",sep = ""), header = TRUE)
rownames(GEX_protein) <- GEX_protein$rowname
GEX_protein <- GEX_protein[,-1]
GEX_protein_EBV <- GEX_protein[match(EBV_cells, rownames(GEX_protein)),]

Tet_metadata <- read.table(paste(raw_data,"TCR_and_metadata.txt",sep = ""),sep = "\t",  header = TRUE)
rownames(Tet_metadata) <- Tet_metadata$CellID
Tet_metadata_EBV <- Tet_metadata[match(EBV_cells, rownames(Tet_metadata)),]
Tet_metadata_EBV$condition <- gsub("_.*.","",Tet_metadata_EBV$donor)
Tet_metadata_EBV$gene.final_condition <- paste(Tet_metadata_EBV$gene.final, Tet_metadata_EBV$condition, sep = "_")

#### GEX EBV ######
GEX_EBV <- t(GEX_protein_EBV[60:ncol(GEX_protein_EBV)])
rownames(GEX_EBV) <- gsub("\\.NM.*.|\\.ENS.*.","",rownames(GEX_EBV))
rownames(GEX_EBV) <-  make.unique(as.character(rownames(GEX_EBV)))
all(rownames(Tet_metadata_EBV)== colnames(GEX_EBV))

GEX_EBV_obj <- CreateSeuratObject(GEX_EBV, meta.data = Tet_metadata_EBV)
GEX_EBV_obj <- NormalizeData(GEX_EBV_obj)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_sctransform_V2.R")
objname = "GEX"
Assay = "RNA"
process = "sctransform"
GEX_EBV_obj_sctransformed <- sctransform_V2_integration(obj = GEX_EBV_obj, saveDir = savedir, ngenes = 300,
                                                        regress = c("nCount_RNA"),
                                                        dims = 30,
                                                        Assay = Assay, process = process, objname = objname,
                                                        split_by = "chip",
                                                        reference = NULL,
                                                        sample_tree = NULL)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "GEX_EBV"
process = "integration"
Assay = "RNA"
GEX_EBV_integrated_RNA <- RNA_integration(GEX_EBV_obj_sctransformed, savedir, dims = 20, RNA_features = c("CD3E","CD8A"),
                                          Assay=Assay, process=process, objname=objname, ncol = 3, ndims = 50)

GEX_EBV_integrated_RNA_NN <- FindNeighbors(GEX_EBV_integrated_RNA, dims = 1:30)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8,1)
res = 0.6
for (i in 1:length(res)) {
  process = paste("RNA_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  GEX_EBV_integrated_RNA_NN_cluster <- cluster_UMAP_and_QC(obj_path = GEX_EBV_integrated_RNA_NN, dims = 20, res = res[i], 
                                                           saveDir=savedir, Assay = Assay,
                                                           QC_features = c("nCount_RNA","nFeature_RNA"),
                                                           objname = "GEX_EBV_integrated_RNA_NN_cluster",
                                                           process = process, col_sel = c("chip","orig.ident"))
}

DefaultAssay(GEX_EBV_integrated_RNA_NN_cluster) <- "RNA"
GEX_EBV_integrated_RNA_NN_cluster <- NormalizeData(GEX_EBV_integrated_RNA_NN_cluster)
GEX_EBV_integrated_RNA_NN_cluster <- ScaleData(GEX_EBV_integrated_RNA_NN_cluster, features = rownames(GEX_EBV_integrated_RNA_NN_cluster))

dir.create(paste(savedir,"saveRDS_obj",sep = ""))
saveRDS(GEX_EBV_integrated_RNA_NN_cluster,paste(savedir,"saveRDS_obj/GEX_EBV_integrated_RNA_NN_cluster.RDS",sep = ""))

#### ADT EBV ####
ADT_EBV <- t(GEX_protein_EBV[1:59])
rownames(ADT_EBV) <- gsub(".AH.*.|.Ab.*.","",rownames(ADT_EBV))
all(colnames(ADT_EBV) == rownames(Tet_metadata_EBV))

ADT_EBV_obj <- CreateSeuratObject(ADT_EBV, assay = "ADT", meta.data = Tet_metadata_EBV)
ADT_EBV_obj <- NormalizeData(ADT_EBV_obj, normalization.method = 'CLR', margin = 2)
ADT_EBV_obj <- ScaleData(ADT_EBV_obj,features = rownames(ADT_EBV_obj))

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")
objname = "obj_EBV"
process = "integrated"
Assay = "ADT"
obj_ADT_EBV_integrated <- ADT_merging(ADT_EBV_obj, savedir, dims = 10, numfeatures=40,
                                      Assay=Assay, process=process, objname=objname,
                                      sample_tree = NULL, split_by = "chip",
                                      reference=NULL, margin = 2)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "obj_EBV"
process = "integration"
Assay = "ADT"
obj_ADT_EBV_integrated <- RNA_integration(obj_path = obj_ADT_EBV_integrated, saveDir = savedir, dims = 10,
                                          RNA_features = c("CD45RO.PTPRC","CD45RA.PTPRC"),Assay = Assay,
                                          process = process, objname = objname, ncol = 2)

### Nearest Neighbour and Clustering
obj_ADT_EBV_integrated_NN <- FindNeighbors(obj_ADT_EBV_integrated, dims = 1:10)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8)
res = 0.8
for (i in 1:length(res)) {
  process = paste("ADT_EBV_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  obj_ADT_EBV_integrated_NN_cluster <- cluster_UMAP_and_QC(obj_path = obj_ADT_EBV_integrated_NN, dims = 10, res = res[i], saveDir=savedir, Assay = Assay,
                                                           QC_features = c("nCount_ADT","nFeature_ADT"), objname = "obj_ADT_EBV_integrated_NN_cluster",
                                                           process = process,col_sel = c("chip","orig.ident"))
}

saveRDS(obj_ADT_EBV_integrated_NN_cluster,paste(savedir,"saveRDS_obj/obj_ADT_EBV_integrated_NN_cluster.RDS",sep = ""))

#### EBV Modality Integration ######
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITEseq_modality_integration.R")
# Find the last save object for RNA and ADT
# we have RNA already in the path
library(Seurat)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/"
# savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_analysis/"
# CD4_ADT_integrated_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD4_ADT_integrated_NN_cluster.RDS",sep = ""))
# CD4_integrated_RNA_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD4_integrated_RNA_NN_cluster.RDS",sep = ""))
RNA_obj_path <- GEX_EBV_integrated_RNA_NN_cluster
ADT_obj_path <- obj_ADT_EBV_integrated_NN_cluster

objname <- "obj_EBV"
process <- "modality_integrate"
Assay <- "integrated"

obj_EBV_modality_integrate <- modality_integration(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path,
                                                   RNA_dims = 20, ADT_dims = 10,
                                                   saveDir = savedir, Assay = Assay, process = process, objname = objname)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/Modal_cluster_and_UMAP.R")
obj_path <- obj_EBV_modality_integrate
RNA_features = c("CD3E","CD8A")
ADT_features = c("CD45RO.PTPRC","CD45RA.PTPRC")

res = c(0.2,0.4,0.6,0.8,1)
res=0.8
for (i in 1:length(res)) {
  objname = "obj_EBV"
  Assay = "integrated"
  process = paste("modality_cluster_and_QC",res[i],sep = "_")
  obj_EBV_modality_integrate_cluster <- Modal_Cluster_and_UMAP_QC(modal_integrate_obj_path = obj_path, saveDir = savedir,
                                                                  res = res[i], RNA_features = RNA_features,
                                                                  ADT_features = ADT_features,Assay = Assay, protein_assay = "ADT",
                                                                  process = process, objname=objname,
                                                                  col_sel=c("chip","orig.ident","specificity.final","donor","condition"))
}

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/ADT_featureplot.R")
ADT_featureplot(obj = obj_EBV_modality_integrate_cluster, savedir = savedir,objname = "obj_EBV_modality_integrate_cluster_2",
                process = "featureplots",reduction = "wnn.umap",x = "wnnUMAP_1",y = "wnnUMAP_2", slot="data")

obj_EBV_samplewise_counts <- table(obj_EBV_modality_integrate_cluster@meta.data$orig.ident) %>% as.data.frame()
write.table(obj_EBV_samplewise_counts, paste(savedir,"Table/obj_EBV_chipwise_counts.txt",sep = ""), sep = "\t",
            row.names = F, col.names = T, quote = F)

obj_EBV_samplewise_counts <- table(obj_EBV_modality_integrate_cluster@meta.data$donor) %>% as.data.frame()
write.table(obj_EBV_samplewise_counts, paste(savedir,"Table/obj_EBV_samplewise_counts.txt",sep = ""), sep = "\t",
            row.names = F, col.names = T, quote = F)

obj_EBV_samplewise_counts <- table(obj_EBV_modality_integrate_cluster@meta.data$donor,obj_EBV_modality_integrate_cluster@meta.data$seurat_clusters) 
write.table(obj_EBV_samplewise_counts, paste(savedir,"Table/obj_EBV_samplewise_counts_seurat_cluster.txt",sep = ""), sep = "\t",
            row.names = T, col.names = T, quote = F)

obj_EBV_modality_integrate_cluster@meta.data$EBV_Ag_condition <- paste(obj_EBV_modality_integrate_cluster@meta.data$specificity.final, obj_EBV_modality_integrate_cluster@meta.data$condition,sep = "_")

pdf(paste(savedir,"UMAP/obj_EBV_modality_integrate_cluster_antigen.pdf",sep = ""), width = 18)
DimPlot(obj_EBV_modality_integrate_cluster, group.by = "EBV_Ag_condition", reduction = "wnn.umap") 
dev.off()

pdf(paste(savedir,"UMAP/obj_EBV_modality_integrate_cluster_antigen_splitted.pdf",sep = ""), width = 8, height = 12)
DimPlot(obj_EBV_modality_integrate_cluster, group.by = "EBV_Ag_condition", reduction = "wnn.umap", 
        split.by = "EBV_Ag_condition", ncol = 4) + NoLegend()
dev.off()

obj_EBV_modality_integrate_cluster@meta.data$donor_Ag_specific <- paste(obj_EBV_modality_integrate_cluster@meta.data$donor, obj_EBV_modality_integrate_cluster@meta.data$specificity.final,sep="_")
obj_EBV_samplewise_counts <- table(obj_EBV_modality_integrate_cluster@meta.data$donor_Ag_specific,obj_EBV_modality_integrate_cluster@meta.data$seurat_clusters) 
write.table(obj_EBV_samplewise_counts, paste(savedir,"Table/obj_EBV_Ag_samplewise_counts_seurat_cluster.txt",sep = ""), sep = "\t",
            row.names = T, col.names = T, quote = F)

saveRDS(obj_EBV_modality_integrate_cluster,paste(savedir,"saveRDS_obj/obj_EBV_modality_integrate_cluster_0.8.RDS",sep = ""))

### featureplot for each gene
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
DefaultAssay(obj_EBV_modality_integrate_cluster) <- "RNA"
genes <- rownames(obj_EBV_modality_integrate_cluster)
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p <- featureplot_front(obj_EBV_modality_integrate_cluster, genes[i],
                         reduction = "wnn.umap", x="wnnUMAP_1", y="wnnUMAP_2",size=0.5) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- p
}

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/All_genes_unimpute.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list)
dev.off()


source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
DefaultAssay(obj_EBV_modality_integrate_cluster) <- "MAGIC_RNA"
genes <- rownames(obj_EBV_modality_integrate_cluster)
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p <- featureplot_front(obj_EBV_modality_integrate_cluster, genes[i],
                         reduction = "wnn.umap", x="wnnUMAP_1", y="wnnUMAP_2",size=0.5) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- p
}

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/All_genes_imputed.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list)
dev.off()


antigens <- obj_EBV_modality_integrate_cluster@meta.data[grep("^EBV$",obj_EBV_modality_integrate_cluster@meta.data$gene.final),"specificity.final"]  %>% unique() %>% grep("undetermined",.,value=TRUE,invert=TRUE)

for (i in 1:length(antigens)) {
  healthy_cellnames <- rownames(obj_EBV_modality_integrate_cluster@meta.data[grep(paste(antigens[i],"_healthy",sep = ""),
                                                                              obj_EBV_modality_integrate_cluster@meta.data$EBV_Ag_condition),])
  p1 <- DimPlot(obj_EBV_modality_integrate_cluster,cells.highlight = healthy_cellnames,pt.size = 0.05, sizes.highlight = 0.5, reduction = "wnn.umap") + ggtitle(paste("healthy_EBV_",antigens[i],sep = "")) +
    NoLegend() + 
    theme(plot.title = element_text(hjust = 0.5))
  
  T1D_cellnames <- rownames(obj_EBV_modality_integrate_cluster@meta.data[grep(paste(antigens[i],"_T1D",sep = ""),
                                                                          obj_EBV_modality_integrate_cluster@meta.data$EBV_Ag_condition),])
  p2 <- DimPlot(obj_EBV_modality_integrate_cluster,cells.highlight = T1D_cellnames,sizes.highlight = 0.5, pt.size = 0.05, reduction = "wnn.umap") + 
    ggtitle(paste("T1D_EBV_",antigens[i],sep = "")) +
    NoLegend() + 
    theme(plot.title = element_text(hjust = 0.5))
  pdf(paste(savedir,"featureplot/",antigens[i],"_healthy_T1D.pdf",sep = ""), width = 10, height = 5)
  print(p1+p2)
  dev.off()
}

### EBV Combined ####
library(Seurat)
obj_EBV_modality_integrate_cluster <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/saveRDS_obj/obj_EBV_modality_integrate_cluster_0.8.RDS")
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/"

library(stringr)
patterns = c(10,2,8,0,1,9,3,4,5,6,7,11,12)
replacements = c("Naive","CM/EM","CM/EM","EM","EM","EM","TEMRA","TEMRA","TEMRA","TEMRA","TEMRA","unknown","unknown")
# names(replacement) <- patterns
obj_EBV_modality_integrate_cluster@meta.data$Celltypes <- obj_EBV_modality_integrate_cluster@meta.data$seurat_clusters
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  obj_EBV_modality_integrate_cluster@meta.data$Celltypes <- str_replace_all(obj_EBV_modality_integrate_cluster@meta.data$Celltypes, pattern, replacements[i])
}

UMAP <- obj_EBV_modality_integrate_cluster@reductions$wnn.umap@cell.embeddings %>% as.data.frame()

all(rownames(UMAP) == rownames(obj_EBV_modality_integrate_cluster@meta.data))
UMAP$Cluster <- obj_EBV_modality_integrate_cluster@meta.data$seurat_clusters
UMAP$Celltypes <- obj_EBV_modality_integrate_cluster@meta.data$Celltypes
UMAP$Antigen <- obj_EBV_modality_integrate_cluster@meta.data$specificity.final
library(stringr)
patterns = c(0:12)
replacements = c("steelblue1","skyblue1","limegreen","firebrick1","tomato","red2","darkred","red4",
                 "seagreen3","royalblue","goldenrod1","blue1","gray80")
# names(replacement) <- patterns
UMAP$cluster_color <- UMAP$Cluster
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  UMAP$cluster_color <- str_replace_all(UMAP$cluster_color, pattern, replacements[i])
}
col <- as.character(UMAP$cluster_color)
names(col) <- as.character(UMAP$Cluster)

Ag <- grep("undetermined",unique(UMAP$Antigen),invert=TRUE,value=TRUE)
for (i in 1:length(Ag)) {
  UMAP_sub <- UMAP[grep(paste("^",Ag[i],"$",sep=""),UMAP$Antigen),]
  p1 <- ggplot(UMAP_sub,aes(x=wnnUMAP_1,y=wnnUMAP_2,label=cluster)) +
    geom_point(aes(colour = cluster), size=1, show.legend = TRUE) +
    scale_color_manual(values=col) +
    xlim(min(UMAP$wnnUMAP_1), max(UMAP$wnnUMAP_1)) +
    ylim(min(UMAP$wnnUMAP_2), max(UMAP$wnnUMAP_2)) +
    # NoLegend() +
    # ggrepel::geom_label_repel(data = label.df_2, aes(label = cluster), size=5) +
    ggtitle(paste(Ag[i])) +
    theme(panel.background = element_rect(fill = 'white',color='black'),
          plot.title = element_text(hjust = 0.5)) +
    
  pdf(paste(savedir,"UMAP/",Ag[i],"_subset.pdf",sep = ""), width = 4.3, height = 4)
  print(p1)
  dev.off()
}

library(stringr)
patterns = c("CM/EM","EM","Naive","TEMRA","unknown")
replacements = c("limegreen","skyblue1","goldenrod1","darkred","gray80")
# names(replacement) <- patterns
UMAP$celltypes_color <- UMAP$celltypes
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  UMAP$celltypes_color <- str_replace_all(UMAP$celltypes_color, pattern, replacements[i])
}
col <- as.character(UMAP$celltypes_color)
names(col) <- as.character(UMAP$celltypes)

Ag <- grep("undetermined",unique(UMAP$Antigen),invert=TRUE,value=TRUE)
for (i in 1:length(Ag)) {
  UMAP_sub <- UMAP[grep(paste("^",Ag[i],"$",sep=""),UMAP$Antigen),]
  p1 <- ggplot(UMAP_sub,aes(x=wnnUMAP_1,y=wnnUMAP_2,label=celltypes)) +
    geom_point(aes(colour = celltypes), size=1, show.legend = TRUE) +
    scale_color_manual(values=col) +
    xlim(min(UMAP$wnnUMAP_1), max(UMAP$wnnUMAP_1)) +
    ylim(min(UMAP$wnnUMAP_2), max(UMAP$wnnUMAP_2)) +
    # NoLegend() +
    # ggrepel::geom_label_repel(data = label.df_2, aes(label = cluster), size=5) +
    ggtitle(paste(Ag[i])) +
    theme(panel.background = element_rect(fill = 'white',color='black'),
          plot.title = element_text(hjust = 0.5)) +
    
    pdf(paste(savedir,"UMAP/",Ag[i],"_subset_celltypes.pdf",sep = ""), width = 4.7, height = 4)
  print(p1)
  dev.off()
}

# obj_EBV_modality_integrate_cluster_2 <- obj_EBV_modality_integrate_cluster
# Idents(obj_EBV_modality_integrate_cluster_2) <- obj_EBV_modality_integrate_cluster@meta.data$specificity.final
# Ag <- grep("undetermined",unique(obj_EBV_modality_integrate_cluster_2@meta.data$specificity.final),invert=TRUE,value=TRUE)
# for (i in 1:length(Ag)) {
#   obj <- subset(obj_EBV_modality_integrate_cluster_2, idents = c(Ag[i]))
#   pdf(paste(savedir,"UMAP/",Ag[i],"_subset.pdf",sep = ""), width = 4.3, height = 4)
#   print(DimPlot(obj, group.by = "seurat_clusters", reduction = "wnn.umap") + geom_point(aes(color = col) + scale_color_manual(values=col)))
#   dev.off()
# }
# 
# obj_EBV_modality_integrate_cluster_2 <- obj_EBV_modality_integrate_cluster
# Idents(obj_EBV_modality_integrate_cluster_2) <- obj_EBV_modality_integrate_cluster@meta.data$specificity.final
# Ag <- grep("undetermined",unique(obj_EBV_modality_integrate_cluster_2@meta.data$specificity.final),invert=TRUE,value=TRUE)
# for (i in 1:length(Ag)) {
#   obj <- subset(obj_EBV_modality_integrate_cluster_2, idents = c(Ag[i]))
#   
#   pdf(paste(savedir,"UMAP/",Ag[i],"_subset_celltypes.pdf",sep = ""),width = 4.7, height = 4)
#   print(DimPlot(obj, group.by = "Celltypes", reduction = "wnn.umap", 
#                 cols = c("limegreen","skyblue1","goldenrod1","darkred","gray80")))
#   dev.off()
# }

saveRDS(obj_EBV_modality_integrate_cluster,paste(savedir,"saveRDS_obj/obj_EBV_modality_integrate_cluster_0.8_2.RDS",sep = ""))

### TCR alpha and beta UMAP ####
### CD4
# library(taRifx)
library(ggplot2)
library(Seurat)
library(dplyr)
library(Seurat)
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/TCR/"
dir.create(savedir,showWarnings = FALSE)
setwd(savedir)

### In here we will do TRA CDR3 and TRB CDR3 separately
TRB_CDR3 <- table(obj_EBV_modality_integrate_cluster@meta.data$TCR) %>% as.data.frame()
TRB_clonal <- TRB_CDR3[order(-TRB_CDR3$Freq),]
TRB_clonal_gt_2 <- TRB_clonal[TRB_clonal$Freq > 1,]
TRB_clonal_gt_2$Var1 <- factor(TRB_clonal_gt_2$Var1, levels = TRB_clonal_gt_2$Var1)
TRB_clonal_gt_2 <- TRB_clonal_gt_2[grep("^_$",TRB_clonal_gt_2$Var1,invert = TRUE),]

pdf(paste(savedir,"TRB_clonal_expanded_gt_1.pdf",sep = ""), width = 70, height = 6)
ggplot(TRB_clonal_gt_2, aes(x=Var1,y=Freq, fill = Var1)) + geom_bar(stat="identity", color="black") +  
  ggtitle("obj_EBV_modality_integrate_cluster TRB CDR3 aa gt 1") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey")) + NoLegend()
dev.off()

rm(plot_list)
plot_list <- list()
for(i in 1:length(TRB_clonal_gt_2$Var1)){
  clonal_expand <- rownames(obj_EBV_modality_integrate_cluster@meta.data[grep(TRB_clonal_gt_2$Var1[i],obj_EBV_modality_integrate_cluster@meta.data$TCR),])
  p <- DimPlot(obj_EBV_modality_integrate_cluster,cells.highlight = clonal_expand, reduction = "wnn.umap", label = TRUE, label.size = 3, sizes.highlight = 0.5) + 
    ggtitle(TRB_clonal_gt_2$Var1[i]) + 
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5))
  # pdf(paste(savedir,"obj_EBV_modality_integrate_cluster/obj_EBV_modality_integrate_cluster_TRB_",TRB_clonal_gt_2$Var1[i],".pdf",sep = ""))
  # print(p)
  # dev.off()
  plot_list[[i]] <- p
}

pdf(paste(savedir,"obj_EBV_modality_integrate_cluster_TCR_all_UMAP.pdf",sep = ""), width = 4.3, height = 4)
plot_list
dev.off()

## TRB greater than 1
clonal_expand <- rownames(obj_EBV_modality_integrate_cluster@meta.data[grep(paste(TRB_clonal_gt_2$Var1,collapse = "|"),obj_EBV_modality_integrate_cluster@meta.data$TCR),])
p <- DimPlot(obj_EBV_modality_integrate_cluster,cells.highlight = clonal_expand, reduction = "wnn.umap", label = TRUE, sizes.highlight = 0.3) + ggtitle("TRB clones gt 1") + 
  NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5))
pdf(paste(savedir,"obj_EBV_modality_integrate_cluster_TCR_all_UMAP_combined.pdf",sep = ""), width = 5, height = 5)
print(p)
dev.off()

n = levels(obj_EBV_modality_integrate_cluster@meta.data$seurat_clusters)
df <- data.frame(matrix(ncol = length(n), nrow = length(n)))
colnames(df) <- paste("C",n,sep = "")
rownames(df) <- paste("C",n,sep = "")

obj_EBV_modality_integrate_cluster_TCR_cluster <- table(obj_EBV_modality_integrate_cluster@meta.data$TCR, obj_EBV_modality_integrate_cluster@meta.data$seurat_clusters)
# head(obj_EBV_modality_integrate_cluster_TCR_cluster)
obj_EBV_modality_integrate_cluster_TCR_cluster <- obj_EBV_modality_integrate_cluster_TCR_cluster[grep("^_$",rownames(obj_EBV_modality_integrate_cluster_TCR_cluster),invert=TRUE),]
obj_EBV_modality_integrate_cluster_TCR_cluster[obj_EBV_modality_integrate_cluster_TCR_cluster > 0] <- 1
# test <- head(obj_EBV_modality_integrate_cluster_TCR_cluster)

for (i in 1:length(df)) {
  for (j in 1:length(df)) {
    df[i,j] <- as.data.frame(table(obj_EBV_modality_integrate_cluster_TCR_cluster[,i] + obj_EBV_modality_integrate_cluster_TCR_cluster[,j] == 2)[2])[,1]
  }
}

df$cluster <- rownames(df)
hm <- melt(df)
colnames(hm) <- c("Cluster1","Cluster2",'Clones')

hm$Cluster1 <- factor(hm$Cluster1,levels = paste("C",c(0:12),sep = ""))
hm$Cluster2 <- factor(hm$Cluster2,levels = paste("C",c(0:12),sep = ""))

library(ggplot2)
p <- ggplot(hm, aes(x=Cluster1, y=Cluster2, fill=Clones)) +
  geom_tile(color="black") + theme_bw() + coord_equal() +
  scale_fill_distiller(palette="RdYlBu", direction=1) +
  scale_fill_gradientn(limits = c(0,100),colours=c(ArchRPalettes$solarExtra)) +
  # guides(fill=F) + # removing legend for `fill`
  labs(title = "Value distribution") + # using a title instead
  geom_text(aes(label=Clones), color="black") # printing values

dir.create(paste(savedir,"heatmap",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"heatmap/TCR_number.pdf",sep=""))
p
dev.off()

## Celltypes
n = unique(obj_EBV_modality_integrate_cluster@meta.data$Celltypes)
df <- data.frame(matrix(ncol = length(n), nrow = length(n)))
colnames(df) <- n
rownames(df) <- n

obj_EBV_modality_integrate_cluster_TCR_cluster <- table(obj_EBV_modality_integrate_cluster@meta.data$TCR, obj_EBV_modality_integrate_cluster@meta.data$Celltypes)
obj_EBV_modality_integrate_cluster_TCR_cluster <- obj_EBV_modality_integrate_cluster_TCR_cluster[grep("^_$",rownames(obj_EBV_modality_integrate_cluster_TCR_cluster),invert=TRUE),]
head(obj_EBV_modality_integrate_cluster_TCR_cluster)
obj_EBV_modality_integrate_cluster_TCR_cluster[obj_EBV_modality_integrate_cluster_TCR_cluster > 0] <- 1
# test <- head(obj_EBV_modality_integrate_cluster_TCR_cluster)

for (i in 1:length(df)) {
  for (j in 1:length(df)) {
    df[i,j] <- as.data.frame(table(obj_EBV_modality_integrate_cluster_TCR_cluster[,i] + obj_EBV_modality_integrate_cluster_TCR_cluster[,j] == 2)[2])[,1]
  }
}

# colnames(df) <- c("CM_EM","EM","Naive","TEMRA","Cycling_2")
# rownames(df)  <- c("CM_EM","EM","Naive","TEMRA","Cycling_2")

df$cluster <- rownames(df)
hm <- melt(df)
colnames(hm) <- c("Cluster1","Cluster2",'Clones')
hm$Cluster1 <- factor(hm$Cluster1,levels = c("Naive","CM/EM","EM","TEMRA","unknown"))
hm$Cluster2 <- factor(hm$Cluster2,levels = c("Naive","CM/EM","EM","TEMRA","unknown"))

library(ggplot2)
p <- ggplot(hm, aes(x=Cluster1, y=Cluster2, fill=Clones)) +
  geom_tile(color="black") + theme_bw() + coord_equal() +
  scale_fill_distiller(palette="RdYlBu", direction=1) +
  # guides(fill=F) + # removing legend for `fill`
  scale_fill_gradientn(limits = c(0,230),colours=c(ArchRPalettes$solarExtra)) +
  labs(title = "TRB Value distribution") + # using a title instead
  geom_text(aes(label=Clones), color="black") # printing values


dir.create(paste(savedir,"heatmap",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"heatmap/celltypes_number.pdf",sep=""))
p
dev.off()


### Using Bhattacharya cofficient 
# from this paper https://www.sciencedirect.com/science/article/pii/S0092867421013337?via%3Dihub#bib28
library(stringr)

n = levels(obj_EBV_modality_integrate_cluster@meta.data$seurat_clusters)
df <- data.frame(matrix(ncol = length(n), nrow = length(n)))
colnames(df) <- paste("C",n,sep = "")
rownames(df) <- paste("C",n,sep = "")

obj_EBV_modality_integrate_cluster_TCR_cluster <- table(obj_EBV_modality_integrate_cluster@meta.data$TCR, obj_EBV_modality_integrate_cluster@meta.data$seurat_clusters)
obj_EBV_modality_integrate_cluster_TCR_cluster <- obj_EBV_modality_integrate_cluster_TCR_cluster[grep("^_$",rownames(obj_EBV_modality_integrate_cluster_TCR_cluster),invert=TRUE),]
head(obj_EBV_modality_integrate_cluster_TCR_cluster)
# obj_EBV_modality_integrate_cluster_TCR_cluster[obj_EBV_modality_integrate_cluster_TCR_cluster > 0] <- 1
# test <- head(obj_EBV_modality_integrate_cluster_TCR_cluster,20) %>% tail(10) ## this is more interesting
column_sum <- colSums(obj_EBV_modality_integrate_cluster_TCR_cluster) %>% as.vector()
vec <- c()

for (i in 1:length(df)) {
  for (j in 1:length(df)) {
    req <- obj_EBV_modality_integrate_cluster_TCR_cluster[obj_EBV_modality_integrate_cluster_TCR_cluster[,i] + obj_EBV_modality_integrate_cluster_TCR_cluster[,j] > 1,]
    if(length(req) > 0){
      if (length(req) > 13) {
        req2 <- req[,c(i,j)]
        vec <- c()
        for (k in 1:nrow(req2)) {
          vec[k] <- sqrt((req2[k,1]/column_sum[i])*(req2[k,2]/column_sum[j]))
        }
      }
      else if (length(req) == 13) {
        req_df <- as.data.frame(req) 
        req2 <- req_df[c(i,j),]
        vec <- c()
        vec <- sqrt((req2[1]/column_sum[i])*(req2[2]/column_sum[j]))
      }
      df[i,j] <- sum(vec)
    }
  }
}

# colnames(df) <- c("TCF1hi","Tfh_like","Cycling_1","T_Cyto","Cycling_2")
# rownames(df)  <- c("TCF1hi","Tfh_like","Cycling_1","T_Cyto","Cycling_2")

df$cluster <- rownames(df)
hm <- melt(df)
colnames(hm) <- c("Cluster1","Cluster2",'TRSS')
hm$TRSS <- round(hm$TRSS, digits = 2)
hm$Cluster1 <- factor(hm$Cluster1,levels = paste("C",c(0:12),sep = ""))
hm$Cluster2 <- factor(hm$Cluster2,levels = paste("C",c(0:12),sep = ""))

library(ggplot2)
p <- ggplot(hm, aes(x=Cluster1, y=Cluster2, fill=TRSS)) +
  geom_tile(color="black") + theme_bw() + coord_equal() +
  scale_fill_distiller(palette="RdYlBu", direction=1) +
  # guides(fill=F) + # removing legend for `fill`
  scale_fill_gradientn(limits = c(0,0.60),colours=c(ArchRPalettes$solarExtra)) +
  labs(title = "TRB Value distribution") + # using a title instead
  geom_text(aes(label=TRSS), color="black") # printing values


dir.create(paste(savedir,"heatmap",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"heatmap/celltypes_TCR_only_proportion_TRSS.pdf",sep=""))
p
dev.off()

write.table(df,paste(savedir,"heatmap/celltypes_TRB_only_proportion_TRSS_table.txt",sep=""), quote = F, 
            row.names = T, col.names = T, sep = "\t")


### Shannon Diversity ####
library(Seurat)
library(ggplot2)
library(dplyr)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/"
EBV_obj <- readRDS(paste(savedir,"saveRDS_obj/obj_EBV_modality_integrate_cluster_0.8_2.RDS",sep = ""))
metadata <- EBV_obj@meta.data
metadata_req <- metadata[grep("^_$",metadata$TCR,invert=TRUE),]

donor_antigen_count <- table(paste(metadata_req$donor, metadata_req$specificity.final,sep = "_")) %>% as.data.frame()
donor_antigen_count_req <- donor_antigen_count[donor_antigen_count[,"Freq"] > 10,]
metadata_req$donor_antigen <- paste(metadata_req$donor, metadata_req$specificity.final,sep = "_")
metadata_req_sample <- metadata_req[grep(paste("^",donor_antigen_count_req$Var1,"$",sep = "",collapse = "|"), metadata_req$donor_antigen),]

cellnames=rownames(metadata_req_sample)
EBV_obj$donor_antigen <- paste(EBV_obj$donor, EBV_obj$specificity.final,sep = "#")
EBV_obj_subset <- subset(EBV_obj, cells = cellnames)

DefaultAssay(EBV_obj_subset) <- "RNA"

# Clonal Sharing
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/VZV_clonal_sharing.R")
clus_length = length(unique(EBV_obj_subset@meta.data$seurat_clusters))
donor_antigen <- unique(EBV_obj_subset@meta.data$donor_antigen)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/TCR/"
splitting = "#"

for (i in 1:length(donor_antigen)) {
  TCR_AJ_VZV(object = EBV_obj_subset, savedir = paste(savedir2,sep = ""),
             clus_col = "seurat_clusters", TCR_col = "TCR",
             group_col = "donor_antigen", group_val = donor_antigen[i],
             split_col = "donor_antigen",
             column_name = c("ID","Ag"), 
             total_clusters = clus_length, splitting = splitting)
}

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Shannon_Diversity.R")
rm(df2)
sample_id <- EBV_obj_subset@meta.data$donor_antigen %>% unique()
cluster_num <- paste("C",unique(EBV_obj_subset@meta.data$seurat_clusters)[order(unique(EBV_obj_subset@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"diversity",sep = ""),showWarnings = FALSE)
rm(sha_ent)
splitting = "#"

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Shannon_Diversity.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Shannon_AJ(object = EBV_obj_subset, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TCR",
                             group_col = "donor_antigen", group_val = sample_id[i], split_col = "donor_antigen",
                             column_name = c("ID","Ag"), total_clusters = length(cluster_num), splitting = splitting)
  })
}

dir.create(paste(savedir,"diversity/",sep = ""), showWarnings = FALSE)
write.table(df2, paste(savedir,"diversity/Shannon_sample_diversity_all_cluster.txt",sep = ""), quote = F,
            row.names = T, col.names = T, sep = "\t")
# df2 <- read.table(paste(savedir,"diversity/same_size/Shannon_sample_diversity_all_cluster.txt",sep = ""), header = TRUE)
df2$sample <- row.names(df2)
df2_melted <- melt(df2)
colnames(df2_melted) <- c("Sample","Cluster","Shannon_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "#")[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name

for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Shannon_diversity)) + 
  geom_boxplot(aes(fill = Ag)) + geom_point(position=position_dodge(width=0.75),aes(group=Ag)) +
  ggtitle(paste("Shannon Diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"diversity/Ag_donor_shannon_diversity_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### Disease group
df2_melted$group <- gsub("_.*","",df2_melted$ID)

p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Shannon_diversity)) + 
  geom_boxplot(aes(fill = group)) + 
  geom_point(position=position_dodge(width=0.75),aes(group=group)) +
  ggtitle(paste("Shannon diversity Disease group")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"diversity/Disease_group_shannon_diversity_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### Simpson Diversity index #####
rm(df2)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/TCR/"
sample_id <- EBV_obj_subset@meta.data$donor_antigen %>% unique()
cluster_num <- paste("C",unique(EBV_obj_subset@meta.data$seurat_clusters)[order(unique(EBV_obj_subset@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"diversity/",sep = ""),showWarnings = FALSE)
rm(sha_ent)
splitting = "#"

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Simpson_Diversity.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Simpson_AJ(object = EBV_obj_subset, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TCR",
                             group_col = "donor_antigen", group_val = sample_id[i], split_col = "donor_antigen",
                             column_name = c("ID","Ag"), total_clusters = length(cluster_num), splitting = splitting)
  })
}

write.table(df2, paste(savedir,"diversity/Simpson_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df2 <- read.delim(paste(savedir,"diversity/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df2$sample <- row.names(df2)
df2_melted <- melt(df2, id = "sample")
colnames(df2_melted) <- c("Sample","Cluster","Simpson")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "#")[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name

if("condition" %in% colnames(df2_melted_split) & "vaccine" %in% colnames(df2_melted_split)){
  df2_melted_split[grep("DMSO",df2_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

#### Gini Simpson Diversity index ####
df2 <- read.delim(paste(savedir,"diversity/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1-df2

write.table(df3, paste(savedir,"diversity/Gini_Simpson_diversity_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Inverse Simpson index #####
df2 <- read.delim(paste(savedir,"diversity/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1/df2

write.table(df3, paste(savedir,"diversity/inverse_Simpson_diversity_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df3$sample <- row.names(df3)
df3_melted <- melt(df3, id="sample")
colnames(df3_melted) <- c("Sample","Cluster","inverse_Simpson_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "#")[[1]]
  data.frame(t(split_elements))
}

df3_melted_split <- do.call(rbind, lapply(df3_melted[,"Sample"], split_to_dataframe))
colnames(df3_melted_split) <- column_name

if("condition" %in% colnames(df3_melted_split) & "vaccine" %in% colnames(df3_melted_split)){
  df3_melted_split[grep("DMSO",df3_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df3_melted[,column_name[i]] <- df3_melted_split[,column_name[i]]
}

### Age 
p <- ggplot(na.omit(df3_melted), aes(x=Cluster,y=inverse_Simpson_diversity)) + 
  geom_boxplot(aes(fill = Ag)) + 
  geom_point(position=position_dodge(width=0.75),aes(group=Ag)) +
  ggtitle(paste("Inverse Simpson Age")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"diversity/Ag_inverse_simpson_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

df3_melted$group <- gsub("_.*","", df3_melted$ID)
p <- ggplot(na.omit(df3_melted), aes(x=Cluster,y=inverse_Simpson_diversity)) + 
  geom_boxplot(aes(fill = group)) + 
  geom_point(position=position_dodge(width=0.75),aes(group=group)) +
  ggtitle(paste("Inverse Simpson Disease group")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"diversity/disease_group_inverse_Simpson_diversity_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

#### Gini Clonality Index #####
rm(df2)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/TCR/"
sample_id <- EBV_obj_subset@meta.data$donor_antigen %>% unique()
cluster_num <- paste("C",unique(EBV_obj_subset@meta.data$seurat_clusters)[order(unique(EBV_obj_subset@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"expansion",sep = ""),showWarnings = FALSE)
rm(sha_ent)
splitting = "#"

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Gini_coefficient_index.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Gini_coef_AJ(object = EBV_obj_subset, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TCR",
                               group_col = "donor_antigen", group_val = sample_id[i], split_col = "donor_antigen",
                               column_name = c("ID","Ag"), total_clusters = length(cluster_num),splitting = splitting)
  })
}

write.table(df2, paste(savedir,"expansion/Gini_coef_clonality_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df2 <- read.delim(paste(savedir,"expansion/Gini_coef_clonality_sample_all_cluster.txt",sep = ""), header = TRUE)
df2$sample <- row.names(df2)
df2_melted <- melt(df2, id = "sample")
colnames(df2_melted) <- c("Sample","Cluster","Gini_coef")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, splitting)[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name

if("condition" %in% colnames(df2_melted_split) & "vaccine" %in% colnames(df2_melted_split)){
  df2_melted_split[grep("DMSO",df2_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}
### Age 
p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Gini_coef)) + 
  geom_boxplot(aes(fill = Ag)) + 
  geom_point(position=position_dodge(width=0.75),aes(group=Ag)) +
  ggtitle(paste("Gini coefficient Clonality")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"expansion/Ag_Gini_coef_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

## group
df2_melted$group <- gsub("_.*","",df2_melted$ID)
p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Gini_coef)) + 
  geom_boxplot(aes(fill = group)) + 
  geom_point(position=position_dodge(width=0.75),aes(group=group)) +
  ggtitle(paste("Gini coefficient Clonality")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"expansion/Disease_group_Gini_coef_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

all_EBV_clusters <- as.matrix(table(EBV_obj_subset@meta.data$donor_antigen,EBV_obj_subset@meta.data$seurat_clusters))
write.table(all_EBV_clusters, paste(savedir,"all_EBV_clusters.txt",sep = ""), sep = "\t", row.names = T,col.names = T, quote = F)

### BC Celltypes ####
library(stringr)

n = unique(obj_EBV_modality_integrate_cluster@meta.data$Celltypes)
df <- data.frame(matrix(ncol = length(n), nrow = length(n)))
colnames(df) <- n
rownames(df) <- n

obj_EBV_modality_integrate_cluster_TCR_cluster <- table(obj_EBV_modality_integrate_cluster@meta.data$TCR, obj_EBV_modality_integrate_cluster@meta.data$Celltypes)
obj_EBV_modality_integrate_cluster_TCR_cluster <- obj_EBV_modality_integrate_cluster_TCR_cluster[grep("^_$",rownames(obj_EBV_modality_integrate_cluster_TCR_cluster),invert=TRUE),]
head(obj_EBV_modality_integrate_cluster_TCR_cluster)
# obj_EBV_modality_integrate_cluster_TCR_cluster[obj_EBV_modality_integrate_cluster_TCR_cluster > 0] <- 1
# test <- head(obj_EBV_modality_integrate_cluster_TCR_cluster,20) %>% tail(10) ## this is more interesting
column_sum <- colSums(obj_EBV_modality_integrate_cluster_TCR_cluster) %>% as.vector()
vec <- c()

for (i in 1:length(df)) {
  for (j in 1:length(df)) {
    req <- obj_EBV_modality_integrate_cluster_TCR_cluster[obj_EBV_modality_integrate_cluster_TCR_cluster[,i] + obj_EBV_modality_integrate_cluster_TCR_cluster[,j] > 1,]
    if(length(req) > 0){
      if (length(req) > 5) {
        req2 <- req[,c(i,j)]
        vec <- c()
        for (k in 1:nrow(req2)) {
          vec[k] <- sqrt((req2[k,1]/column_sum[i])*(req2[k,2]/column_sum[j]))
        }
      }
      else if (length(req) == 5) {
        req_df <- as.data.frame(req)
        req2 <- req_df[c(i,j),]
        vec <- c()
        vec <- sqrt((req2[1]/column_sum[i])*(req2[2]/column_sum[j]))
      }
      df[i,j] <- sum(vec)
    }
  }
}

# colnames(df) <- c("TCF1hi","Tfh_like","Cycling_1","T_Cyto","Cycling_2")
# rownames(df)  <- c("TCF1hi","Tfh_like","Cycling_1","T_Cyto","Cycling_2")

df$cluster <- rownames(df)
hm <- melt(df)
colnames(hm) <- c("Cluster1","Cluster2",'TRSS')
hm$TRSS <- round(hm$TRSS, digits = 2)
hm$Cluster1 <- factor(hm$Cluster1,levels = c("Naive","CM/EM","EM","TEMRA","unknown"))
hm$Cluster2 <- factor(hm$Cluster2,levels = c("Naive","CM/EM","EM","TEMRA","unknown"))

library(ggplot2)
p <- ggplot(hm, aes(x=Cluster1, y=Cluster2, fill=TRSS)) +
  geom_tile(color="black") + theme_bw() + coord_equal() +
  scale_fill_distiller(palette="RdYlBu", direction=1) +
  # guides(fill=F) + # removing legend for `fill`
  scale_fill_gradientn(limits = c(0,0.50),colours=c(ArchRPalettes$solarExtra)) +
  labs(title = "TRB Value distribution") + # using a title instead
  geom_text(aes(label=TRSS), color="black") # printing values


dir.create(paste(savedir,"heatmap",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"heatmap/celltypes_TCR_only_proportion_TRSS.pdf",sep=""))
p
dev.off()

write.table(df,paste(savedir,"heatmap/celltypes_TCR_only_proportion_TRSS_table.txt",sep=""), quote = F,
            row.names = T, col.names = T, sep = "\t")

#### Automation for each Antigen #####
Antigen <- unique(obj_EBV_modality_integrate_cluster@meta.data$specificity.final)
Antigen <- grep("undetermined",Antigen,invert = TRUE,value = TRUE)
obj <- obj_EBV_modality_integrate_cluster
Idents(obj) <- obj_EBV_modality_integrate_cluster@meta.data$specificity.final

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/EBV_TCR_TRSS.R")
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/TCR/"
for (i in 1:length(Antigen)) {
  obj_subset <- subset(obj, idents = Antigen[i])
  TCR_AJ_EBV(object = obj_subset, savedir = paste(savedir,Antigen[i],"/",sep = ""),objname = Antigen[i])
}

### For All
TCR_AJ_EBV(object = obj, savedir = paste(savedir,"All/",sep = ""),objname = "All")

#### Pseudobulk Analysis ####
library(Seurat)
EBV_obj <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/Analysis/saveRDS_obj/obj_modality_integrate_cluster_1.RDS")
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/"
dir.create(savedir, showWarnings = FALSE)

EBV_lytic <- grep("BLMF1|brlf1|bzlf1",EBV_obj@meta.data$specificity.final, ignore.case = T,value=TRUE) %>% grep("undetermined",.,invert=TRUE,value=TRUE)  %>% unique()
EBV_latent <- grep("BLMF1|brlf1|bzlf1",EBV_obj@meta.data$specificity.final, ignore.case = T,value=TRUE,invert=TRUE) %>% grep("undetermined",.,invert=TRUE,value=TRUE) %>% unique()
EBV_obj@meta.data$specificity.final -> EBV_obj@meta.data$antigen
EBV_obj@meta.data$antigen <- gsub(paste("^",EBV_lytic,"$", collapse = "|",sep = ""),"lytic",EBV_obj@meta.data$antigen)
EBV_obj@meta.data$antigen <- gsub(paste("^",EBV_latent,"$", collapse = "|",sep = ""),"latent",EBV_obj@meta.data$antigen)
EBV_obj@meta.data$antigen <- gsub(paste('.*.undetermined', collapse = "|",sep = ""),"undetermined",EBV_obj@meta.data$antigen)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/EBV_pseudobulk_PCA_within.R")
EBV_obj_subset_dds <- pseudobulk_within_cluster_AJ(obj = EBV_obj, savedir = savedir,group1 = "lytic", group2 = "latent", 
                                               grouping_by = "antigen", cluster = "all",cell_freq = 10,remove_samples = "undetermined",
                                               gene_min_counts =  1, batch_col = "chip",sample_col = "donor",cluster_group = "seurat_clusters")


savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/pseudobulk/clus_all_removed_undetermined_lytic_vs_latent/"
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/RNAseq/pipeline_function/RNAseq_limma_EdgeR.R")
design0 <- model.matrix(~ 0 + antigen + Run, data = colData(EBV_obj_subset_dds))
colnames(design0) <- c("latent","lytic","Run02","Run03","Run04")
cm <- makeContrasts(lytic_VS_latent = lytic-latent,levels = design0)
desl_clus <- LimmaEdgeR_differential(dds = EBV_obj_subset_dds,
                                     design0 = design0,
                                     cm = cm, 
                                     savedir = savedir,
                                     logfc = 0.5,
                                     p_value_adj = 0.05)


#### Clusters Differential ####
clus <- levels(EBV_obj@meta.data$seurat_cluster)
for (i in 1:length(clus)) {
  savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/"
  dir.create(savedir, showWarnings = FALSE)
  source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/EBV_pseudobulk_PCA_within.R")
  EBV_obj_subset_dds <- pseudobulk_within_cluster_AJ(obj = EBV_obj, savedir = savedir,group1 = "lytic", group2 = "latent", 
                                                     grouping_by = "antigen", cluster = clus[i],cell_freq = 5,remove_samples = "undetermined",
                                                     gene_min_counts =  1, batch_col = "chip",sample_col = "donor",cluster_group = "seurat_clusters")
  
  
  savedir <- paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/pseudobulk/clus_",clus[i],"_removed_undetermined_lytic_vs_latent/",sep = "")
  source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/RNAseq/pipeline_function/RNAseq_limma_EdgeR.R")
  design0 <- model.matrix(~ 0 + antigen + Run, data = colData(EBV_obj_subset_dds))
  colnames(design0) <- c("latent","lytic","Run02","Run03","Run04")
  cm <- makeContrasts(lytic_VS_latent = lytic-latent,levels = design0)
  desl_clus <- LimmaEdgeR_differential(dds = EBV_obj_subset_dds,
                                       design0 = design0,
                                       cm = cm, 
                                       savedir = savedir,
                                       logfc = 0.5,
                                       p_value_adj = 0.05)
}

#### Celltype Differential ####
EBV_obj@meta.data$Celltypes <- gsub("CM/EM","CM_EM",EBV_obj@meta.data$Celltypes)
clus <- unique(EBV_obj@meta.data$Celltypes)
for (i in 1:length(clus)) {
  savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/"
  dir.create(savedir, showWarnings = FALSE)
  source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/EBV_pseudobulk_PCA_within.R")
  EBV_obj_subset_dds <- pseudobulk_within_cluster_AJ(obj = EBV_obj, savedir = savedir,group1 = "lytic", group2 = "latent", 
                                                     grouping_by = "antigen", cluster = clus[i],cell_freq = 5,remove_samples = "undetermined",
                                                     gene_min_counts =  1, batch_col = "chip",sample_col = "donor",cluster_group = "Celltypes")
  
  
  savedir <- paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/EBV_healthy_T1D/pseudobulk/clus_",clus[i],"_removed_undetermined_lytic_vs_latent/",sep = "")
  source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/RNAseq/pipeline_function/RNAseq_limma_EdgeR.R")
  design0 <- model.matrix(~ 0 + antigen + Run, data = colData(EBV_obj_subset_dds))
  colnames(design0) <- c("latent","lytic","Run02","Run03","Run04")
  cm <- makeContrasts(lytic_VS_latent = lytic-latent,levels = design0)
  desl_clus <- LimmaEdgeR_differential(dds = EBV_obj_subset_dds,
                                       design0 = design0,
                                       cm = cm, 
                                       savedir = savedir,
                                       logfc = 0.5,
                                       p_value_adj = 0.05)
}



#### Lutian Analysis ####
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/Lutian_analysis/"
dir.create(savedir, showWarnings = FALSE)
Ag_mem_Y <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/Lutian_analysis/Ag_freq", header = TRUE, sep = "\t", row.names = 1)

metadata <- data.frame(samples = rownames(Ag_mem_Y))
rownames(metadata) <- metadata$samples
metadata$Ag <- gsub(".*._","",metadata$samples)
metadata$samplename <- gsub("_.*.","",metadata$samples)

### Running EdgeR
set.seed(123)
library(stringr)
library(DESeq2)
library(edgeR)
#Note that normalization in edgeR is model-based, and the original read counts are not themselves
#transformed. This means that users should not transform the read counts in any way before
#inputing them to edgeR.For example, users should not enter RPKM or FPKM val- ues to edgeR in place of read counts.
#edgeR can work with expected counts as output by RSEM, but raw counts are still preferred.
library(RColorBrewer)
# library(Glimma)
library(ggplot2)
library(ggrepel)
# library(session)
# library(DEFormats)
# library(EnhancedVolcano)
library(limma)
library(stringr)
all(metadata$samples == rownames(Ag_mem_Y)) # if TRUE move forward
# Th17_metadata$samplename <- gsub("gE_|DMSO_|_DMSO|_gE|B22015","",Th17_metadata$samples)
## Run will be treated as a covariate in the regression model,
Ag_mem_Y_t <- t(Ag_mem_Y)
pc <- prcomp(Ag_mem_Y_t)
attributes(pc)
pc$center
summary(pc)

pc$rotation

rownames(metadata) <- metadata$samples
pc_with_metadata <- merge(pc$rotation,metadata,by = 'row.names', all = TRUE)

write.table(pc_with_metadata, file = paste(savedir,"PCA_coordinate_matched_samples_with_metadata.txt",sep = ""),
            quote = F, col.names = T, row.names = F, sep = "\t")

# pc_with_metadata <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/all_run_analysis_CR7/PCA/PCA/absolute_PCA_coordinate_with_metadata.txt",
#                                header = TRUE)
# pc_with_metadata$vaccine <- gsub("S$","Shingrix",pc_with_metadata$vaccine)
# pc_with_metadata$vaccine <- gsub("Z$","Zostavax",pc_with_metadata$vaccine)

pc_summary <- summary(pc)
pc_summary_imp <- as.data.frame(pc_summary$importance)
pc_with_metadata
rm(plot_list)
plot_list <- list()
plot_list[[1]] <-ggplot(pc_with_metadata, aes(PC1, PC2, label=samples, color=Ag, shape=samplename)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(pc_summary_imp[2,1])*100,"% variance")) +
  ylab(paste0("PC2: ",(pc_summary_imp[2,2])*100,"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

plot_list[[2]] <-ggplot(pc_with_metadata, aes(PC3, PC4, label=samples, color=Ag, shape=samplename)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(pc_summary_imp[2,3])*100,"% variance")) +
  ylab(paste0("PC2: ",(pc_summary_imp[2,4])*100,"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

plot_list[[3]] <-ggplot(pc_with_metadata, aes(PC5, PC6, label=samples, color=Ag, shape=samplename)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(pc_summary_imp[2,5])*100,"% variance")) +
  ylab(paste0("PC2: ",(pc_summary_imp[2,6])*100,"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))


pdf(paste(savedir,"PCA_norm_samples.pdf",sep = ""),width = 9, height = 8)
plot_list
dev.off()

### Antigen Lutian Analysis
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/Lutian_analysis/"
dir.create(savedir, showWarnings = FALSE)
Ag_mem_Y <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/Lutian_analysis/Ag_freq", header = TRUE, sep = "\t", row.names = 1)
metadata <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/Nat_Immunology_Jiang_EBV/Lutian_analysis/Ag_freq_metadata", header = TRUE, sep = "\t")
samplename <- gsub("_EBV.*.|_EBNA.*.|_BZLF1.*.","",rownames(Ag_mem_Y))

Ag_mem_Y$samplename <- metadata[match(samplename, metadata$Samplename),1]
Ag_mem_Y$Sex <- metadata[match(samplename, metadata$Samplename),"Sex"]
Ag_mem_Y$Age_num <- metadata[match(samplename, metadata$Samplename),"Age_num"]
Ag_mem_Y$Ags <- gsub(".*._","",rownames(Ag_mem_Y))

metadata_2 <- Ag_mem_Y[,c("samplename","Sex","Age_num","Ags")]
rownames(metadata_2) <- paste(metadata_2$samplename, metadata_2$Sex, metadata_2$Age_num, metadata_2$Ags, sep = "~") 

### combining the samplename for the metadata
### Only the Samplename
Ag_mem_Y <- Ag_mem_Y[,1:13]

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "~")[[1]]
  data.frame(t(split_elements))
}

# Convert split elements to dataframe
metadata_3 <- do.call(rbind, lapply(rownames(metadata_2), split_to_dataframe))
colnames(metadata_3) <- c("samplename","Sex","Age_num","Ags")
metadata_3$samples <- metadata_2$samples

Ag <- unique(metadata_3$Ags)

df <- data.frame(matrix(nrow = length(Ag), ncol = length(Ag)))
rownames(df) <- Ag
colnames(df) <- Ag

for (i in 1:(length(Ag))) {
  for (j in 1:length(Ag)) {
    print(paste("Performing Lutian analysis for ",Ag[i]," and ",Ag[j]))
    data <- Ag_mem_Y[grep(paste(Ag[i],"|",Ag[j],sep = ""),rownames(Ag_mem_Y)),]
    metadata_4 <- metadata_3[grep(paste(Ag[i],"|",Ag[j],sep = ""),metadata_3$Ags),]
    
    if (Ag[i] == Ag[j]) {
      print(paste(Ag[i]," is same as ",Ag[j]))
    } else if (Ag[i] != Ag[j]){
      Sex <- metadata_4$Sex
      Ag_var <- metadata_4$Ags
      Ag_var <- as.numeric(factor(Ag_var))
      female=1*(Sex=="Female")
      
      print(dim(data))
      
      datax <- data[,(colSums(data)) > 0]
      clus <- dim(datax)[2]
      
      z=rep(NA,clus)
      for(k in 1:clus){
        fit=summary(lm(data[,k]~Ag_var+female))$coef
        z[k]=fit[2,3]
      }
      z.obs=sum(z^2)
      
      B=10000
      z.null=rep(NA, B)
      n1=sum(female==1)
      n0=sum(female==0)
      for(b in 1:B){
        Ag_var0=Ag_var
        Ag_var0[female==1]=sample(Ag_var[female==1], n1, F)
        Ag_var0[female==0]=sample(Ag_var[female==0], n0, F)
        z=rep(NA,clus)
        for(k in 1:clus)
        {
          fit=summary(lm(datax[,k]~Ag_var0+female))$coef
          z[k]=fit[2,3]
        }
        z.null[b]=sum(z^2)
      }
      # print(paste("Removed this vaccine",vaccine[i],"for this time",time[j],sep = " "))
      print(paste("pvalue for ",Ag_var[i],"and",Ag_var[j],sep = " "))
      print(mean(z.null>=z.obs))
      df[i,j] <- mean(z.null>=z.obs)
    }
  }
}

write.table(df, paste(savedir,"Lutian_Y_only_Ag.txt",sep = ""), sep = "\t", quote = F, row.names = T, col.names = T)

savedir











