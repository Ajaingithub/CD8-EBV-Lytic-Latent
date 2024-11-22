#### Nan Ping scRNA ####
# We have downloaded the data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136184
# qlogin -q lg-mem -l h_vmem=120G -pe threaded 10
library(Seurat)
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/"
p1p4 <- readRDS(paste(savedir,"GSE136184_seu_obj.Rds",sep = ""))
metadata <- read.table(paste(savedir,"individual_cell_counts_metadata.txt",sep = ""), header = TRUE)
metadata$study <- metadata$donor_ID
metadata[grep("sc",metadata$study),"study"] <- "cross"
metadata[grep("L",metadata$study),"study"] <- "long"
p1p4@meta.data$Code_visit <- paste(p1p4@meta.data$Code,p1p4@meta.data$visit,sep="_")

library(stringi)
p1p4@meta.data$Donor_id <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                                  pattern=metadata$Code_visit,
                                                  replacement=metadata$donor_ID,
                                                  vectorize=FALSE)

p1p4@meta.data$Age <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                             pattern=metadata$Code_visit,
                                             replacement=metadata$age,
                                             vectorize=FALSE)

p1p4@meta.data$Sex <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                             pattern=metadata$Code_visit,
                                             replacement=metadata$sex,
                                             vectorize=FALSE)

p1p4@meta.data$study <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                               pattern=metadata$Code_visit,
                                               replacement=metadata$study,
                                               vectorize=FALSE)

p1p4@meta.data$orig.ident <- paste(p1p4@meta.data$Code_visit,p1p4@meta.data$Age,p1p4@meta.data$Sex,p1p4@meta.data$Donor_id,sep="_")

p1p4[["percent.mt"]] <- PercentageFeatureSet(p1p4, pattern = "^MT-")
p1p4$mitoRatio <- p1p4@meta.data$percent.mt / 100

# Changing the Idents
all(names(Idents(p1p4))==rownames(p1p4@meta.data))
p1p4 <- SetIdent(p1p4, value = p1p4@meta.data$orig.ident)

pdf(paste(savedir,"p1p4_QC_mt_percent.pdf",sep = ""), width = 20, height = 7)
VlnPlot(p1p4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.000005, ncol = 3)
dev.off()

pdf(paste(savedir,"p1p4_QC_mt_percent_study.pdf",sep = ""), width = 20, height = 7)
VlnPlot(p1p4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "study", pt.size = 0.000005, ncol = 3)
dev.off()

### saveRDS obj
# saveRDS(p1p4,paste(savedir,"p1p4.RDS",sep = ""))
p1p4 <- readRDS(paste(savedir,"p1p4.RDS",sep = ""))

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_sctransform_V2.R")
objname = "p1p4"
Assay = "RNA"
process = "sctransform"
ngenes=1000
dims=30
regress="mitoRatio"
split_by = "study"
obj=p1p4
reference = 1
p1p4_sctransformed <- sctransform_V2_integration(obj = p1p4, saveDir = savedir, ngenes = 1000,
                                                 regress = c("mitoRatio"),
                                                 dims = 30,
                                                 Assay = Assay, process = process, objname = objname,
                                                 split_by = "study",
                                                 reference = 1,
                                                 sample_tree = NULL)

# ### Following Nan-ping code
p1p4.list <- SplitObject(p1p4, split.by = "study")

for (i in 1:length(p1p4.list)) {
  # p1p4.list[[i]] <- NormalizeData(p1p4.list[[i]], verbose = T)
  p1p4.list[[i]] <- FindVariableFeatures(p1p4.list[[i]], verbose = T, nfeatures = 1000)
}

reference_dataset <- which(names(p1p4.list) == "cross")
pbmc.anchors <- FindIntegrationAnchors(object.list = p1p4.list, dims = 1:50, reference = reference_dataset)
# p1p4.list <- PrepSCTIntegration(object.list = p1p4.list, anchor.features = pbmc.anchors@anchor.features)
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, dims = 1:50)
pbmc.integrated <- ScaleData(pbmc.integrated)
pbmc.integrated <- RunPCA(object = pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunUMAP(pbmc.integrated, reduction = "pca", dims = 1:30)

## Plot Data
pdf(paste(savedir,"dimplot_PCA.pdf",sep = ""),width = 9, height = 6)
DimPlot(pbmc.integrated, reduction = 'pca', group.by = 'study', split.by = "study")
dev.off()

dims <- c(1:10)
perplexities <- c(30)
seed <- c(7456)
iters <- c(3000)
pbmc.integrated <- RunTSNE(method = 'FIt-SNE', object = pbmc.integrated, dims = dims, perplexity = perplexities, seed.use = seed, max_iter = iters)

pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:30)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 0.8)

pdf(paste(savedir,"dimplot_tnse.pdf",sep = ""),width = 9, height = 6)
DimPlot(pbmc.integrated, reduction = 'tsne', group.by = 'study', split.by = "study")
dev.off()

pdf(paste(savedir,"dimplot_umap.pdf",sep = ""),width = 9, height = 6)
DimPlot(pbmc.integrated, reduction = 'umap', group.by = 'study', split.by = "study")
dev.off()

pdf(paste(savedir,"dimplot.pdf",sep = ""),width = 6, height = 5.5)
DimPlot(pbmc.integrated, reduction = 'tsne')
dev.off()

rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated) = "RNA"
genes <- c("CCR7","IKZF2")
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated,genes[i], reduction = "umap", x="UMAP_1", y="UMAP_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

p=featureplot_front(pbmc.integrated,"IKZF2", reduction = "umap", x="UMAP_1", y="UMAP_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits=c(0,0.3))

dir.create(paste(savedir,"featureplot",sep = ""))
pdf(paste(savedir,"featureplot/IKZF2_MAGIC_umap.pdf",sep = ""))
p
dev.off()

require(gridExtra)
pdf(paste(savedir,"featureplot/genes.pdf",sep = ""), width = 9, height = 5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()

rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated) = "RNA"
genes <- c("CCR7","IKZF2")
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated,genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

require(gridExtra)
pdf(paste(savedir,"featureplot/genes_tSNE.pdf",sep = ""), width = 9, height = 5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()

DefaultAssay(pbmc.integrated) <- "RNA"
gene_zero <- as.vector(which(rowSums(pbmc.integrated@assays$RNA@data)==0))
pbmc.integrated_2 <- subset(pbmc.integrated, features = rownames(pbmc.integrated[-gene_zero,]))
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/python/python-3.9.2/bin/python3.9")
library(Rmagic)
py_discover_config("magic") # to check
pbmc.integrated_impute <- magic(pbmc.integrated_2, npca=30) ## imputing the vasc data
DefaultAssay(pbmc.integrated_impute) <- "MAGIC_RNA"
saveRDS(pbmc.integrated_impute,"/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/analyzed/scRNA_nan_ping.RDS")


rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated_impute) = "MAGIC_RNA"
genes <- c("CCR7","IKZF2")
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute,genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

require(gridExtra)
pdf(paste(savedir,"featureplot/genes_tSNE_MAGIC.pdf",sep = ""), width = 9, height = 5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()

rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated_impute) = "MAGIC_RNA"
genes <- c("CCR7","IKZF2")
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute,genes[i], reduction = "umap", x="UMAP_1", y="UMAP_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

require(gridExtra)
pdf(paste(savedir,"featureplot/genes_UMAP_MAGIC.pdf",sep = ""), width = 9, height = 5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()


rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated_impute) = "MAGIC_RNA"
genes <- c("CD247", "FCER1G", "LPP", "LYN","NCAM1","FYN","B3GAT1","IKZF1","IKZF3","IKZF5","TRGC1","TRDC")
genes <- c("CD28")
match(genes,rownames(pbmc.integrated_impute))
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute,genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

require(gridExtra)
pdf(paste(savedir,"genes_tSNE_MAGIC_CD28.pdf",sep = ""), width = 15, height = 12)
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],
             plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],
             plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],ncol=4, nrow=3)
dev.off()

pdf(paste(savedir,"genes_tSNE_MAGIC_CD28.pdf",sep = ""), width = 5.5, height = 5)
p
dev.off()

#### Adding the scCITESeq CD45RA and CD28 Dataset #####
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/analyzed/"
CD45_CD28 <- read.table(paste(savedir,"CD45RA_CD28_data",sep = ""), header = TRUE, sep = "\t")

# meta_cd8 <- meta_cd8 %>% inner_join(meta.data[,c(4,5)], by = 'barcode')
# colnames(meta_cd8)[colnames(meta_cd8) == 'umap1'] <- 'tsne_1'
# colnames(meta_cd8)[colnames(meta_cd8) == 'umap2'] <- 'tsne_2'
# Normalizing the data
CD45_sum <- sum(CD45_CD28$CD45RA)
CD45_CD28$CD45RA_data <- log1p((CD45_CD28$CD45RA/CD45_sum)*10000)
CD28_sum <- sum(CD45_CD28$CD28)
CD45_CD28$CD28_data <- log1p((CD45_CD28$CD28/CD28_sum)*10000)
CD45_CD28$cell_barcodes <- paste(CD45_CD28$Donor_id,CD45_CD28$Age,CD45_CD28$Sex,CD45_CD28$visit,CD45_CD28$Cell.barcode,sep="_")

pbmc.integrated_impute@meta.data$cell_codes <- gsub("-.*","",rownames(pbmc.integrated_impute@meta.data))
pbmc.integrated_impute@meta.data$cell_barcodes <- paste(pbmc.integrated_impute@meta.data$Donor_id,pbmc.integrated_impute@meta.data$Age,pbmc.integrated_impute@meta.data$Sex,pbmc.integrated_impute@meta.data$visit,pbmc.integrated_impute@meta.data$cell_codes,sep="_")
cell_index <- match(CD45_CD28$cell_barcodes, pbmc.integrated_impute@meta.data$cell_barcodes)
pbmc.integrated_impute@meta.data$CD45RA_normalized <- 0
pbmc.integrated_impute@meta.data$CD28_normalized <- 0
pbmc.integrated_impute@meta.data$celltype <- "unknown"

pbmc.integrated_impute@meta.data[cell_index,"CD45RA_normalized"] <- CD45_CD28$CD45RA_data
pbmc.integrated_impute@meta.data[cell_index,"CD28_normalized"] <- CD45_CD28$CD28_data
pbmc.integrated_impute@meta.data[cell_index,"celltype"] <- CD45_CD28$Subset

genes <-  c("CD45RA_normalized", "CD28_normalized")
pdf(paste(savedir,"CD45RA_CD28_featureplot_umap.pdf",sep = ""), width = 10, height = 6)
FeaturePlot(pbmc.integrated_impute,
            reduction = "umap",
            features = genes,
            pt.size = 0.8,
            order = TRUE,
            label = TRUE,
            ncol = 2)
dev.off()

genes <-  c("CD45RA_normalized", "CD28_normalized")
pdf(paste(savedir,"CD45RA_CD28_featureplot_tsne.pdf",sep = ""), width = 10, height = 6)
FeaturePlot(pbmc.integrated_impute,
            reduction = "tsne",
            features = genes,
            pt.size = 0.8,
            order = TRUE,
            label = TRUE,
            ncol = 2)
dev.off()

rm(plot_list)
plot_list <- list()
library(ArchR)
genes <-  c("CD45RA_normalized", "CD28_normalized")
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute,genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}


require(gridExtra)
pdf(paste(savedir,"CD45RA_CD28_featureplot_tsne.pdf",sep = ""), width = 14, height = 6.5)
grid.arrange(plot_list[[1]],plot_list[[2]],ncol=2)
dev.off()


p=featureplot_front(pbmc.integrated_impute,"CD28_normalized", reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits=c(0,2))

pdf(paste(savedir,"CD28.pdf",sep = ""), width = 6, height = 5.5)
p
dev.off()

p2=featureplot_front(pbmc.integrated_impute,"CD45RA_normalized", reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits=c(0,2))

pdf(paste(savedir,"CD45RA.pdf",sep = ""), width = 6, height = 5.5)
p2
dev.off()

require(gridExtra)
pdf(paste(savedir,"CD45RA_CD28_featureplot_tsne.pdf",sep = ""), width = 14, height = 6.5)
grid.arrange(p,p2,ncol=2)
dev.off()

pdf(paste(savedir,"celltypes.pdf",sep = ""), width = 6, height = 5.5)
DimPlot(pbmc.integrated_impute, group.by = "celltype", reduction = "tsne")
dev.off()

cellnames <- rownames(pbmc.integrated_impute@meta.data[grep("unknown",pbmc.integrated_impute@meta.data$celltype, invert=TRUE),])
pbmc.integrated_impute_subset <- subset(pbmc.integrated_impute,cells=cellnames)

pdf(paste(savedir,"celltypes_subset.pdf",sep = ""), width = 6, height = 5.5)
DimPlot(pbmc.integrated_impute_subset, group.by = "celltype", reduction = "tsne", label = TRUE)
dev.off()

pdf(paste(savedir,"scRNA_sc_long_dimplot_tsne.pdf",sep = ""), width = 6, height = 5.5)
DimPlot(pbmc.integrated_impute_subset, reduction = "tsne", label = TRUE)
dev.off()

pdf(paste(savedir,"scRNA_sc_long_dimplot_grey_tsne.pdf",sep = ""), width = 6, height = 5.5)
DimPlot(pbmc.integrated_impute, reduction = "tsne", label = TRUE, cols = c("0" = "red","1" = "grey","2"="grey","3"="grey","4"="grey","5"="grey",
                                                                           "6"="grey","7"="grey","8"="blue","9"="grey","10"="grey","11"="grey",
                                                                           "12"="grey","13"="grey","14"="black","15"="grey","16"="green","17"="grey",
                                                                           "18"="grey","19"="grey"))
dev.off()

# IKZF2 Correlation for cluster 8 and 9
library(Seurat)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/analyzed/IKZF_correlation/"
dir.create(savedir,showWarnings = F)
setwd(savedir)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/scRNA_correlation.R")

### CD4
pbmc.integrated_impute <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/analyzed/scRNA_nan_ping.RDS")
magic_transpose <- t(pbmc.integrated_impute@assays$MAGIC_RNA@data)
genes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/IKZF2_correlation/genes.txt")[,1]
for (i in 1:length(genes)) {
  scpearson_correlation(magic_transpose,genes[i],savedir,filename = "CD4_scCITESeq")
}

# for TEMRA and EM cluster 8,9
cluster <- c(8,9)
celltype <- c("EM","TEMRA")
for (j in 1:length(cluster)) {
  pbmc.integrated_impute_subset <- subset(pbmc.integrated_impute, idents=cluster[j])
  pbmc.integrated_impute_subset_magic_transpose <- t(pbmc.integrated_impute_subset@assays$MAGIC_RNA@data)
  subset_savedir <- paste(savedir,"cluster_",cluster[j],"/",sep = "")
  dir.create(subset_savedir,showWarnings = FALSE)
  scpearson_correlation(pbmc.integrated_impute_subset_magic_transpose,"IKZF2",
                        savedir=subset_savedir,
                        filename = paste("scRNA_",celltype[j],"_",cluster[j],sep = ""))
}

genes <- c("PHACTR2","KLRF1")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute, genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- print(p)
}

require(gridExtra)
pdf(paste(savedir,"PHACTR2_KLRF1.pdf",sep = ""), width = 12, height = 5.5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()

pdf("scRNA_donorid_separated_tsne.pdf",width = 12, height = 12)
DimPlot(pbmc.integrated_impute, split.by = "Donor_id",group.by = "Donor_id", ncol = 6,pt.size = 0.8, reduction = "tsne") + NoLegend()
dev.off()

DefaultAssay(pbmc.integrated_impute) <- "RNA"
EM3_vs_EM1 <- FindMarkers(pbmc.integrated_impute,ident.1="8",ident.2="0")
EM3_vs_EM2 <- FindMarkers(pbmc.integrated_impute,ident.1="8",ident.2="14")

write.table(EM3_vs_EM1,"EM3_vs_EM1.txt",sep = "\t", row.names = T, col.names = T, quote = F)
write.table(EM3_vs_EM2,"EM3_vs_EM2.txt",sep = "\t", row.names = T, col.names = T, quote = F)

#### Raw Counts #####
library(Seurat)
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/"
rawdata = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/"
raw_counts <- readRDS(paste(rawdata,"GSE136184_rsc_Gene_Expression_Aging_UMI_Counts2.Rds",sep = ""))
metadata <- read.table(paste(rawdata,"individual_cell_counts_metadata.txt",sep = ""), header = TRUE)
metadata$study <- metadata$donor_ID
metadata[grep("sc",metadata$study),"study"] <- "cross"
metadata[grep("L",metadata$study),"study"] <- "long"
seu_NP <- readRDS(paste(rawdata,"GSE136184_seu_obj2.Rds",sep = ""))
p1p4 <- CreateSeuratObject(counts = raw_counts, project="Nanping_Object")
p1p4@meta.data$Code <- seu_NP@meta.data$Code
p1p4@meta.data$visit <- seu_NP@meta.data$visit
p1p4@meta.data$bc <- seu_NP@meta.data$bc
p1p4@meta.data$short_bc <- seu_NP@meta.data$short_bc
p1p4@meta.data$Code_visit <- paste(p1p4@meta.data$Code,p1p4@meta.data$visit,sep = "_")
metadata$DonorID <- gsub("-.*","",metadata$donor_ID)
metadata$visit <- gsub(".*_","",metadata$Code_visit)
metadata$DonorID_visit <- paste(metadata$DonorID,metadata$visit,sep = "_")

library(stringi)
p1p4@meta.data$Donor_id <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                                  pattern=metadata$DonorID_visit,
                                                  replacement=metadata$Code_visit,
                                                  vectorize=FALSE)

p1p4@meta.data$Age <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                             pattern=metadata$DonorID_visit,
                                             replacement=metadata$age,
                                             vectorize=FALSE)

p1p4@meta.data$Sex <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                             pattern=metadata$DonorID_visit,
                                             replacement=metadata$sex,
                                             vectorize=FALSE)

p1p4@meta.data$study <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                               pattern=metadata$DonorID_visit,
                                               replacement=metadata$study,
                                               vectorize=FALSE)

p1p4@meta.data$orig.ident <- paste(p1p4@meta.data$Code_visit,p1p4@meta.data$Age,p1p4@meta.data$Sex,p1p4@meta.data$Donor_id,sep="_")

p1p4[["percent.mt"]] <- PercentageFeatureSet(p1p4, pattern = "^MT-")
p1p4$mitoRatio <- p1p4@meta.data$percent.mt / 100

# Changing the Idents
pdf(paste(savedir,"p1p4_QC_mt_percent.pdf",sep = ""), width = 20, height = 7)
VlnPlot(p1p4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.000005, ncol = 3)
dev.off()

pdf(paste(savedir,"p1p4_QC_mt_percent_study.pdf",sep = ""), width = 20, height = 7)
VlnPlot(p1p4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "study", pt.size = 0.000005, ncol = 3)
dev.off()

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_sctransform_V2.R")
objname = "p1p4"
Assay = "RNA"
process = "sctransform"
ngenes=1000
dims=30
regress="mitoRatio"
split_by = "study"
obj=p1p4
reference = 1
p1p4_sctransformed <- sctransform_V2_integration(obj = p1p4, saveDir = savedir, ngenes = 1000,
                                                 regress = c("nFeature_RNA","mitoRatio"),
                                                 dims = 30,
                                                 Assay = Assay, process = process, objname = objname,
                                                 split_by = "study",
                                                 reference = 1,
                                                 sample_tree = NULL)
pbmc.integrated <- p1p4_sctransformed

pbmc.integrated <- RunUMAP(pbmc.integrated, reduction = "pca", dims = 1:20)

## Plot Data
dir.create(paste(savedir,"UMAP/",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"UMAP/dimplot_PCA.pdf",sep = ""),width = 9, height = 6)
DimPlot(pbmc.integrated, reduction = 'pca', group.by = 'study', split.by = "study")
dev.off()

dims <- c(1:20)
perplexities <- c(30)
seed <- c(7456)
iters <- c(3000)
pbmc.integrated <- RunTSNE(method = 'FIt-SNE', object = pbmc.integrated, dims = dims, perplexity = perplexities, seed.use = seed, max_iter = iters)

pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:20)

DefaultAssay(pbmc.integrated) <- "integrated"
dir.create(paste(savedir,"UMAP",sep = ""),showWarnings = FALSE)
res = c(0.2,0.4,0.6,0.8)
res = c(1,1.2,1.4,1.6,1.8)
res = c(1.4)
for (i in 1:length(res)) {
  pbmc.integrated <- FindClusters(pbmc.integrated, resolution = res[i])
  pdf(paste(savedir,"UMAP/dimplot_tsne_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(pbmc.integrated, reduction = 'tsne', label = TRUE, label.size = 5))
  dev.off()

  pdf(paste(savedir,"UMAP/dimplot_umap_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(pbmc.integrated, reduction = 'umap', label = TRUE, label.size = 5))
  dev.off()
}

pdf(paste(savedir,"UMAP/dimplot_tnse.pdf",sep = ""),width = 9, height = 6)
DimPlot(pbmc.integrated, reduction = 'tsne', group.by = 'study', split.by = "study")
dev.off()

pdf(paste(savedir,"UMAP/dimplot_umap.pdf",sep = ""),width = 9, height = 6)
DimPlot(pbmc.integrated, reduction = 'umap', group.by = 'study', split.by = "study")
dev.off()


rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated) = "RNA"
genes <- c("CCR7","IKZF2")
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated,genes[i], reduction = "umap", x="UMAP_1", y="UMAP_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

dir.create(paste(savedir,"featureplot",sep = ""),showWarnings = FALSE)
require(gridExtra)
pdf(paste(savedir,"featureplot/CCR7_IKZF2_genes.pdf",sep = ""), width = 9, height = 5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()

rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated) = "RNA"
genes <- c("CCR7","IKZF2")
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated,genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

require(gridExtra)
pdf(paste(savedir,"featureplot/CCR7_IKZF2_genes_tSNE.pdf",sep = ""), width = 9, height = 5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()

DefaultAssay(pbmc.integrated) <- "RNA"
gene_zero <- as.vector(which(rowSums(pbmc.integrated@assays$RNA@data)==0))
pbmc.integrated_2 <- subset(pbmc.integrated, features = rownames(pbmc.integrated[-gene_zero,]))
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/python/python-3.9.2/bin/python3.9")
library(Rmagic)
py_discover_config("magic") # to check
pbmc.integrated_impute <- magic(pbmc.integrated_2, npca=20) ## imputing the vasc data
DefaultAssay(pbmc.integrated_impute) <- "MAGIC_RNA"

DefaultAssay(pbmc.integrated_impute) <- "integrated"
res = c(1.4)
for (i in 1:length(res)) {
  pbmc.integrated_impute <- FindClusters(pbmc.integrated_impute, resolution = res[i])
  pdf(paste(savedir,"UMAP/dimplot_tsne_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(pbmc.integrated_impute, reduction = 'tsne', label = TRUE, label.size = 5))
  dev.off()

  pdf(paste(savedir,"UMAP/dimplot_umap_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(pbmc.integrated_impute, reduction = 'umap', label = TRUE, label.size = 5))
  dev.off()
}

dir.create(paste(savedir,"saveRDS_obj",sep = ""),showWarnings = FALSE)
saveRDS(pbmc.integrated_impute,paste(savedir,"saveRDS_obj/Nan_ping_CD8_impute.RDS",sep = ""))

rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated) = "MAGIC_RNA"
genes <- c("CCR7","IKZF2")
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute,genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

require(gridExtra)
pdf(paste(savedir,"featureplot/CCR7_IKZF2_tSNE_MAGIC.pdf",sep = ""), width = 9, height = 5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()

rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated) = "MAGIC_RNA"
genes <- c("CCR7","IKZF2")
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute,genes[i], reduction = "umap", x="UMAP_1", y="UMAP_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

require(gridExtra)
pdf(paste(savedir,"featureplot/CCR7_IKZF2_UMAP_MAGIC.pdf",sep = ""), width = 9, height = 5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()

rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated_impute) = "MAGIC_RNA"
genes <- c("CD247", "FCER1G", "LPP", "LYN","NCAM1","FYN","B3GAT1","IKZF1","IKZF3","IKZF5","TRGC1","TRDC")
genes <- c("CD28")
match(genes,rownames(pbmc.integrated_impute))
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute,genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

require(gridExtra)
pdf(paste(savedir,"genes_tSNE_MAGIC_CD28.pdf",sep = ""), width = 15, height = 12)
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],
             plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],
             plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],ncol=4, nrow=3)
dev.off()

pdf(paste(savedir,"genes_tSNE_MAGIC_CD28.pdf",sep = ""), width = 5.5, height = 5)
p
dev.off()

#### Adding the scCITESeq CD45RA and CD28 Dataset #####
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/"
CD45_CD28 <- read.table(paste(savedir,"CD45RA_CD28_data",sep = ""), header = TRUE, sep = "\t")

# meta_cd8 <- meta_cd8 %>% inner_join(meta.data[,c(4,5)], by = 'barcode')
# colnames(meta_cd8)[colnames(meta_cd8) == 'umap1'] <- 'tsne_1'
# colnames(meta_cd8)[colnames(meta_cd8) == 'umap2'] <- 'tsne_2'
# Normalizing the data
CD45_sum <- sum(CD45_CD28$CD45RA)
CD45_CD28$CD45RA_data <- log1p((CD45_CD28$CD45RA/CD45_sum)*10000)
CD28_sum <- sum(CD45_CD28$CD28)
CD45_CD28$CD28_data <- log1p((CD45_CD28$CD28/CD28_sum)*10000)
CD45_CD28$cell_barcodes <- paste(CD45_CD28$Donor_id,
                                 CD45_CD28$Age,
                                 CD45_CD28$Sex,
                                 CD45_CD28$visit,
                                 CD45_CD28$Cell.barcode,sep="_")

CD45_CD28$cell_barcodes <- gsub("-","_",CD45_CD28$cell_barcodes)

pbmc.integrated_impute@meta.data$cell_codes <- gsub("-.*","",rownames(pbmc.integrated_impute@meta.data))
pbmc.integrated_impute@meta.data$cell_barcodes <- paste(pbmc.integrated_impute@meta.data$Code,
                                                        pbmc.integrated_impute@meta.data$visit,
                                                        pbmc.integrated_impute@meta.data$Age,
                                                        pbmc.integrated_impute@meta.data$Sex,
                                                        pbmc.integrated_impute@meta.data$visit,
                                                        pbmc.integrated_impute@meta.data$cell_codes,sep="_")
cell_index <- match(CD45_CD28$cell_barcodes, pbmc.integrated_impute@meta.data$cell_barcodes)
pbmc.integrated_impute@meta.data$CD45RA_normalized <- 0
pbmc.integrated_impute@meta.data$CD28_normalized <- 0
pbmc.integrated_impute@meta.data$celltype <- "unknown"

pbmc.integrated_impute@meta.data[cell_index,"CD45RA_normalized"] <- CD45_CD28$CD45RA_data
pbmc.integrated_impute@meta.data[cell_index,"CD28_normalized"] <- CD45_CD28$CD28_data
pbmc.integrated_impute@meta.data[cell_index,"celltype"] <- CD45_CD28$Subset

genes <-  c("CD45RA_normalized", "CD28_normalized")
pdf(paste(savedir,"featureplot/CD45RA_CD28_featureplot_umap.pdf",sep = ""), width = 10, height = 6)
FeaturePlot(pbmc.integrated_impute,
            reduction = "umap",
            features = genes,
            pt.size = 0.8,
            order = TRUE,
            label = TRUE,
            ncol = 2)
dev.off()

genes <-  c("CD45RA_normalized", "CD28_normalized")
pdf(paste(savedir,"CD45RA_CD28_featureplot_tsne.pdf",sep = ""), width = 10, height = 6)
FeaturePlot(pbmc.integrated_impute,
            reduction = "tsne",
            features = genes,
            pt.size = 0.8,
            order = TRUE,
            label = TRUE,
            ncol = 2)
dev.off()

rm(plot_list)
plot_list <- list()
library(ArchR)
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute,genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}


require(gridExtra)
pdf(paste(savedir,"featureplot/CD45RA_CD28_featureplot_tsne.pdf",sep = ""), width = 14, height = 6.5)
grid.arrange(plot_list[[1]],plot_list[[2]],ncol=2)
dev.off()

p=featureplot_front(pbmc.integrated_impute,"CD28_normalized", reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits=c(0,2))

pdf(paste(savedir,"featureplot/CD28.pdf",sep = ""), width = 6, height = 5.5)
p
dev.off()

p2=featureplot_front(pbmc.integrated_impute,"CD45RA_normalized", reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits=c(0,2))

pdf(paste(savedir,"featureplot/CD45RA.pdf",sep = ""), width = 6, height = 5.5)
p2
dev.off()

require(gridExtra)
pdf(paste(savedir,"CD45RA_CD28_featureplot_tsne.pdf",sep = ""), width = 14, height = 6.5)
grid.arrange(p,p2,ncol=2)
dev.off()

cellnames <- rownames(pbmc.integrated_impute@meta.data[grep("unknown",pbmc.integrated_impute@meta.data$celltype, invert=TRUE),])
pbmc.integrated_impute_subset <- subset(pbmc.integrated_impute,cells=cellnames)

pdf(paste(savedir,"celltypes_subset.pdf",sep = ""), width = 6, height = 5.5)
DimPlot(pbmc.integrated_impute_subset, group.by = "celltype", reduction = "tsne", label = TRUE)
dev.off()

pdf(paste(savedir,"scRNA_sc_long_dimplot_tsne.pdf",sep = ""), width = 6, height = 5.5)
DimPlot(pbmc.integrated_impute_subset, reduction = "tsne", label = TRUE)
dev.off()

DefaultAssay(pbmc.integrated_impute) <- "MAGIC_RNA"
pdf(paste(savedir,"Vlnplot/IKZF2_MAGIC_boxplot.pdf",sep = ""),width = 10, height = 6)
VlnPlot(pbmc.integrated_impute, features = "IKZF2",pt.size = 0) + geom_boxplot()
dev.off()

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/"
dir.create(paste(savedir,"saveRDS_obj",sep = ""),showWarnings = FALSE)
saveRDS(pbmc.integrated_impute,paste(savedir,"saveRDS_obj/Nan_ping_CD8_impute.RDS",sep = ""))

# IKZF2 Correlation for cluster 8 17, 18
library(Seurat)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/IKZF_correlation/"
dir.create(savedir,showWarnings = F)
setwd(savedir)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/scRNA_correlation.R")

### CD4
# pbmc.integrated_impute <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/analyzed/scRNA_nan_ping.RDS")
# magic_transpose <- t(pbmc.integrated_impute@assays$MAGIC_RNA@data)
# genes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/IKZF2_correlation/genes.txt")[,1]
# for (i in 1:length(genes)) {
#   scpearson_correlation(magic_transpose,genes[i],savedir,filename = "CD4_scCITESeq")
# }

# for TEMRA and EM cluster 8,17,18,7,23
cluster <- c(8,17,18,7,23)
celltype <- c("TEMRA","TEMRA","TEMRA","EM","EM")
for (j in 1:length(cluster)) {
  pbmc.integrated_impute_subset <- subset(pbmc.integrated_impute, idents=cluster[j])
  pbmc.integrated_impute_subset_magic_transpose <- t(pbmc.integrated_impute_subset@assays$MAGIC_RNA@data)
  subset_savedir <- paste(savedir,"cluster_",cluster[j],"/",sep = "")
  dir.create(subset_savedir,showWarnings = FALSE)
  scpearson_correlation(pbmc.integrated_impute_subset_magic_transpose,"IKZF2",
                        savedir=subset_savedir,
                        filename = paste("scRNA_",celltype[j],"_",cluster[j],sep = ""))
}

genes <- c("PHACTR2","KLRF1")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute, genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- print(p)
}

require(gridExtra)
pdf(paste(savedir,"PHACTR2_KLRF1.pdf",sep = ""), width = 12, height = 5.5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()

pdf("scRNA_donorid_separated_tsne.pdf",width = 12, height = 12)
DimPlot(pbmc.integrated_impute, split.by = "Donor_id",group.by = "Donor_id", ncol = 6,pt.size = 0.8, reduction = "tsne") + NoLegend()
dev.off()

DefaultAssay(pbmc.integrated_impute) <- "RNA"
EM3_vs_EM1 <- FindMarkers(pbmc.integrated_impute,ident.1="8",ident.2="0")
EM3_vs_EM2 <- FindMarkers(pbmc.integrated_impute,ident.1="8",ident.2="14")

write.table(EM3_vs_EM1,"EM3_vs_EM1.txt",sep = "\t", row.names = T, col.names = T, quote = F)
write.table(EM3_vs_EM2,"EM3_vs_EM2.txt",sep = "\t", row.names = T, col.names = T, quote = F)

### Raw Count following their code ####
## We are following their
library(Seurat)
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/Nan_ping_code/"
dir.create(savedir, showWarnings = FALSE)
rawdata = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/"
raw_counts <- readRDS(paste(rawdata,"GSE136184_rsc_Gene_Expression_Aging_UMI_Counts2.Rds",sep = ""))
seu <- readRDS(paste(rawdata,"GSE136184_seu_obj2.Rds",sep = ""))
gene_req <- rownames(seu)
raw_counts_gene_req <- raw_counts[match(gene_req,rownames(raw_counts)),]

metadata <- read.table(paste(rawdata,"individual_cell_counts_metadata.txt",sep = ""), header = TRUE)
metadata$study <- metadata$donor_ID
metadata[grep("sc",metadata$study),"study"] <- "cross"
metadata[grep("L",metadata$study),"study"] <- "long"

seu_NP <- seu
p1p4 <- CreateSeuratObject(counts = raw_counts_gene_req, project="Nanping_Object")
p1p4@meta.data$Code <- seu_NP@meta.data$Code
p1p4@meta.data$visit <- seu_NP@meta.data$visit
p1p4@meta.data$bc <- seu_NP@meta.data$bc
p1p4@meta.data$short_bc <- seu_NP@meta.data$short_bc
p1p4@meta.data$Code_visit <- paste(p1p4@meta.data$Code,p1p4@meta.data$visit,sep = "_")
metadata$DonorID <- gsub("-.*","",metadata$donor_ID)
metadata$visit <- gsub(".*_","",metadata$Code_visit)
metadata$DonorID_visit <- paste(metadata$DonorID,metadata$visit,sep = "_")

library(stringi)
p1p4@meta.data$Donor_id <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                                  pattern=metadata$DonorID_visit,
                                                  replacement=metadata$Code_visit,
                                                  vectorize=FALSE)

p1p4@meta.data$Age <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                             pattern=metadata$DonorID_visit,
                                             replacement=metadata$age,
                                             vectorize=FALSE)

p1p4@meta.data$Sex <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                             pattern=metadata$DonorID_visit,
                                             replacement=metadata$sex,
                                             vectorize=FALSE)

p1p4@meta.data$study <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                               pattern=metadata$DonorID_visit,
                                               replacement=metadata$study,
                                               vectorize=FALSE)

p1p4@meta.data$orig.ident <- paste(p1p4@meta.data$Code_visit,p1p4@meta.data$Age,p1p4@meta.data$Sex,p1p4@meta.data$Donor_id,sep="_")

p1p4[["percent.mt"]] <- PercentageFeatureSet(p1p4, pattern = "^MT-")
p1p4$mitoRatio <- p1p4@meta.data$percent.mt / 100

# Changing the Idents
pdf(paste(savedir,"p1p4_QC_mt_percent.pdf",sep = ""), width = 20, height = 7)
VlnPlot(p1p4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.000005, ncol = 3)
dev.off()

pdf(paste(savedir,"p1p4_QC_mt_percent_study.pdf",sep = ""), width = 20, height = 7)
VlnPlot(p1p4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "study", pt.size = 0.000005, ncol = 3)
dev.off()

p1p4.list <- SplitObject(p1p4, split.by = "study")

for (i in 1:length(p1p4.list)) {
  p1p4.list[[i]] <- NormalizeData(p1p4.list[[i]], verbose = T)
  p1p4.list[[i]] <- FindVariableFeatures(p1p4.list[[i]], verbose = T, nfeatures = 2000)
}

# savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/Nan_ping_code_log_normalize/"
# reference_dataset <- which(names(p1p4.list) == "cross")
# features <- SelectIntegrationFeatures(object.list = p1p4.list)
# pbmc.anchors <- FindIntegrationAnchors(object.list = p1p4.list, dims = 1:30,
#                                        reference = reference_dataset, anchor.features = features)
# pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, dims = 1:30)
# pbmc.integrated <- ScaleData(pbmc.integrated)
# pbmc.integrated <- RunPCA(object = pbmc.integrated, verbose = FALSE)
# pbmc.integrated <- RunUMAP(pbmc.integrated, reduction = "pca", dims = 1:30)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_sctransform_V2.R")
objname = "p1p4"
Assay = "RNA"
process = "sctransform"
ngenes=2000
dims=30
regress="mitoRatio"
split_by = "study"
obj=p1p4
reference = 1
p1p4_sctransformed <- sctransform_V2_integration(obj = p1p4, saveDir = savedir, ngenes = 2000,
                                                 regress = c("nFeature_RNA","mitoRatio"),
                                                 dims = 30,
                                                 Assay = Assay, process = process,
                                                 objname = objname,
                                                 split_by = "study",
                                                 reference = 1,
                                                 sample_tree = NULL)
pbmc.integrated <- p1p4_sctransformed
pbmc.integrated <- RunUMAP(pbmc.integrated, reduction = "pca", dims = 1:30)

## Plot Data
dir.create("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/Nan_ping_code_log_normalize/",showWarnings = FALSE)
dir.create(paste(savedir,"UMAP/",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"UMAP/dimplot_UMAP_sctransform.pdf",sep = ""),width = 9, height = 6)
DimPlot(pbmc.integrated, reduction = 'umap', group.by = 'study', split.by = "study")
dev.off()

dims <- c(1:30)
perplexities <- c(30)
seed <- c(7456)
iters <- c(3000)
pbmc.integrated <- RunTSNE(method = 'FIt-SNE', object = pbmc.integrated, dims = dims, perplexity = perplexities, seed.use = seed, max_iter = iters)
dir.create(paste(savedir,"UMAP/",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"UMAP/dimplot_tSNE_sctransform.pdf",sep = ""),width = 9, height = 6)
DimPlot(pbmc.integrated, reduction = 'tsne', group.by = 'study', split.by = "study")
dev.off()

pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:30)

dir.create(paste(savedir,"UMAP/",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"UMAP/dimplot_tSNE_sctransform.pdf",sep = ""),width = 9, height = 6)
DimPlot(pbmc.integrated, reduction = 'tsne', group.by = 'study', split.by = "study")
dev.off()

DefaultAssay(pbmc.integrated) <- "integrated"
dir.create(paste(savedir,"UMAP",sep = ""),showWarnings = FALSE)
res = c(0.2,0.4,0.6,0.8,1,1.2)
for (i in 1:length(res)) {
  pbmc.integrated <- FindClusters(pbmc.integrated, resolution = res[i])
  pdf(paste(savedir,"UMAP/dimplot_tsne_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(pbmc.integrated, reduction = 'tsne', label = TRUE, label.size = 5))
  dev.off()

  pdf(paste(savedir,"UMAP/dimplot_umap_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(pbmc.integrated, reduction = 'umap', label = TRUE, label.size = 5))
  dev.off()
}

rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated) = "MAGIC_RNA"
genes <- c("CCR7","IKZF2")
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated,genes[i], reduction = "umap", x="UMAP_1", y="UMAP_2",
                      size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

dir.create(paste(savedir,"featureplot",sep = ""),showWarnings = FALSE)
require(gridExtra)
pdf(paste(savedir,"featureplot/CCR7_IKZF2_genes_MAGIC.pdf",sep = ""), width = 9, height = 5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()

rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated) = "RNA"
genes <- c("CCR7","IKZF2")
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated,genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

require(gridExtra)
pdf(paste(savedir,"featureplot/CCR7_IKZF2_genes_tSNE.pdf",sep = ""), width = 9, height = 5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()

DefaultAssay(pbmc.integrated) <- "RNA"
gene_zero <- as.vector(which(rowSums(pbmc.integrated@assays$RNA@data)==0))
pbmc.integrated_2 <- subset(pbmc.integrated, features = rownames(pbmc.integrated[-gene_zero,]))
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/python/python-3.9.2/bin/python3.9")
library(Rmagic)
py_discover_config("magic") # to check
pbmc.integrated_impute <- magic(pbmc.integrated_2, npca=20) ## imputing the vasc data
DefaultAssay(pbmc.integrated_impute) <- "MAGIC_RNA"

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/Nan_ping_code/"
dir.create(paste(savedir,"saveRDS_obj",sep = ""),showWarnings = FALSE)
saveRDS(pbmc.integrated_impute,paste(savedir,"saveRDS_obj/Nan_ping_CD8_impute_sctransform.RDS",sep = ""))

rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated_impute) = "MAGIC_RNA"
genes <- c("CCR7","IKZF2")
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute,genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

require(gridExtra)
pdf(paste(savedir,"featureplot/CCR7_IKZF2_tSNE_MAGIC.pdf",sep = ""), width = 9, height = 5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()

rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated_impute) = "MAGIC_RNA"
genes <- c("CCR7","IKZF2")
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute,genes[i], reduction = "umap", x="UMAP_1", y="UMAP_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

require(gridExtra)
pdf(paste(savedir,"featureplot/CCR7_IKZF2_UMAP_MAGIC.pdf",sep = ""), width = 9, height = 5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/Nan_ping_code/"
DefaultAssay(pbmc.integrated_impute) <- "MAGIC_RNA"
pdf(paste(savedir,"featureplot/IKZF2_UMAP.pdf",sep = ""), width = 5.5, height = 5)
featureplot_front(pbmc.integrated_impute,"IKZF2", reduction = "umap", x="UMAP_1", y="UMAP_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits=c(0,0.20))
dev.off()

DefaultAssay(pbmc.integrated_impute) <- "MAGIC_RNA"
pdf(paste(savedir,"featureplot/IKZF2.pdf",sep = ""), width = 5.5, height = 5)
featureplot_front(pbmc.integrated_impute,"IKZF2", reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits=c(0,0.20))
dev.off()


rm(plot_list)
plot_list = list()
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(pbmc.integrated_impute) = "MAGIC_RNA"
genes <- c("CD247", "FCER1G", "LPP", "LYN","NCAM1","FYN","B3GAT1","IKZF1","IKZF3","IKZF5","TRGC1","TRDC")
genes <- c("CD28")
match(genes,rownames(pbmc.integrated_impute))
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute,genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

require(gridExtra)
pdf(paste(savedir,"genes_tSNE_MAGIC_CD28.pdf",sep = ""), width = 15, height = 12)
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],
             plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],
             plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],ncol=4, nrow=3)
dev.off()

pdf(paste(savedir,"genes_tSNE_MAGIC_CD28.pdf",sep = ""), width = 5.5, height = 5)
p
dev.off()

#### Adding the scCITESeq CD45RA and CD28 Dataset #####
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/Nan_ping_code/"
CD45_CD28 <- read.table(paste(savedir,"CD45RA_CD28_data",sep = ""), header = TRUE, sep = "\t")

# meta_cd8 <- meta_cd8 %>% inner_join(meta.data[,c(4,5)], by = 'barcode')
# colnames(meta_cd8)[colnames(meta_cd8) == 'umap1'] <- 'tsne_1'
# colnames(meta_cd8)[colnames(meta_cd8) == 'umap2'] <- 'tsne_2'
# Normalizing the data
CD45_sum <- sum(CD45_CD28$CD45RA)
CD45_CD28$CD45RA_data <- log1p((CD45_CD28$CD45RA/CD45_sum)*10000)
CD28_sum <- sum(CD45_CD28$CD28)
CD45_CD28$CD28_data <- log1p((CD45_CD28$CD28/CD28_sum)*10000)
CD45_CD28$cell_barcodes <- paste(CD45_CD28$Donor_id,
                                 CD45_CD28$Age,
                                 CD45_CD28$Sex,
                                 CD45_CD28$visit,
                                 CD45_CD28$Cell.barcode,sep="_")

CD45_CD28$cell_barcodes <- gsub("-","_",CD45_CD28$cell_barcodes)

pbmc.integrated_impute@meta.data$cell_codes <- gsub("-.*","",rownames(pbmc.integrated_impute@meta.data))
pbmc.integrated_impute@meta.data$cell_barcodes <- paste(pbmc.integrated_impute@meta.data$Code,
                                                        pbmc.integrated_impute@meta.data$visit,
                                                        pbmc.integrated_impute@meta.data$Age,
                                                        pbmc.integrated_impute@meta.data$Sex,
                                                        pbmc.integrated_impute@meta.data$visit,
                                                        pbmc.integrated_impute@meta.data$cell_codes,sep="_")
cell_index <- match(CD45_CD28$cell_barcodes, pbmc.integrated_impute@meta.data$cell_barcodes)
pbmc.integrated_impute@meta.data$CD45RA_normalized <- 0
pbmc.integrated_impute@meta.data$CD28_normalized <- 0
pbmc.integrated_impute@meta.data$celltype <- "unknown"

pbmc.integrated_impute@meta.data[cell_index,"CD45RA_normalized"] <- CD45_CD28$CD45RA_data
pbmc.integrated_impute@meta.data[cell_index,"CD28_normalized"] <- CD45_CD28$CD28_data
pbmc.integrated_impute@meta.data[cell_index,"celltype"] <- CD45_CD28$Subset

genes <-  c("CD45RA_normalized", "CD28_normalized")
pdf(paste(savedir,"featureplot/CD45RA_CD28_featureplot_umap.pdf",sep = ""), width = 10, height = 6)
FeaturePlot(pbmc.integrated_impute,
            reduction = "umap",
            features = genes,
            pt.size = 0.8,
            order = TRUE,
            label = TRUE,
            ncol = 2)
dev.off()

genes <-  c("CD45RA_normalized", "CD28_normalized")
pdf(paste(savedir,"featureplot/CD45RA_CD28_featureplot_tsne.pdf",sep = ""), width = 10, height = 6)
FeaturePlot(pbmc.integrated_impute,
            reduction = "tsne",
            features = genes,
            pt.size = 0.8,
            order = TRUE,
            label = TRUE,
            ncol = 2)
dev.off()

rm(plot_list)
plot_list <- list()
library(ArchR)
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute,genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}


require(gridExtra)
pdf(paste(savedir,"featureplot/CD45RA_CD28_featureplot_tsne.pdf",sep = ""), width = 14, height = 6.5)
grid.arrange(plot_list[[1]],plot_list[[2]],ncol=2)
dev.off()

p=featureplot_front(pbmc.integrated_impute,"CD28_normalized", reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits=c(0,2))

pdf(paste(savedir,"featureplot/CD28.pdf",sep = ""), width = 6, height = 5.5)
p
dev.off()

p2=featureplot_front(pbmc.integrated_impute,"CD45RA_normalized", reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits=c(0,2))

pdf(paste(savedir,"featureplot/CD45RA.pdf",sep = ""), width = 6, height = 5.5)
p2
dev.off()

require(gridExtra)
pdf(paste(savedir,"CD45RA_CD28_featureplot_tsne.pdf",sep = ""), width = 14, height = 6.5)
grid.arrange(p,p2,ncol=2)
dev.off()

cellnames <- rownames(pbmc.integrated_impute@meta.data[grep("unknown",pbmc.integrated_impute@meta.data$celltype, invert=TRUE),])
pbmc.integrated_impute_subset <- subset(pbmc.integrated_impute,cells=cellnames)

pdf(paste(savedir,"celltypes_subset.pdf",sep = ""), width = 6, height = 5.5)
DimPlot(pbmc.integrated_impute_subset, group.by = "celltype", reduction = "tsne", label = TRUE)
dev.off()

pdf(paste(savedir,"scRNA_sc_long_dimplot_tsne.pdf",sep = ""), width = 6, height = 5.5)
DimPlot(pbmc.integrated_impute_subset, reduction = "tsne", label = TRUE)
dev.off()

dir.create(paste(savedir,"Vlnplot/",sep = ""))
DefaultAssay(pbmc.integrated_impute) <- "MAGIC_RNA"


DefaultAssay(pbmc.integrated_impute) <- "integrated"
for (i in 1:length(res)) {
  DefaultAssay(pbmc.integrated_impute) <- "integrated"
  pbmc.integrated_impute <- FindClusters(pbmc.integrated_impute, resolution = res[i])
  DefaultAssay(pbmc.integrated_impute) <- "MAGIC_RNA"
  pdf(paste(savedir,"Vlnplot/IKZF2_MAGIC_boxplot_",res[i],".pdf",sep = ""),width = 10, height = 6)
  print(VlnPlot(pbmc.integrated_impute, features = "IKZF2",pt.size = 0) + geom_boxplot())
  dev.off()
}


savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/"
dir.create(paste(savedir,"saveRDS_obj",sep = ""),showWarnings = FALSE)
saveRDS(pbmc.integrated_impute_subset,paste(savedir,"saveRDS_obj/Nan_ping_CD8_impute_subset_sctransform.RDS",sep = ""))

saveRDS(pbmc.integrated_impute,paste(savedir,"saveRDS_obj/Nan_ping_CD8_impute_subset.RDS",sep = ""))

# IKZF2 Correlation for cluster 8 17, 18
library(Seurat)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/IKZF_correlation/"
dir.create(savedir,showWarnings = F)
setwd(savedir)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/scRNA_correlation.R")

### CD4
# pbmc.integrated_impute <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/analyzed/scRNA_nan_ping.RDS")
# magic_transpose <- t(pbmc.integrated_impute@assays$MAGIC_RNA@data)
# genes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/IKZF2_correlation/genes.txt")[,1]
# for (i in 1:length(genes)) {
#   scpearson_correlation(magic_transpose,genes[i],savedir,filename = "CD4_scCITESeq")
# }

# for TEMRA and EM cluster 8,17,18,7,23
cluster <- c(8,17,18,7,23)
celltype <- c("TEMRA","TEMRA","TEMRA","EM","EM")
for (j in 1:length(cluster)) {
  pbmc.integrated_impute_subset <- subset(pbmc.integrated_impute, idents=cluster[j])
  pbmc.integrated_impute_subset_magic_transpose <- t(pbmc.integrated_impute_subset@assays$MAGIC_RNA@data)
  subset_savedir <- paste(savedir,"cluster_",cluster[j],"/",sep = "")
  dir.create(subset_savedir,showWarnings = FALSE)
  scpearson_correlation(pbmc.integrated_impute_subset_magic_transpose,"IKZF2",
                        savedir=subset_savedir,
                        filename = paste("scRNA_",celltype[j],"_",cluster[j],sep = ""))
}

genes <- c("PHACTR2","KLRF1")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p=featureplot_front(pbmc.integrated_impute, genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- print(p)
}

require(gridExtra)
pdf(paste(savedir,"PHACTR2_KLRF1.pdf",sep = ""), width = 12, height = 5.5)
grid.arrange(plot_list[[1]],plot_list[[2]], ncol =2)
dev.off()

pdf("scRNA_donorid_separated_tsne.pdf",width = 12, height = 12)
DimPlot(pbmc.integrated_impute, split.by = "Donor_id",group.by = "Donor_id", ncol = 6,pt.size = 0.8, reduction = "tsne") + NoLegend()
dev.off()

DefaultAssay(pbmc.integrated_impute) <- "RNA"
EM3_vs_EM1 <- FindMarkers(pbmc.integrated_impute,ident.1="8",ident.2="0")
EM3_vs_EM2 <- FindMarkers(pbmc.integrated_impute,ident.1="8",ident.2="14")

write.table(EM3_vs_EM1,"EM3_vs_EM1.txt",sep = "\t", row.names = T, col.names = T, quote = F)
write.table(EM3_vs_EM2,"EM3_vs_EM2.txt",sep = "\t", row.names = T, col.names = T, quote = F)


### Extracting out the cross sectional ######
## We will perform the analysis on the cross section and will reference mapped to the longitudinal so that CD45RA label or TEMRA is transferred
# on to the RNA
library(Seurat)
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/"
dir.create(savedir, showWarnings = FALSE)
rawdata = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/"
raw_counts <- readRDS(paste(rawdata,"GSE136184_rsc_Gene_Expression_Aging_UMI_Counts2.Rds",sep = ""))
seu <- readRDS(paste(rawdata,"GSE136184_seu_obj2.Rds",sep = ""))
gene_req <- rownames(seu)
raw_counts_gene_req <- raw_counts[match(gene_req,rownames(raw_counts)),]

metadata <- read.table(paste(rawdata,"individual_cell_counts_metadata.txt",sep = ""), header = TRUE)
metadata$study <- metadata$donor_ID
metadata[grep("sc",metadata$study),"study"] <- "cross"
metadata[grep("L",metadata$study),"study"] <- "long"

seu_NP <- seu
p1p4 <- CreateSeuratObject(counts = raw_counts_gene_req, project="Nanping_Object")
p1p4@meta.data$Code <- seu_NP@meta.data$Code
p1p4@meta.data$visit <- seu_NP@meta.data$visit
p1p4@meta.data$bc <- seu_NP@meta.data$bc
p1p4@meta.data$short_bc <- seu_NP@meta.data$short_bc
p1p4@meta.data$Code_visit <- paste(p1p4@meta.data$Code,p1p4@meta.data$visit,sep = "_")
metadata$DonorID <- gsub("-.*","",metadata$donor_ID)
metadata$visit <- gsub(".*_","",metadata$Code_visit)
metadata$DonorID_visit <- paste(metadata$DonorID,metadata$visit,sep = "_")

library(stringi)
p1p4@meta.data$Donor_id <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                                  pattern=metadata$DonorID_visit,
                                                  replacement=metadata$Code_visit,
                                                  vectorize=FALSE)

p1p4@meta.data$Age <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                             pattern=metadata$DonorID_visit,
                                             replacement=metadata$age,
                                             vectorize=FALSE)

p1p4@meta.data$Sex <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                             pattern=metadata$DonorID_visit,
                                             replacement=metadata$sex,
                                             vectorize=FALSE)

p1p4@meta.data$study <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                               pattern=metadata$DonorID_visit,
                                               replacement=metadata$study,
                                               vectorize=FALSE)

p1p4@meta.data$orig.ident <- paste(p1p4@meta.data$Code_visit,p1p4@meta.data$Age,p1p4@meta.data$Sex,
                                   p1p4@meta.data$Donor_id,sep="_")

p1p4[["percent.mt"]] <- PercentageFeatureSet(p1p4, pattern = "^MT-")
p1p4$mitoRatio <- p1p4@meta.data$percent.mt / 100

# Changing the Idents
pdf(paste(savedir,"p1p4_QC_mt_percent.pdf",sep = ""), width = 20, height = 7)
VlnPlot(p1p4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.000005, ncol = 3)
dev.off()

pdf(paste(savedir,"p1p4_QC_mt_percent_study.pdf",sep = ""), width = 20, height = 7)
VlnPlot(p1p4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "study", pt.size = 0.000005, ncol = 3)
dev.off()

cross_cellnames = rownames(p1p4@meta.data)[grep("cross",p1p4@meta.data$study)]
cross_seu <- subset(p1p4,cells=cross_cellnames)

cross_seu <- NormalizeData(cross_seu, normalization.method = "LogNormalize", scale.factor = 10000)
cross_seu <- FindVariableFeatures(cross_seu, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cross_seu)
cross_seu <- ScaleData(cross_seu, features = all.genes)
cross_seu <- RunPCA(cross_seu, features = VariableFeatures(object = cross_seu))
pdf(paste(savedir,"Elbow_plot.pdf",sep = ""))
ElbowPlot(cross_seu, ndims = 50)
dev.off()

dims <- c(1:20)
seed <- c(7456)
cross_seu <- RunTSNE(object = cross_seu, dims = dims,seed.use = seed)

cross_seu <- RunUMAP(cross_seu, dims = 1:20)
cross_seu <- FindNeighbors(cross_seu, dims = 1:20)

DefaultAssay(cross_seu) <- "RNA"
gene_zero <- as.vector(which(rowSums(cross_seu@assays$RNA@data)==0))
cross_seu_2 <- subset(cross_seu, features = rownames(cross_seu[-gene_zero,]))
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/python/python-3.9.2/bin/python3.9")
library(Rmagic)
py_discover_config("magic") # to check
cross_seu_impute <- magic(cross_seu_2, npca=20) ## imputing the vasc data
DefaultAssay(cross_seu_impute) <- "MAGIC_RNA"
dir.create(paste(savedir,"saveRDS_obj/",sep = ""))
saveRDS(cross_seu_impute,paste(savedir,"saveRDS_obj/cross_sectional_impute.RDS",sep = ""))

dir.create(paste(savedir,"featureplot",sep = ""),showWarnings = FALSE)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(cross_seu_impute) = "MAGIC_RNA"
p1 <- featureplot_front(cross_seu_impute,"IKZF2", reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits = c(0,0.3))
p2 <- featureplot_front(cross_seu_impute,"IKZF2", reduction = "umap", x="UMAP_1", y="UMAP_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits = c(0,0.3))
pdf(paste(savedir,"featureplot/IKZF2_tSNE_UMAP.pdf",sep = ""), width = 9, height = 5)
p1+p2
dev.off()

dir.create(paste(savedir,"UMAP",sep = ""),showWarnings = FALSE)
dir.create(paste(savedir,"Vlnplot",sep = ""),showWarnings = FALSE)

res = c(0.2,0.4,0.6,0.8,1,1.2)
for (i in 1:length(res)) {
  DefaultAssay(cross_seu_impute) <- "RNA"
  cross_seu_impute <- FindClusters(cross_seu_impute, resolution = res[i])
  pdf(paste(savedir,"UMAP/dimplot_tsne_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(cross_seu_impute, reduction = 'tsne', label = TRUE, label.size = 5))
  dev.off()

  pdf(paste(savedir,"UMAP/dimplot_umap_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(cross_seu_impute, reduction = 'umap', label = TRUE, label.size = 5))
  dev.off()

  DefaultAssay(cross_seu_impute) <- "RNA"
  pdf(paste(savedir,"Vlnplot/IKZF2_boxplot_",res[i],".pdf",sep = ""),width = 10, height = 6)
  print(VlnPlot(cross_seu_impute, features = "IKZF2",pt.size = 0) + geom_boxplot())
  dev.off()

  DefaultAssay(cross_seu_impute) <- "MAGIC_RNA"
  pdf(paste(savedir,"Vlnplot/IKZF2_MAGIC_boxplot_",res[i],".pdf",sep = ""),width = 10, height = 6)
  print(VlnPlot(cross_seu_impute, features = "IKZF2",pt.size = 0) + geom_boxplot())
  dev.off()
}

res = c(0.2,0.4,0.6,0.8,1,1.2)
for (i in 1:length(res)) {
  DefaultAssay(cross_seu_impute) <- "RNA"
  cross_seu_impute <- FindClusters(cross_seu_impute, resolution = res[i])

  pdf(paste(savedir,"UMAP/dimplot_pc_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(cross_seu_impute, reduction = 'pca', label = TRUE, label.size = 5))
  dev.off()
}

### Makling longitudinal as reference ####
long_cellnames = rownames(p1p4@meta.data)[grep("long",p1p4@meta.data$study)]
# need to remove cell that does not have any subset
long <- subset(p1p4,cells=long_cellnames)
## Adding the subset detail
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/Nan_ping_code/"
CD45_CD28 <- read.table(paste(maindir,"CD45RA_CD28_data",sep = ""), header = TRUE, sep = "\t")

# Normalizing the data
CD45_sum <- sum(CD45_CD28$CD45RA)
CD45_CD28$CD45RA_data <- log1p((CD45_CD28$CD45RA/CD45_sum)*10000)
CD28_sum <- sum(CD45_CD28$CD28)
CD45_CD28$CD28_data <- log1p((CD45_CD28$CD28/CD28_sum)*10000)
CD45_CD28$cell_barcodes <- paste(CD45_CD28$Donor_id,
                                 CD45_CD28$Age,
                                 CD45_CD28$Sex,
                                 CD45_CD28$visit,
                                 CD45_CD28$Cell.barcode,sep="_")

CD45_CD28$cell_barcodes <- gsub("-","_",CD45_CD28$cell_barcodes)

long@meta.data$cell_codes <- gsub("-.*","",rownames(long@meta.data))
long@meta.data$cell_barcodes <- paste(long@meta.data$Code,long@meta.data$visit,
                                      long@meta.data$Age,long@meta.data$Sex,
                                      long@meta.data$visit,long@meta.data$cell_codes,
                                      sep="_")
cell_index <- match(CD45_CD28$cell_barcodes, long@meta.data$cell_barcodes)
long@meta.data$CD45RA_normalized <- 0
long@meta.data$CD28_normalized <- 0
long@meta.data$celltype <- "unknown"

long@meta.data[cell_index,"CD45RA_normalized"] <- CD45_CD28$CD45RA_data
long@meta.data[cell_index,"CD28_normalized"] <- CD45_CD28$CD28_data
long@meta.data[cell_index,"celltype"] <- CD45_CD28$Subset

long_cellnames_2 <- rownames(long@meta.data)[grep("unknown",long@meta.data$celltype, invert=TRUE)]
long_2 <- subset(long,cells=long_cellnames_2)

long_2 <- NormalizeData(long_2, normalization.method = "LogNormalize", scale.factor = 10000)
long_2 <- FindVariableFeatures(long_2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(long_2)
long_2 <- ScaleData(long_2, features = all.genes)
long_2 <- RunPCA(long_2, features = VariableFeatures(object = long_2))
pdf(paste(savedir,"Elbow_long_2_plot.pdf",sep = ""))
ElbowPlot(long_2, ndims = 50)
dev.off()

long_2 <- RunUMAP(long_2, dims = 1:20)
dims <- c(1:20)
seed <- c(7456)
long_2 <- RunTSNE(object = long_2, dims = dims,seed.use = seed)
long_2 <- FindNeighbors(long_2, dims = 1:20)
long_2 <- FindClusters(long_2, resolution=0.6)

pdf(paste(savedir,"UMAP/dimplot_tsne_long_2.pdf",sep = ""),width = 6, height = 5.5)
print(DimPlot(long_2, reduction = 'tsne', label = TRUE, label.size = 5))
dev.off()

pdf(paste(savedir,"UMAP/dimplot_umap_long_2.pdf",sep = ""),width = 6, height = 5.5)
print(DimPlot(long_2, reduction = 'umap', label = TRUE, label.size = 5))
dev.off()

pdf(paste(savedir,"UMAP/dimplot_tsne_long_2_cell_subsets.pdf",sep = ""),width = 6, height = 5.5)
print(DimPlot(long_2, reduction = 'tsne', label = TRUE, label.size = 5, group.by = "celltype"))
dev.off()

pdf(paste(savedir,"UMAP/dimplot_umap_long_2_cell_subsets.pdf",sep = ""),width = 6, height = 5.5)
print(DimPlot(long_2, reduction = 'umap', label = TRUE, label.size = 5, group.by = "celltype"))
dev.off()

rm(plot_list)
plot_list <- list()
library(ArchR)
genes <-  c("CD45RA_normalized", "CD28_normalized")
for (i in 1:length(genes)) {
  p=featureplot_front(long_2,genes[i], reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01,color=ArchRPalettes$solarExtra)
  plot_list[[i]]=print(p)
}

require(gridExtra)
pdf(paste(savedir,"CD45RA_CD28_featureplot_tsne.pdf",sep = ""), width = 14, height = 6.5)
grid.arrange(plot_list[[1]],plot_list[[2]],ncol=2)
dev.off()

p=featureplot_front(long_2,"CD28_normalized", reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits=c(0,2))

pdf(paste(savedir,"featureplot/CD28.pdf",sep = ""), width = 6, height = 5.5)
p
dev.off()

p2=featureplot_front(long_2,"CD45RA_normalized", reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits=c(0,2))

pdf(paste(savedir,"featureplot/CD45RA.pdf",sep = ""), width = 6, height = 5.5)
p2
dev.off()

require(gridExtra)
pdf(paste(savedir,"CD45RA_CD28_featureplot_tsne.pdf",sep = ""), width = 14, height = 6.5)
grid.arrange(p,p2,ncol=2)
dev.off()

### Reference mapping longitudinal on cross sectional ####
reference <- long_2
query <- cross_seu_impute
DefaultAssay(query) = "RNA"

query.anchors <- FindTransferAnchors(reference = reference, query = query,dims = 1:20, reference.reduction = "pca")
predictions <- TransferData(anchorset = query.anchors, refdata = reference$celltype,dims = 1:20)
query <- AddMetaData(query, metadata = predictions)
# query <- MapQuery(anchorset = query.anchors, reference = reference, query = query,
#                       refdata = list(seurat_clusters = "seurat_clusters"), reference.reduction = "pca", reduction.model = "umap")
# dir.create(paste(savedir,"saveRDS_obj",sep = ""),showWarnings = FALSE)
# saveRDS(query, paste(savedir,"saveRDS_obj/",tetramer_name[i],"_50k_bulk_reference_mapped_obj.RDS",sep = ""))

p1 <- DimPlot(query, reduction = "tsne", group.by = "predicted.id", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(query, reduction = "umap", group.by = "predicted.id", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Reference annotation UMAP")

pdf(paste(savedir,"UMAP/Cross_sectional_reference_mapping_long.pdf",sep = ""),width = 10, height = 6)
p1 + p2
dev.off()

### Based on the discussion with the Ines
## Selected the clusters for TEMRA 16 and 7 Helios +ve, Cluster 11 and 4 Helios -ve
## Selected the clusters for EM 9 and 13 Helios +ve, Cluster 1 Helios -ve

dir.create(paste(savedir,"UMAP/celltypes/",sep = ""), showWarnings = FALSE)
cross_celltypes <- unique(query@meta.data$predicted.id)
query@meta.data$predicted.id <- factor(query@meta.data$predicted.id, levels = cross_celltypes)
for (i in 1:length(cross_celltypes)) {
  col <- rep("white",length(cross_celltypes))
  col[i] <- "red"
  p <- DimPlot(query, group.by = "predicted.id", label = TRUE, cols=alpha(col,0.3), reduction = "umap")
  pdf(paste(savedir,"UMAP/celltypes/",cross_celltypes[i],"_UMAP",".pdf",sep = ""))
  print(p)
  dev.off()

  p <- DimPlot(query, group.by = "predicted.id", label = TRUE, cols=alpha(col,0.3), reduction = "tsne")
  pdf(paste(savedir,"UMAP/celltypes/",cross_celltypes[i],"_tSNE",".pdf",sep = ""))
  print(p)
  dev.off()
}

dir.create(paste(savedir,"UMAP/clusters/",sep = ""), showWarnings = FALSE)
clusters = 0:20
query@meta.data$seurat_clusters <- factor(query@meta.data$seurat_clusters, levels = 0:20)
for (i in 1:length(clusters)) {
  col <- rep("white",length(clusters))
  col[i] <- "red"
  p <- DimPlot(query, group.by = "seurat_clusters", label = TRUE, cols=alpha(col,0.3), reduction = "umap")
  pdf(paste(savedir,"UMAP/clusters/",clusters[i],"_UMAP",".pdf",sep = ""))
  print(p)
  dev.off()

  p <- DimPlot(query, group.by = "seurat_clusters", label = TRUE, cols=alpha(col,0.3), reduction = "tsne")
  pdf(paste(savedir,"UMAP/clusters/",clusters[i],"_tSNE",".pdf",sep = ""))
  print(p)
  dev.off()
}

DefaultAssay(query) <- "RNA"
Temra_helios_pos_vs_neg <- FindMarkers(query, ident.1 = c(7,16), ident.2=c(11,4))
Temra_helios_pos_7_16_vs_neg_8 <- FindMarkers(query, ident.1 = c(7,16), ident.2=c(8))
EM_helios_pos_vs_neg <- FindMarkers(query, ident.1 = c(9,13), ident.2=c(1))

dir.create(paste(savedir,"differential/",sep = ""), showWarnings = FALSE)
write.table(Temra_helios_pos_vs_neg, paste(savedir,"differential/Temra_helios_pos_vs_neg.txt",sep = ""),
            quote = FALSE, sep = "\t", row.names = T, col.names = T)

write.table(EM_helios_pos_vs_neg, paste(savedir,"differential/EM_helios_pos_vs_neg.txt",sep = ""),
            quote = FALSE, sep = "\t", row.names = T, col.names = T)

write.table(Temra_helios_pos_7_16_vs_neg_8, paste(savedir,"differential/Temra_helios_pos_7_16_vs_neg_8.txt",sep = ""),
            quote = FALSE, sep = "\t", row.names = T, col.names = T)

dir.create(paste(savedir,"saveRDS_obj/",sep = ""))

#### pyscenic ######
####### Gene Regulatory Network #######
# https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/scenic.html
### Running SCENIC in R
library(Seurat)
library(tidyverse)
library(magrittr)
library(SCENIC)
library(SingleCellExperiment)
library(SCopeLoomR)

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/"
cross_seu_impute_ct <- query
DefaultAssay(cross_seu_impute_ct) <- 'RNA'

exprMat <- cross_seu_impute_ct@assays$RNA@data
cellInfo <- cross_seu_impute_ct@meta.data
# loci1 <- which(rowSums(exprMat) > 1*.01*ncol(exprMat))
# exprMat_filter <- exprMat[loci1, ]

add_cell_annotation <- function(loom, cellAnnotation)
{
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }

  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")

  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }

  invisible(loom)
}

dir.create(paste(savedir,"GRN",sep = ""),showWarnings = FALSE)
loom <- build_loom(paste(savedir,"GRN/cross_seu_impute_ct_normal54934.pts-96.mforgehn1ized_data.loom",sep=""), dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

## In Shell
# source ~/.bashrc
# conda activate pyscenic
# pyscenic grn /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/GRN/cross_seu_impute_ct_normalized_data.loom  /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/SCENIC/TFs/hs_hgnc_tfs.txt -m grnboost2 --seed 2 --num_workers 20 --output /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/GRN/gene_TFs_GRNBoost2.tsv
# qsub -l h_vmem=32G -pe threaded 20 -N cross_seu_pyscenic -q 1-day -o cross_GRNboost2.log -e cross_GRNboost2.err -m ae -M jain.abhinav@mayo.edu -cwd cross_GRNboost2.sh
## Nan-ping has done alignment on hg38
## pyscenic ctx gene_TFs_GRNBoost2.tsv /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/SCENIC/ranking_database/cisTarget/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/SCENIC/ranking_database/cisTarget/motif2tf_annotation/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname cross_seu_impute_ct_normalized_data.loom --mode "dask_multiprocessing" --output cross_seu_impute_ct_ctx.csv --num_workers 20 --mask_dropouts
# qsub -l h_vmem=32G -pe threaded 20 -N cross_seu_pyscenic -q 1-day -o cross_ctx.log -e cross_ctx.err -m ae -M jain.abhinav@mayo.edu -cwd pyscenic_CTX.sh
## pyscenic aucell cross_seu_impute_ct_normalized_data.loom cross_seu_impute_ct_ctx.csv --output cross_seu_impute_ct_normalized_data_scenic.loom --num_workers 20
# qsub -l h_vmem=32G -pe threaded 20 -N cross_seu_pyscenic -q 1-day -o cross_aucell.log -e cross_aucell.err -m ae -M jain.abhinav@mayo.edu -cwd aucell.sh

loom <- open_loom(paste(savedir,"GRN/cross_seu_impute_ct_normalized_data_scenic.loom",sep=""))
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
saveRDS(list('regulons' = regulons, 'regulonAUC' = regulonAUC),
        file = paste(savedir,"GRN/cross_seu_impute_ct_normalized_data_scenic_aucell_output.rds",sep=""))

AUC_mat <- as.data.frame(regulonAUC@assays@data$AUC)
write.table(AUC_mat, file = paste(savedir,"GRN/cross_seu_impute_ct_AUC_mat.csv",sep=""), quote = FALSE, sep = ",")

### Binarize the regulon AUC matrix
## Difficult to do in the R so we will do it in the Python
# conda activate pyscenic
# python
# from pyscenic.binarization import binarize
# import pandas as pd
# maindir="/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/GRN/"
# aux_mtx = pd.read_csv(maindir+"cross_seu_impute_ct_AUC_mat.csv",index_col=[0])
# auc_mtx_binary = binarize(aux_mtx)
# auc_mtx_binary_2 = pd.DataFrame(auc_mtx_binary[0])
# auc_mtx_binary_2.to_csv(maindir+'cross_seu_impute_ct_auc_mtx_binary.csv',header=True,index=True)

#### SCENIC Differential Regulons #####
library(Seurat)
library(tidyverse)
library(magrittr)
library(viridis)
library(SCENIC)
library(network)
library(igraph)
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/"
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/GRN/"
cross_seu_impute_ct = readRDS(paste(maindir,"saveRDS_obj/cross_sectional_impute.RDS",sep = ""))
AUC_mat <- read.csv(paste(savedir,"cross_seu_impute_ct_AUC_mat.csv",sep = ""), header = TRUE)
rownames(AUC_mat) <- gsub("[(+)]", "", rownames(AUC_mat))
# AUC_binary <- read.csv(paste(savedir,"cross_seu_impute_ct/cross_seu_impute_ct_auc_mtx_binary.csv",sep = ""), header = TRUE)
# rownames(AUC_binary) <- AUC_binary$X
# AUC_binary$X <- NULL
# rownames(AUC_binary) <- gsub("[(+)]", "", rownames(AUC_binary))
# colnames(AUC_binary) <- gsub(".1$","-1",colnames(AUC_binary))
# all(colnames(AUC_binary) == colnames(cross_seu_impute_ct@assays$RNA@counts))
colnames(AUC_mat) <- gsub("*\\.","-",colnames(AUC_mat))
all(colnames(AUC_mat) == colnames(cross_seu_impute_ct@assays$RNA@counts))

cross_seu_impute_ct[['AUC']] <- CreateAssayObject(data = as.matrix(AUC_mat))
# cross_seu_impute_ct[['AUCBinary']] <- CreateAssayObject(data = as.matrix(AUC_binary))

# Identify differential regulons based AUC scores
DefaultAssay(cross_seu_impute_ct) <- "AUC"
cross_seu_impute_ct_markers <- FindAllMarkers(cross_seu_impute_ct, logfc.threshold = 0.02)
cross_seu_impute_ct_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top10

write.table(cross_seu_impute_ct_markers, paste(savedir,"cross_seu_impute_ct_cluster_markers.txt",sep = ""), row.names = TRUE,
            col.names = TRUE, sep = "\t", quote = FALSE)

DefaultAssay(cross_seu_impute_ct) <- 'AUC'
cross_seu_impute_ct <- ScaleData(cross_seu_impute_ct, assay = 'AUC')

pdf(paste(maindir,"heatmap/cross_seu_impute_ct_heatmap_TF_top10_all.pdf",sep = ""), width = 10, height = 10)
print(DoHeatmap(cross_seu_impute_ct, features = top10$gene) + NoLegend())
dev.off()
cross_seu_impute_ct <- ScaleData(cross_seu_impute_ct, assay = 'RNA')

cross_seu_impute_ct_markers <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run02_03_04_age_65_80_analysis/GRN/cross_seu_impute_ct/cross_seu_impute_ct_cluster_markers.txt")
clus <- unique(cross_seu_impute_ct_markers$cluster)

for (i in 1:length(clus)) {
  features <- cross_seu_impute_ct_markers[grep(paste("^",clus[i],"$",sep = ""),cross_seu_impute_ct_markers$cluster),"gene"]
  DefaultAssay(cross_seu_impute_ct) <- "AUC"
  pdf(paste(maindir,"heatmap/cross_seu_impute_ct_heatmap_TF_",clus[i],"_pyscenic_score.pdf",sep = ""), width = 8, height = 8)
  print(DoHeatmap(cross_seu_impute_ct, features = features) + NoLegend())
  dev.off()
}

saveRDS(cross_seu_impute_ct,file = paste(maindir,"saveRDS_obj/cross_seu_impute_ct_modality_integrate_cluster_0.8_pyscenic.RDS",sep=""))

cross_seu_impute_ct <- readRDS(paste(maindir,"saveRDS_obj/cross_seu_impute_ct_modality_integrate_cluster_0.8_pyscenic.RDS",sep=""))

#### Extracting TEMRA #####
## Cluster 11,4,7,16,8
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/Temra/"
cross_seu_impute_ct = readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/saveRDS_obj/cross_sectional_impute.RDS")
DefaultAssay(cross_seu_impute_ct) = "RNA"
temra_seu <- subset(cross_seu_impute_ct, idents=c(4,7,8,11,16))

temra_seu <- NormalizeData(temra_seu, normalization.method = "LogNormalize", scale.factor = 10000)
temra_seu <- FindVariableFeatures(temra_seu, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(temra_seu)
temra_seu <- ScaleData(temra_seu, features = all.genes)
temra_seu <- RunPCA(temra_seu, features = VariableFeatures(object = temra_seu))
pdf(paste(savedir,"Elbow_plot.pdf",sep = ""))
ElbowPlot(temra_seu, ndims = 50)
dev.off()

dims <- c(1:20)
seed <- c(7456)
temra_seu <- RunTSNE(object = temra_seu, dims = dims,seed.use = seed)

temra_seu <- RunUMAP(temra_seu, dims = 1:20)
temra_seu <- FindNeighbors(temra_seu, dims = 1:20)

DefaultAssay(temra_seu) <- "RNA"
gene_zero <- as.vector(which(rowSums(temra_seu@assays$RNA@data)==0))
temra_seu_2 <- subset(temra_seu, features = rownames(temra_seu[-gene_zero,]))
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/python/python-3.9.2/bin/python3.9")
library(Rmagic)
py_discover_config("magic") # to check
temra_seu_impute <- magic(temra_seu_2, npca=20) ## imputing the vasc data
DefaultAssay(temra_seu_impute) <- "MAGIC_RNA"

dir.create(paste(savedir,"featureplot",sep = ""),showWarnings = FALSE)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(temra_seu_impute) = "MAGIC_RNA"
p1 <- featureplot_front(temra_seu_impute,"IKZF2", reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits = c(0,0.3))
p2 <- featureplot_front(temra_seu_impute,"IKZF2", reduction = "umap", x="UMAP_1", y="UMAP_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits = c(0,0.3))
pdf(paste(savedir,"featureplot/IKZF2_tSNE_UMAP.pdf",sep = ""), width = 9, height = 5)
p1+p2
dev.off()

dir.create(paste(savedir,"UMAP",sep = ""),showWarnings = FALSE)
dir.create(paste(savedir,"Vlnplot",sep = ""),showWarnings = FALSE)

res = c(0.2,0.4,0.6,0.8,1,1.2)
res = 0.6
for (i in 1:length(res)) {
  DefaultAssay(temra_seu_impute) <- "RNA"
  temra_seu_impute <- FindClusters(temra_seu_impute, resolution = res[i])
  pdf(paste(savedir,"UMAP/dimplot_tsne_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(temra_seu_impute, reduction = 'tsne', label = TRUE, label.size = 5))
  dev.off()
  
  pdf(paste(savedir,"UMAP/dimplot_umap_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(temra_seu_impute, reduction = 'umap', label = TRUE, label.size = 5))
  dev.off()
  
  DefaultAssay(temra_seu_impute) <- "RNA"
  pdf(paste(savedir,"Vlnplot/IKZF2_boxplot_",res[i],".pdf",sep = ""),width = 10, height = 6)
  print(VlnPlot(temra_seu_impute, features = "IKZF2",pt.size = 0) + geom_boxplot())
  dev.off()
  
  DefaultAssay(temra_seu_impute) <- "MAGIC_RNA"
  pdf(paste(savedir,"Vlnplot/IKZF2_MAGIC_boxplot_",res[i],".pdf",sep = ""),width = 10, height = 6)
  print(VlnPlot(temra_seu_impute, features = "IKZF2",pt.size = 0) + geom_boxplot())
  dev.off()
}

DefaultAssay(temra_seu_impute) <- "RNA"
Helios_positive = FindMarkers(temra_seu_impute, ident.1 = c(0,8,9), min.pct = 0.05)
temra_seu_impute = temra 
n_cells_2 <- FetchData(temra_seu_impute,vars = c("ident", "orig.ident"))
n_cells <- table(n_cells_2[,"orig.ident"],n_cells_2[,"ident"])
n_cells_sum <- as.vector(rowSums(n_cells))

dir.create(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/Temra/Table/",sep = ""),showWarnings = F)
write.table(n_cells,
            paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/Temra/Table/TEMRA_cluster_and_orig.ident.txt",sep = ""),
            quote = F,
            row.names = T,
            col.names = T)


write.table(Helios_positive, paste(savedir,"Helios_positive_TEMRA_0_8_9_cluster.txt",sep = ""), sep = "\t", quote = FALSE)
dir.create(paste(savedir,"saveRDS_obj/",sep = ""), showWarnings = FALSE)
saveRDS(temra_seu_impute,paste(savedir,"saveRDS_obj/temra_sectional_impute_res_0.6.RDS",sep = ""))

temra_list <- list()
a=1
for (i in 1:length(temra_scenic$regulons)) {
  p <- grep("IKZF2",temra_scenic$regulons[[i]])
  if(length(p)>0){
    temra_list[[a]]<-i
    a=a+1
  }
}

TF_IKZF2_TEMRA <- names(temra_scenic$regulons)[unlist(temra_list)]
TF_IKZF2_TEMRA <- gsub("[(+)]", "", TF_IKZF2_TEMRA)
DefaultAssay(temra_seu) <- "AUC"
FeaturePlot(temra_seu,TF_IKZF2_TEMRA, reduction = "tsne")

## Find Helios IKZF2 


### Extracting scale data
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/heatmap_scRNA.R")
AJ_heatmap_RNA(obj = temra_seu_impute, savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/Temra/",
                  factorize = c(0:9), objname = "TEMRA_cluster",sel_col = "seurat_clusters")
AJ_heatmap_RNA_markers <- FindMarkers(temra_seu, ident.1 = c(0,8,9), ident.2 = c(1,2,4,5,6,7), logfc.threshold = 0.01,min.pct = 0.01)
AJ_heatmap_RNA_markers <- FindAllMarkers(temra_seu, logfc.threshold = 0.01,min.pct = 0.01)
genes <- unique(AJ_heatmap_RNA_markers$gene)
DoHeatmap(temra_seu, features = genes)
### GRN for TEMRA Helios ####
# https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/scenic.html
### Running SCENIC in R
library(Seurat)
library(tidyverse)
library(magrittr)
library(SCENIC)
library(SingleCellExperiment)
library(SCopeLoomR)

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/Temra/"
temra_seu_impute = readRDS(paste(savedir,"saveRDS_obj/temra_sectional_impute_res_0.6.RDS",sep = ""))
DefaultAssay(temra_seu_impute) <- 'RNA'

exprMat <- temra_seu_impute@assays$RNA@data
cellInfo <- temra_seu_impute@meta.data
# loci1 <- which(rowSums(exprMat) > 1*.01*ncol(exprMat))
# exprMat_filter <- exprMat[loci1, ]

add_cell_annotation <- function(loom, cellAnnotation)
{
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}

dir.create(paste(savedir,"GRN",sep = ""),showWarnings = FALSE)
loom <- build_loom(paste(savedir,"GRN/temra_seu_impute_normalized_data.loom",sep=""), dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

## In Shell
# source ~/.bashrc
# conda activate pyscenic
# pyscenic grn /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/GRN/temra_seu_impute_normalized_data.loom  /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/SCENIC/TFs/hs_hgnc_tfs.txt -m grnboost2 --seed 2 --num_workers 20 --output /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/GRN/gene_TFs_GRNBoost2.tsv
# qsub -l h_vmem=32G -pe threaded 20 -N temra_seu_pyscenic -q 1-day -o temra_GRNboost2.log -e temra_GRNboost2.err -m ae -M jain.abhinav@mayo.edu -cwd temra_GRNboost2.sh
## Nan-ping has done alignment on hg38
## pyscenic ctx gene_TFs_GRNBoost2.tsv /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/SCENIC/ranking_database/cisTarget/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/SCENIC/ranking_database/cisTarget/motif2tf_annotation/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname temra_seu_impute_normalized_data.loom --mode "dask_multiprocessing" --output temra_seu_impute_ctx.csv --num_workers 20 --mask_dropouts
# qsub -l h_vmem=32G -pe threaded 20 -N temra_seu_pyscenic -q 1-day -o temra_ctx.log -e temra_ctx.err -m ae -M jain.abhinav@mayo.edu -cwd pyscenic_CTX.sh
# pyscenic aucell temra_seu_impute_normalized_data.loom temra_seu_impute_ctx.csv --output temra_seu_impute_normalized_data_scenic.loom --num_workers 20
# qsub -l h_vmem=32G -pe threaded 20 -N temra_seu_pyscenic -q 1-day -o temra_aucell.log -e temra_aucell.err -m ae -M jain.abhinav@mayo.edu -cwd aucell.sh

loom <- open_loom(paste(savedir,"GRN/temra_seu_impute_normalized_data_scenic.loom",sep=""))
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
saveRDS(list('regulons' = regulons, 'regulonAUC' = regulonAUC),
        file = paste(savedir,"GRN/temra_seu_impute_normalized_data_scenic_aucell_output.rds",sep=""))

AUC_mat <- as.data.frame(regulonAUC@assays@data$AUC)
write.table(AUC_mat, file = paste(savedir,"GRN/temra_seu_impute_ct_AUC_mat.csv",sep=""), quote = FALSE, sep = ",")

### Binarize the regulon AUC matrix
## Difficult to do in the R so we will do it in the Python
# conda activate pyscenic
# python
# from pyscenic.binarization import binarize
# import pandas as pd
# maindir="/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/GRN/"
# aux_mtx = pd.read_csv(maindir+"temra_seu_impute_ct_AUC_mat.csv",index_col=[0])
# auc_mtx_binary = binarize(aux_mtx)
# auc_mtx_binary_2 = pd.DataFrame(auc_mtx_binary[0])
# auc_mtx_binary_2.to_csv(maindir+'temra_seu_impute_ct_auc_mtx_binary.csv',header=True,index=True)

#### SCENIC Differential Regulons #####
library(Seurat)
library(tidyverse)
library(magrittr)
library(viridis)
library(SCENIC)
library(network)
library(igraph)
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/Temra/"
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/Temra/GRN/"
temra_seu_impute = readRDS(paste(maindir,"saveRDS_obj/temra_sectional_impute_res_0.6.RDS",sep = ""))
AUC_mat <- read.csv(paste(savedir,"temra_seu_impute_ct_AUC_mat.csv",sep = ""), header = TRUE)
rownames(AUC_mat) <- gsub("[(+)]", "", rownames(AUC_mat))
# AUC_binary <- read.csv(paste(savedir,"temra_seu_impute_ct/temra_seu_impute_ct_auc_mtx_binary.csv",sep = ""), header = TRUE)
# rownames(AUC_binary) <- AUC_binary$X
# AUC_binary$X <- NULL
# rownames(AUC_binary) <- gsub("[(+)]", "", rownames(AUC_binary))
# colnames(AUC_binary) <- gsub(".1$","-1",colnames(AUC_binary))
# all(colnames(AUC_binary) == colnames(temra_seu_impute_ct@assays$RNA@counts))
colnames(AUC_mat) <- gsub("*\\.","-",colnames(AUC_mat))
all(colnames(AUC_mat) == colnames(temra_seu_impute@assays$RNA@counts))

temra_seu_impute[['AUC']] <- CreateAssayObject(data = as.matrix(AUC_mat))
# temra_seu_impute_ct[['AUCBinary']] <- CreateAssayObject(data = as.matrix(AUC_binary))

# Identify differential regulons based AUC scores
DefaultAssay(temra_seu_impute) <- "AUC"
temra_seu_impute_markers <- FindAllMarkers(temra_seu_impute, logfc.threshold = 0.01)
temra_seu_impute_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top10

write.table(temra_seu_impute_markers, paste(savedir,"temra_seu_impute_cluster_markers.txt",sep = ""), row.names = TRUE,
            col.names = TRUE, sep = "\t", quote = FALSE)

DefaultAssay(temra_seu_impute) <- 'AUC'
temra_seu_impute <- ScaleData(temra_seu_impute, assay = 'AUC')

pdf(paste(savedir,"temra_seu_impute_heatmap_TF_top10_all.pdf",sep = ""), height = 5, width = 7)
print(DoHeatmap(temra_seu_impute, features = top10$gene) + NoLegend())
dev.off()

temra_seu_impute <- ScaleData(temra_seu_impute, assay = 'RNA')

# temra_seu_impute_markers <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run02_03_04_age_65_80_analysis/GRN/cross_seu_impute_ct/cross_seu_impute_ct_cluster_markers.txt")
# clus <- unique(temra_seu_impute_markers$cluster)
# 
# for (i in 1:length(clus)) {
#   features <- temra_seu_impute_markers[grep(paste("^",clus[i],"$",sep = ""),temra_seu_impute_markers$cluster),"gene"]
#   DefaultAssay(temra_seu_impute) <- "AUC"
#   pdf(paste(maindir,"heatmap/temra_seu_impute_heatmap_TF_",clus[i],"_pyscenic_score.pdf",sep = ""), width = 8, height = 8)
#   print(DoHeatmap(temra_seu_impute, features = features) + NoLegend())
#   dev.off()
# }

saveRDS(temra_seu_impute,file = paste(savedir,"saveRDS_obj/temra_seu_impute_0.6_pyscenic.RDS",sep=""))
temra_seu_impute <- readRDS(paste(maindir,"saveRDS_obj/temra_seu_impute_0.6_pyscenic.RDS",sep=""))

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/Temra/"
temra_seu <- readRDS(paste(savedir,"saveRDS_obj/temra_seu_impute_0.6_pyscenic.RDS",sep = ""))


### making a Network
### Making for the SCope using R
# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# For some of the plots:
#library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
packageVersion("SCENIC")
##  sed 's/\t/,/g' gene_TFs_GRNBoost2.tsv > gene_TFs_GRNBoost2.csv in Bash

visuNetwork <- function(regulon.name,genes){
  adj.ls <- regulon.name %>% purrr::map(~{
    tmp <- adj[which(adj$TF == .x), ]
    tmp <- tmp[order(tmp$importance, decreasing = T),]
    # loci <- unique(c(1:50, grep(regulon.name, tmp$target)))
    loci <- which(genes %in% tmp$target)
    tmp <- tmp[loci,]
    return(tmp)
  })
  
  adj.sub <- Reduce('rbind', adj.ls)
  
  ## generate network
  edge.df <- adj.sub
  colnames(edge.df) <- c('from', 'to', 'weights')
  edge.df$from <- as.character(edge.df$from)
  edge.df$to <- as.character(edge.df$to)
  
  #net1 <- network(edge.df, directed = T, loops = T)
  #plot(net1)
  
  vertex <- unique(c(edge.df$from, edge.df$to))
  net <- graph_from_data_frame(d = edge.df, vertices = vertex, directed = T)
  
  
  ## vertex size scaled to log2FC
  #combined.tmp <- AverageExpression(combined, assays = 'RNA', features = vertex)
  #fc <- combined.tmp$RNA$SC_Hyper/combined.tmp$RNA$SC_AS
  #fc <- combined.tmp$RNA$SC_AS/combined.tmp$RNA$SC_Hyper
  #fc[which(fc > quantile(fc, 0.9))] <- quantile(fc, 0.9)
  #fc[which(fc < quantile(fc, 0.1))] <- quantile(fc, 0.1)
  #vsize <- scale(fc, center = min(fc), scale = max(fc) - min(fc))
  #vsize <- vsize+0.1
  
  ## vertex color: whether DEG in hyper vs as
  # vcl <- ifelse(vertex %in% dn.deg, 'blue', ifelse(vertex %in% up.deg, 'red', 'grey'))
  #  vlabel <- rep("", length(vertex))
  # vertex label
  # loci <- which(vcl != 'grey' | vertex %in% edge.df$from | vertex == 'MUC2')
  #  vlabel[loci] = vertex[loci]
  # vertex size scaled to weights
  vsize <- scale(edge.df$weights, center = min(quantile(edge.df$weights)), scale = max(quantile(edge.df$weights)) - min(quantile(edge.df$weights)))
  
  
  plot(
    net,
    # vertex.label = vlabel,
    vertex.label.cex = 1,
    vertex.label.dist = 1,
    vertex.label.color = 'black',
    vertex.size = c(1, vsize+0.1)*10,
    #vertex.size = vsize*10,
    vertex.frame.color = NA,
    edge.arrow.size = 0.5,
    #  vertex.color = vcl,
    main = regulon.name
  )
}

library(purrr)
library(igraph)
adj <- read.csv(paste(savedir,"GRN/gene_TFs_GRNBoost2.csv",sep = ""))
regulon_list = readRDS(paste(savedir,"GRN/temra_seu_impute_normalized_data_scenic_aucell_output.rds",sep=""))
regulons <- regulon_list$regulons
names(regulons) <- gsub("[(+)]", "", names(regulons))
visuNetwork("IKZF2",regulons$IKZF2)
regulons$EOMES


### Running EM ######
## Cluster 11,4,7,16,8
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/EM/"
cross_seu_impute_ct = readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/saveRDS_obj/cross_sectional_impute.RDS")
DefaultAssay(cross_seu_impute_ct) = "RNA"
dir.create(paste(savedir,"UMAP/",sep = ""))
res = c(1.4,1.6,1.8,2)
res = c(1.2)
for (i in 1:length(res)) {
  DefaultAssay(cross_seu_impute_ct) <- "RNA"
  cross_seu_impute_ct <- FindClusters(cross_seu_impute_ct, resolution = res[i])
  pdf(paste(savedir,"UMAP/dimplot_tsne_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(cross_seu_impute_ct, reduction = 'tsne', label = TRUE, label.size = 5))
  dev.off()
  
  pdf(paste(savedir,"UMAP/dimplot_umap_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(cross_seu_impute_ct, reduction = 'umap', label = TRUE, label.size = 5))
  dev.off()
}


EM <- subset(cross_seu_impute_ct, idents=c(1,9,13))
EM <- NormalizeData(EM, normalization.method = "LogNormalize", scale.factor = 10000)
EM <- FindVariableFeatures(EM, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(EM)
EM <- ScaleData(EM, features = all.genes)
EM <- RunPCA(EM, features = VariableFeatures(object = EM))
pdf(paste(savedir,"Elbow_plot.pdf",sep = ""))
ElbowPlot(EM, ndims = 50)
dev.off()

dims <- c(1:20)
seed <- c(7456)
EM <- RunTSNE(object = EM, dims = dims,seed.use = seed)

EM <- RunUMAP(EM, dims = 1:20)
EM <- FindNeighbors(EM, dims = 1:20)

DefaultAssay(EM) <- "RNA"
gene_zero <- as.vector(which(rowSums(EM@assays$RNA@data)==0))
EM_2 <- subset(EM, features = rownames(EM[-gene_zero,]))
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/python/python-3.9.2/bin/python3.9")
library(Rmagic)
py_discover_config("magic") # to check
EM_impute <- magic(EM_2, npca=20) ## imputing the vasc data
DefaultAssay(EM_impute) <- "MAGIC_RNA"
dir.create(paste(savedir,"saveRDS_obj/",sep = ""))
saveRDS(EM_impute,paste(savedir,"saveRDS_obj/EM_sectional_impute.RDS",sep = ""))

dir.create(paste(savedir,"featureplot",sep = ""),showWarnings = FALSE)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Shiny/multiple_dataset/Data/express_cell_front.R")
DefaultAssay(EM_impute) = "MAGIC_RNA"
p1 <- featureplot_front(EM_impute,"IKZF2", reduction = "tsne", x="tSNE_1", y="tSNE_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits = c(0,0.3))
p2 <- featureplot_front(EM_impute,"IKZF2", reduction = "umap", x="UMAP_1", y="UMAP_2",size=0.01) +
  scale_color_gradientn(colours = ArchRPalettes$solarExtra, limits = c(0,0.3))
pdf(paste(savedir,"featureplot/IKZF2_tSNE_UMAP.pdf",sep = ""), width = 9, height = 5)
p1+p2
dev.off()

dir.create(paste(savedir,"UMAP",sep = ""),showWarnings = FALSE)
dir.create(paste(savedir,"Vlnplot",sep = ""),showWarnings = FALSE)

EM_impute = readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/EM/saveRDS_obj/EM_cross_impute.RDS")
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/EM/"
res = c(0.2,0.4,0.6,0.8,1,1.2)
res = 0.6
for (i in 1:length(res)) {
  DefaultAssay(EM_impute) <- "RNA"
  EM_impute <- FindClusters(EM_impute, resolution = res[i])
  pdf(paste(savedir,"UMAP/dimplot_tsne_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(EM_impute, reduction = 'tsne', label = TRUE, label.size = 5))
  dev.off()
  
  pdf(paste(savedir,"UMAP/dimplot_umap_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(EM_impute, reduction = 'umap', label = TRUE, label.size = 5))
  dev.off()
  
  DefaultAssay(EM_impute) <- "RNA"
  pdf(paste(savedir,"Vlnplot/IKZF2_boxplot_",res[i],".pdf",sep = ""),width = 10, height = 6)
  print(VlnPlot(EM_impute, features = "IKZF2",pt.size = 0) + geom_boxplot())
  dev.off()
  
  DefaultAssay(EM_impute) <- "MAGIC_RNA"
  pdf(paste(savedir,"Vlnplot/IKZF2_MAGIC_boxplot_",res[i],".pdf",sep = ""),width = 10, height = 6)
  print(VlnPlot(EM_impute, features = "IKZF2",pt.size = 0) + geom_boxplot())
  dev.off()
}

n_cells_2 <- FetchData(EM_impute,vars = c("ident", "orig.ident"))
n_cells <- table(n_cells_2[,"orig.ident"],n_cells_2[,"ident"])
n_cells_sum <- as.vector(rowSums(n_cells))

dir.create(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/EM/Table/",sep = ""),showWarnings = F)
write.table(n_cells,
            paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/EM/Table/EM_cluster_and_orig.ident.txt",sep = ""),
            quote = F,
            row.names = T,
            col.names = T)


saveRDS(EM_impute,paste(savedir,"saveRDS_obj/EM_sectional_impute.RDS",sep = ""))

DefaultAssay(EM_impute) <- "RNA"
Helios_positive = FindMarkers(EM_impute, ident.1 = c(0,3,8,9,10,11,12), min.pct = 0.05, logfc.threshold = 0.1)

write.table(Helios_positive, paste(savedir,"Helios_positive_EM_0_3_8_9_10_11_12_cluster.txt",sep = ""), sep = "\t", quote = FALSE)

dir.create(paste(savedir,"saveRDS_obj/",sep = ""), showWarnings = FALSE)
saveRDS(EM_impute,paste(savedir,"saveRDS_obj/EM_cross_impute.RDS",sep = ""))

### GRN for EM Helios ####
# https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/scenic.html
### Running SCENIC in R
library(Seurat)
library(tidyverse)
library(magrittr)
library(SCENIC)
library(SingleCellExperiment)
library(SCopeLoomR)

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/EM/"
EM_seu_impute = readRDS(paste(savedir,"saveRDS_obj/EM_sectional_impute_res_0.6.RDS",sep = ""))
DefaultAssay(EM_seu_impute) <- 'RNA'

exprMat <- EM_seu_impute@assays$RNA@data
cellInfo <- EM_seu_impute@meta.data
# loci1 <- which(rowSums(exprMat) > 1*.01*ncol(exprMat))
# exprMat_filter <- exprMat[loci1, ]

add_cell_annotation <- function(loom, cellAnnotation)
{
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}

dir.create(paste(savedir,"GRN",sep = ""),showWarnings = FALSE)
loom <- build_loom(paste(savedir,"GRN/EM_seu_impute_normalized_data.loom",sep=""), dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

## In Shell
# source ~/.bashrc
# conda activate pyscenic
# pyscenic grn /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/GRN/EM_seu_impute_normalized_data.loom /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/SCENIC/TFs/hs_hgnc_tfs.txt -m grnboost2 --seed 2 --num_workers 20 --output /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/EM/GRN/gene_TFs_GRNBoost2.tsv
# qsub -l h_vmem=32G -pe threaded 20 -N EM_seu_pyscenic -q 1-day -o EM_GRNboost2.log -e EM_GRNboost2.err -m ae -M jain.abhinav@mayo.edu -cwd EM_GRNboost2.sh
## pyscenic ctx gene_TFs_GRNBoost2.tsv /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/SCENIC/ranking_database/cisTarget/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/SCENIC/ranking_database/cisTarget/motif2tf_annotation/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname EM_seu_impute_normalized_data.loom --mode "dask_multiprocessing" --output EM_seu_impute_ctx.csv --num_workers 20 --mask_dropouts
# qsub -l h_vmem=32G -pe threaded 20 -N EM_seu_pyscenic -q 1-day -o EM_ctx.log -e EM_ctx.err -m ae -M jain.abhinav@mayo.edu -cwd pyscenic_CTX.sh
# pyscenic aucell EM_seu_impute_normalized_data.loom EM_seu_impute_ctx.csv --output EM_seu_impute_normalized_data_scenic.loom --num_workers 20
# qsub -l h_vmem=32G -pe threaded 20 -N EM_seu_pyscenic -q 1-day -o EM_aucell.log -e EM_aucell.err -m ae -M jain.abhinav@mayo.edu -cwd aucell.sh

loom <- open_loom(paste(savedir,"GRN/EM_seu_impute_normalized_data_scenic.loom",sep=""))
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
saveRDS(list('regulons' = regulons, 'regulonAUC' = regulonAUC),
        file = paste(savedir,"GRN/EM_seu_impute_normalized_data_scenic_aucell_output.rds",sep=""))

AUC_mat <- as.data.frame(regulonAUC@assays@data$AUC)
write.table(AUC_mat, file = paste(savedir,"GRN/EM_seu_impute_AUC_mat.csv",sep=""), quote = FALSE, sep = ",")

#### SCENIC Differential Regulons #####
library(Seurat)
library(tidyverse)
library(magrittr)
library(viridis)
library(SCENIC)
library(network)
library(igraph)
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/EM/"
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/EM/GRN/"
temra_seu_impute = readRDS(paste(maindir,"saveRDS_obj/EM_sectional_impute_res_0.6.RDS",sep = ""))
AUC_mat <- read.csv(paste(savedir,"temra_seu_impute_ct_AUC_mat.csv",sep = ""), header = TRUE)
rownames(AUC_mat) <- gsub("[(+)]", "", rownames(AUC_mat))
# AUC_binary <- read.csv(paste(savedir,"temra_seu_impute_ct/temra_seu_impute_ct_auc_mtx_binary.csv",sep = ""), header = TRUE)
# rownames(AUC_binary) <- AUC_binary$X
# AUC_binary$X <- NULL
# rownames(AUC_binary) <- gsub("[(+)]", "", rownames(AUC_binary))
# colnames(AUC_binary) <- gsub(".1$","-1",colnames(AUC_binary))
# all(colnames(AUC_binary) == colnames(temra_seu_impute_ct@assays$RNA@counts))
colnames(AUC_mat) <- gsub("*\\.","-",colnames(AUC_mat))
all(colnames(AUC_mat) == colnames(EM_seu_impute@assays$RNA@counts))

EM_seu_impute[['AUC']] <- CreateAssayObject(data = as.matrix(AUC_mat))
# EM_seu_impute_ct[['AUCBinary']] <- CreateAssayObject(data = as.matrix(AUC_binary))

# Identify differential regulons based AUC scores
DefaultAssay(EM_seu_impute) <- "AUC"
EM_seu_impute_markers <- FindAllMarkers(EM_seu_impute, logfc.threshold = 0.01)
EM_seu_impute_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top10

write.table(EM_seu_impute_markers, paste(savedir,"EM_seu_impute_cluster_markers.txt",sep = ""), row.names = TRUE,
            col.names = TRUE, sep = "\t", quote = FALSE)

DefaultAssay(EM_seu_impute) <- 'AUC'
EM_seu_impute <- ScaleData(EM_seu_impute, assay = 'AUC')

pdf(paste(savedir,"EM_seu_impute_heatmap_TF_top10_all.pdf",sep = ""), height = 5, width = 7)
print(DoHeatmap(EM_seu_impute, features = top10$gene) + NoLegend())
dev.off()

EM_seu_impute <- ScaleData(EM_seu_impute, assay = 'RNA')

# EM_seu_impute_markers <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run02_03_04_age_65_80_analysis/GRN/cross_seu_impute_ct/cross_seu_impute_ct_cluster_markers.txt")
# clus <- unique(EM_seu_impute_markers$cluster)
# 
# for (i in 1:length(clus)) {
#   features <- EM_seu_impute_markers[grep(paste("^",clus[i],"$",sep = ""),EM_seu_impute_markers$cluster),"gene"]
#   DefaultAssay(EM_seu_impute) <- "AUC"
#   pdf(paste(maindir,"heatmap/EM_seu_impute_heatmap_TF_",clus[i],"_pyscenic_score.pdf",sep = ""), width = 8, height = 8)
#   print(DoHeatmap(EM_seu_impute, features = features) + NoLegend())
#   dev.off()
# }

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/EM/"
saveRDS(EM_seu_impute,file = paste(savedir,"saveRDS_obj/EM_seu_impute_0.6_pyscenic.RDS",sep=""))
EM_seu_impute <- readRDS(paste(maindir,"saveRDS_obj/EM_seu_impute_0.6_pyscenic.RDS",sep=""))

### making a Network
### Making for the SCope using R
# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# For some of the plots:
#library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
packageVersion("SCENIC")
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/EM/"
##  sed 's/\t/,/g' gene_TFs_GRNBoost2.tsv > gene_TFs_GRNBoost2.csv in Bash

visuNetwork <- function(regulon.name,genes){

  adj.ls <- regulon.name %>% purrr::map(~{
    tmp <- adj[which(adj$TF == .x), ]
    tmp <- tmp[order(tmp$importance, decreasing = T),]
    # loci <- unique(c(1:50, grep(regulon.name, tmp$target)))
    loci <- match(genes,tmp$target)
    tmp <- tmp[loci,]
    return(tmp)
  })

  adj.sub <- Reduce('rbind', adj.ls)

  ## generate network
  edge.df <- adj.sub
  colnames(edge.df) <- c('from', 'to', 'weights')
  edge.df$from <- as.character(edge.df$from)
  edge.df$to <- as.character(edge.df$to)

  #net1 <- network(edge.df, directed = T, loops = T)
  #plot(net1)

  vertex <- unique(c(edge.df$from, edge.df$to))
  net <- graph_from_data_frame(d = edge.df, vertices = vertex, directed = T)


  ## vertex size scaled to log2FC
  #combined.tmp <- AverageExpression(combined, assays = 'RNA', features = vertex)
  #fc <- combined.tmp$RNA$SC_Hyper/combined.tmp$RNA$SC_AS
  #fc <- combined.tmp$RNA$SC_AS/combined.tmp$RNA$SC_Hyper
  #fc[which(fc > quantile(fc, 0.9))] <- quantile(fc, 0.9)
  #fc[which(fc < quantile(fc, 0.1))] <- quantile(fc, 0.1)
  #vsize <- scale(fc, center = min(fc), scale = max(fc) - min(fc))
  #vsize <- vsize+0.1

  ## vertex color: whether DEG in hyper vs as
  # vcl <- ifelse(vertex %in% dn.deg, 'blue', ifelse(vertex %in% up.deg, 'red', 'grey'))
  #  vlabel <- rep("", length(vertex))
  # vertex label
  # loci <- which(vcl != 'grey' | vertex %in% edge.df$from | vertex == 'MUC2')
  #  vlabel[loci] = vertex[loci]
  # vertex size scaled to weights
  vsize <- scale(edge.df$weights, center = min(quantile(edge.df$weights)), scale = max(quantile(edge.df$weights)) - min(quantile(edge.df$weights)))


  plot(
    net,
    # vertex.label = vlabel,
    vertex.label.cex = 1,
    vertex.label.dist = 1,
    vertex.label.color = 'black',
    vertex.size = c(1, vsize+0.1)*10,
    #vertex.size = vsize*10,
    vertex.frame.color = NA,
    edge.arrow.size = 0.5,
    #  vertex.color = vcl,
    main = regulon.name
  )
}

library(purrr)
library(igraph)
# sed 's/\t/,/g' gene_TFs_GRNBoost2.tsv > gene_TFs_GRNBoost2.csv
adj <- read.csv(paste(savedir,"GRN/gene_TFs_GRNBoost2.csv",sep = ""))
em_scenic = readRDS(paste(savedir,"GRN/EM_seu_impute_normalized_data_scenic_aucell_output.rds",sep=""))
regulons <- em_scenic$regulons
names(regulons) <- gsub("[(+)]", "", names(regulons))
visuNetwork("IKZF2",regulons$IKZF2)
regulons$EOMES

em_scenic <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/subclustering/EM/GRN/EM_seu_impute_normalized_data_scenic_aucell_output.rds")
em_list <- list()
a=1
for (i in 1:length(em_scenic$regulons)) {
  p <- grep("IKZF2",em_scenic$regulons[[i]])
  if(length(p)>0){
    em_list[[a]]<-i
    a=a+1
  }
}

TF_IKZF2_em <- names(em_scenic$regulons)[unlist(em_list)]
TF_IKZF2_em <- gsub("[(+)]", "", TF_IKZF2_em)
DefaultAssay(EM_seu) <- "MAGIC_RNA"
FeaturePlot(EM_seu,TF_IKZF2_em, reduction = "tsne")

DefaultAssay(EM_seu) <- "AUC"
AUC_markers<- FindAllMarkers(EM_seu, logfc.threshold = 0.05)
AUC_markers_ident<- FindMarkers(EM_seu, ident.1 = c(0,1,4), logfc.threshold = 0.05)
DoHeatmap(EM_seu, features = AUC_markers$gene, assay = "AUC")
FeaturePlot(EM_seu,"FOS", reduction = "tsne") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)

#### CD8 all #######
library(Seurat)
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/all_genes/"
dir.create(savedir, showWarnings = FALSE)
rawdata = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/"
raw_counts <- readRDS(paste(rawdata,"GSE136184_rsc_Gene_Expression_Aging_UMI_Counts2.Rds",sep = ""))
seu <- readRDS(paste(rawdata,"GSE136184_seu_obj2.Rds",sep = ""))
# gene_req <- rownames(seu)
# raw_counts_gene_req <- raw_counts[match(gene_req,rownames(raw_counts)),]

metadata <- read.table(paste(rawdata,"individual_cell_counts_metadata.txt",sep = ""), header = TRUE)
metadata$study <- metadata$donor_ID
metadata[grep("sc",metadata$study),"study"] <- "cross"
metadata[grep("L",metadata$study),"study"] <- "long"

seu_NP <- seu
p1p4 <- CreateSeuratObject(counts = raw_counts, project="Nanping_Object")
all(rownames(p1p4@meta.data) == all(seu_NP@meta.data)) # if true move ahead
p1p4@meta.data$Code <- seu_NP@meta.data$Code
p1p4@meta.data$visit <- seu_NP@meta.data$visit
p1p4@meta.data$bc <- seu_NP@meta.data$bc
p1p4@meta.data$short_bc <- seu_NP@meta.data$short_bc
p1p4@meta.data$Code_visit <- paste(p1p4@meta.data$Code,p1p4@meta.data$visit,sep = "_")
metadata$DonorID <- gsub("-.*","",metadata$donor_ID)
metadata$visit <- gsub(".*_","",metadata$Code_visit)
metadata$DonorID_visit <- paste(metadata$DonorID,metadata$visit,sep = "_")

library(stringi)
p1p4@meta.data$Donor_id <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                                  pattern=metadata$DonorID_visit,
                                                  replacement=metadata$Code_visit,
                                                  vectorize=FALSE)

p1p4@meta.data$Age <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                             pattern=metadata$DonorID_visit,
                                             replacement=metadata$age,
                                             vectorize=FALSE)

p1p4@meta.data$Sex <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                             pattern=metadata$DonorID_visit,
                                             replacement=metadata$sex,
                                             vectorize=FALSE)

p1p4@meta.data$study <- stri_replace_all_regex(p1p4@meta.data$Code_visit,
                                               pattern=metadata$DonorID_visit,
                                               replacement=metadata$study,
                                               vectorize=FALSE)

p1p4@meta.data$orig.ident <- paste(p1p4@meta.data$Code_visit,p1p4@meta.data$Age,p1p4@meta.data$Sex,p1p4@meta.data$Donor_id,sep="_")

p1p4[["percent.mt"]] <- PercentageFeatureSet(p1p4, pattern = "^MT-")
p1p4$mitoRatio <- p1p4@meta.data$percent.mt / 100

# Changing the Idents
pdf(paste(savedir,"p1p4_QC_mt_percent.pdf",sep = ""), width = 20, height = 7)
VlnPlot(p1p4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.000005, ncol = 3)
dev.off()

pdf(paste(savedir,"p1p4_QC_mt_percent_study.pdf",sep = ""), width = 20, height = 7)
VlnPlot(p1p4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "study", pt.size = 0.000005, ncol = 3)
dev.off()

cross_cellnames = rownames(p1p4@meta.data)[grep("cross",p1p4@meta.data$study)]
cross_seu <- subset(p1p4,cells=cross_cellnames)

cross_seu_all <- NormalizeData(cross_seu, normalization.method = "LogNormalize", scale.factor = 10000)
cross_seu_all <- FindVariableFeatures(cross_seu_all, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cross_seu_all)
cross_seu_all <- ScaleData(cross_seu_all, features = all.genes)
cross_seu_all <- RunPCA(cross_seu_all, features = VariableFeatures(object = cross_seu_all))
pdf(paste(savedir,"Elbow_plot.pdf",sep = ""))
ElbowPlot(cross_seu_all, ndims = 50)
dev.off()

dims <- c(1:20)
seed <- c(7456)
cross_seu_all <- RunTSNE(object = cross_seu_all, dims = dims,seed.use = seed)

cross_seu_all <- RunUMAP(cross_seu_all, dims = 1:20)
cross_seu_all <- FindNeighbors(cross_seu_all, dims = 1:20)

DefaultAssay(cross_seu_all) <- "RNA"
gene_zero <- as.vector(which(rowSums(cross_seu_all@assays$RNA@data)==0))
cross_seu_all_2 <- subset(cross_seu_all, features = rownames(cross_seu_all[-gene_zero,]))
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/python/python-3.9.2/bin/python3.9")
library(Rmagic)
py_discover_config("magic") # to check
cross_seu_all_impute <- magic(cross_seu_all_2, npca=20) ## imputing the vasc data
DefaultAssay(cross_seu_all_impute) <- "MAGIC_RNA"
dir.create(paste(savedir,"saveRDS_obj/",sep = ""))
saveRDS(cross_seu_all_impute,paste(savedir,"saveRDS_obj/cross_sectional_impute_all_genes.RDS",sep = ""))

dir.create(paste(savedir,"Vlnplot",sep = ""),showWarnings = FALSE)

all(rownames(cross_seu_all_impute@meta.data) == rownames(cross_seu_impute_ct@meta.data))
cross_seu_all_impute@meta.data$seurat_clusters <- cross_seu_impute_ct@meta.data$seurat_clusters

DefaultAssay(cross_seu_all_impute) <- "RNA"
pdf(paste(savedir,"Vlnplot/CD4_boxplot.pdf",sep = ""),width = 10, height = 6)
print(VlnPlot(cross_seu_all_impute, features = "CD4",pt.size = 0,group.by = "seurat_clusters") + geom_boxplot())
dev.off()
  
DefaultAssay(cross_seu_all_impute) <- "MAGIC_RNA"
pdf(paste(savedir,"Vlnplot/CD4_MAGIC.pdf",sep = ""),width = 10, height = 6)
print(VlnPlot(cross_seu_all_impute, features = "CD4",pt.size = 0, group.by = "seurat_clusters") + geom_boxplot())
dev.off()

res = c(0.2,0.4,0.6,0.8,1,1.2)
for (i in 1:length(res)) {
  DefaultAssay(cross_seu_all_impute) <- "RNA"
  cross_seu_all_impute <- FindClusters(cross_seu_all_impute, resolution = res[i])
  
  pdf(paste(savedir,"UMAP/dimplot_pc_",res[i],".pdf",sep = ""),width = 6, height = 5.5)
  print(DimPlot(cross_seu_all_impute, reduction = 'pca', label = TRUE, label.size = 5))
  dev.off()
}

#### Reference mapping PBMC ####
## Using this vignette https://satijalab.org/seurat/articles/multimodal_reference_mapping.html 
## Since the multimodal has been normalized on sctransform we will use it.
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
reference <- LoadH5Seurat("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/Resource/pbmc_multimodal.h5seurat")
DefaultAssay(cross_seu_all_impute) = "RNA"
cross_seu_all_impute = SCTransform(cross_seu_all_impute)

anchors <- FindTransferAnchors(
  reference = reference,
  query = cross_seu_all_impute,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

# We then transfer cell type labels and protein data from the reference to the query. Additionally, we project the query data onto the UMAP structure of the reference.

cross_seu_all_impute <- MapQuery(
  anchorset = anchors,
  query = cross_seu_all_impute,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

cross_seu_all_impute@reductions$tsne <- cross_seu_impute_ct@reductions$tsne
cross_seu_all_impute@reductions$umap <- cross_seu_impute_ct@reductions$umap

pdf(paste(savedir,"cross_seu_predicted.pdf",sep = ""), width = 8, height = 6)
DimPlot(cross_seu_all_impute, reduction = "tsne", group.by = "predicted.celltype.l1", label = TRUE)
dev.off()

pdf(paste(savedir,"cross_seu_predicted_2.pdf",sep = ""), width = 8, height = 6)
DimPlot(cross_seu_all_impute, reduction = "tsne", group.by = "predicted.celltype.l2", label = TRUE)
dev.off()

dir.create(paste(savedir,"saveRDS_obj/",sep = ""))
saveRDS(cross_seu_all_impute, paste(savedir,"saveRDS_obj/cross_seu_all_impute.RDS",sep = ""))

### Trajectory Analysis ####
### SlingShot

# Performing the trajectory it is better to take into consideration the global structure rather than the local structure. So we use reduced dimension UMAP
# SlingShot: multiple branching lineages.
# Slingshot breaks the inference problem into two steps, we are able to make use of appropriate methods for each task and avoid the common 
# trade-off between stabilityand the flexibility to detect complex structures. 
# Using a cluster-based MST for lineage inference allows Slingshot to identify potentially complex global patterns in the data without 
# being overly sensitive to individual data points.
# And our novel simultaneous principal curves method for pseudotime inference extends the stability and robustness properties of principal
# curves to the case of multiple branching lineages.
# 
# 2) the inference of pseudotime variables for cells along each lineage (simultaneous principal curves) to smooth the curves. 
# multiple lineage inference into two stages:
# 1. Identification of lineages, i.e., ordered sets of cell clusters, where all lineages share a starting cluster and each leads to a 
# unique terminal cluster. This is achieved by constructing an MST on  clusters of cells.
# 2. For each lineage, identification of pseudotimes, i.e., a one-dimensional variable representing each cell’s transcriptional progression
# toward the terminal state.
# 
# Robustness to noise.
# The Monocle procedure, which constructs an MST on individual cells and orders them according to a PQ tree along the longest path of the 
# MST, was the least stable of the methods we compared. The path drawn by Monocle was highly variable and sensitive to even small amounts 
# of noise; this instability has been.
# 
# the vertices of the piecewise linear path drawn by the cluster-based MST, multiple cells will often be assigned identical pseudotimes, 
# corresponding to the value at the vertex. principal curve approach was the most stable method, but on more complex datasets, it has the 
# obvious limitation of only characterizing a single lineage. It is for this reason that we chose to extend principal curves to accommodate
# multiple branching lineages.
# 
# Slingshot allows the use of a shape-sensitive distance measure inspired by the Mahalanobis distance which scales the distance between 
# cluster centers based on the covariance structure of the two clusters.

# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html 
suppressPackageStartupMessages({
  library(slingshot); library(SingleCellExperiment)
  library(RColorBrewer); library(scales)
  library(viridis); library(UpSetR)
  library(pheatmap); library(msigdbr)
  library(fgsea); library(knitr)
  library(ggplot2); library(gridExtra)
  library(tradeSeq); library(Seurat);
  library(bioc2020trajectories)
})

# Save the objects as separate matrices for input in slingshot
## convert back to singleCellExperiment
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/scRNA_public/updated_with_raw_counts/analysis/cross_sectional/"
cross_seu <- readRDS(paste(savedir,"saveRDS_obj/cross_sectional_impute.RDS",sep = ""))
DefaultAssay(cross_seu) <- "RNA"
sce <- as.SingleCellExperiment(cross_seu, assay = "RNA")

# The question is: should we fit a separate trajectory for each condition? We might expect the trajectory itself to be changed by the 
# treatment if the treatment effect is systematically large. Otherwise, the treatment may impact the expression profile of some genes but 
# the overall trajectory will be preserved.

# We can calculate the imbalance score Regions with a high score indicate that the local cell distribution according to treatment label is unbalanced 
# compared the overall distribution.Here, we see that, while there are some small regions of imbalance, the global path along the development axis is well-balanced.
# This means that we can fit a global trajectory to the full dataset, so this approach ensures that our trajectory accounts for all cell types present in the overall data.

# Since in this there is no condition so we will not use this funciton
# scores <- bioc2020trajectories::imbalance_score(
#   rd = reducedDims(sce)$UMAP, 
#   cl = colData(sce)$pheno$treatment_id,
#   k = 20, smooth = 40)

# The goal of slingshot is to use clusters of cells to uncover global structure and convert this structure into smooth lineages represented
# by one-dimensional variables, called “pseudotime.” We provide tools for learning cluster relationships in an unsupervised or semi-supervised 
# manner and constructing smooth curves representing each lineage,

# The minimal input to slingshot is a matrix representing the cells in a reduced-dimensional space and a vector of cluster labels. 
# With these two inputs, we then:
# 1. Identify the global lineage structure by constructing an minimum spanning tree (MST) on the clusters, with the getLineages function.
# 2. Construct smooth lineages and infer pseudotime variables by fitting simultaneous principal curves with the getCurves function.
# 3. Assess the output of each step with built-in visualization tools.

# However, we recommend clustering the cells even in datasets where only a single lineage is expected, as it allows for the potential discovery of novel branching events.
# The clusters identified in this step will be used to determine the global structure of the underlying lineages (that is, their number, 
# when they branch off from one another, and the approximate locations of those branching events). This is different than the typical goal
# of clustering single-cell data, which is to identify all biologically relevant cell types present in the dataset. For example, when 
# determining global lineage structure, there is no need to distinguish between immature and mature neurons since both cell types will, 
# presumably, fall along the same segment of a lineage.

# These two steps can be run separately with the getLineages and getCurves functions, or together with the wrapper function, slingshot (recommended). 

# sce_auto <- slingshot(sce,
#                  clusterLabels = 'seurat_clusters', 
#                  reducedDim = 'UMAP', approx_points = 100)
# 
# library(grDevices)
# colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
# plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
# 
# pdf(paste(savedir,"Trajectory/slingshot_trajectory.pdf",sep = ""))
# plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
# lines(SlingshotDataSet(sce), lwd=2, col='black')
# dev.off()
# the outgroup is an artificial cluster that is equidistant from all real clusters at some threshold value. If the original MST sans the 
# outgroup contains an edge that is longer than twice the threshold, the addition of the outgroup will cause the MST to instead be routed 
# through the outgroup. 


sce5 <- slingshot(sce,
                 clusterLabels = 'seurat_clusters', 
                 reducedDim = 'UMAP', approx_points = 100,
                 start.clus = 5,
                 omega = TRUE,
                 omega_scale=1.5)

# library(grDevices)
# colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
# plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
# 
# pdf(paste(savedir,"Trajectory/slingshot_trajectory_5.pdf",sep = ""))
# plot(reducedDims(sce5)$UMAP, col = plotcol, pch=16, asp = 1)
# lines(SlingshotDataSet(sce5), lwd=2, col='black')
# dev.off()

sce2 <- slingshot(sce, 
                  clusterLabels = 'seurat_clusters', 
                  reducedDim = 'UMAP', approx_points = 100,
                  start.clus = 2,
                  omega = TRUE,
                  omega_scale=1.5)

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
pdf(paste(savedir,"Trajectory/slingshot_trajectory_2.pdf",sep = ""))
plot(reducedDims(sce2)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce2), lwd=2, col='black')
lines(SlingshotDataSet(sce2), linInd = 2, type='l')
dev.off()

library(ArchR)
library(scater)
embedded_orig <- embedCurves(sce2, "TSNE")
rm(plot_list)
plot_list <- list()

for (i in 1:12) {
  embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
  embedded <- data.frame(embedded$s[embedded$ord,])
  g <- plotTSNE(sce2, colour_by=paste("slingPseudotime_",i,sep = ""))
  colnames(g$data) <- c("tSNE1","tSNE2",paste("Lineage",i,sep = ""))
  p <- ggplot(g$data, aes_string("tSNE1","tSNE2",color=paste("Lineage",i,sep = ""))) + geom_point(size=0.01) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra) + 
    geom_path(data=embedded, aes(x=tSNE_1, y=tSNE_2), color="black", size=1.2) + theme_bw()
  
  plot_list[[i]] <- p
}

require(gridExtra)
pdf(paste(savedir,"Trajectory/sce2_TSNE_splitted.pdf",sep = ""), width = 12, height = 12)
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
             plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
             nrow=4,ncol=3)
dev.off()


library(ArchR)
library(scater)
embedded_orig <- embedCurves(sce2, "UMAP")
rm(plot_list)
plot_list <- list()

for (i in 1:12) {
  embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
  embedded <- data.frame(embedded$s[embedded$ord,])
  g <- plotUMAP(sce2, colour_by=paste("slingPseudotime_",i,sep = ""))
  colnames(g$data) <- c("UMAP1","UMAP2",paste("Lineage",i,sep = ""))
  p <- ggplot(g$data, aes_string("UMAP1","UMAP2",color=paste("Lineage",i,sep = ""))) + geom_point(size=0.01) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra) + 
    geom_path(data=embedded, aes(x=UMAP_1, y=UMAP_2), color="black", size=1.2) + theme_bw()
  
  plot_list[[i]] <- p
}

require(gridExtra)
pdf(paste(savedir,"Trajectory/sce2_UMAP_splitted.pdf",sep = ""), width = 12, height = 12)
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
             plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
             nrow=4,ncol=3)
dev.off()

library(ArchR)
library(scater)
embedded_orig <- embedCurves(sce5, "TSNE")
rm(plot_list)
plot_list <- list()

for (i in 1:12) {
  embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
  embedded <- data.frame(embedded$s[embedded$ord,])
  g <- plotTSNE(sce5, colour_by=paste("slingPseudotime_",i,sep = ""))
  colnames(g$data) <- c("tSNE1","tSNE2",paste("Lineage",i,sep = ""))
  stopifnot(all(rownames(colData(sce5)) == rownames(g$data)))
  p <- ggplot(g$data, aes_string("tSNE1","tSNE2",color=paste("Lineage",i,sep = ""))) + geom_point(size=0.01) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra) + 
    geom_path(data=embedded, aes(x=tSNE_1, y=tSNE_2), color="black", size=1.2) + theme_bw()
  
  plot_list[[i]] <- p
}

require(gridExtra)
pdf(paste(savedir,"Trajectory/sce5_TSNE_splitted.pdf",sep = ""), width = 12, height = 12)
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
             plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
             nrow=4,ncol=3)
dev.off()

library(ArchR)
library(scater)
embedded_orig <- embedCurves(sce5, "UMAP")
rm(plot_list)
plot_list <- list()

for (i in 1:12) {
  embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
  embedded <- data.frame(embedded$s[embedded$ord,])
  g <- plotUMAP(sce5, colour_by=paste("slingPseudotime_",i,sep = ""))
  colnames(g$data) <- c("UMAP1","UMAP2",paste("Lineage",i,sep = ""))
  p <- ggplot(g$data, aes_string("UMAP1","UMAP2",color=paste("Lineage",i,sep = ""))) + geom_point(size=0.01) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra) + 
    geom_path(data=embedded, aes(x=UMAP_1, y=UMAP_2), color="black", size=1.2) + theme_bw()
  
  plot_list[[i]] <- p
}

require(gridExtra)
pdf(paste(savedir,"Trajectory/sce5_UMAP_splitted.pdf",sep = ""), width = 12, height = 12)
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
             plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
             nrow=4,ncol=3)
dev.off()

