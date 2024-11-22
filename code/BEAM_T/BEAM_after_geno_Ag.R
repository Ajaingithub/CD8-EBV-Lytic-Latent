#### only for Ag
### Performing analysis extracting each Ags
CD8_metadata <- CD8@meta.data
Ags <- grep(",",unique(CD8@meta.data$Ag_range_10),invert=TRUE,value=TRUE)

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/clones/"
dir.create(savedir, showWarnings = FALSE)

for (i in 1:length(Ags)) {
  dir.create(paste(savedir,Ags[i],sep = ""), showWarnings = FALSE)
  Ag_meta <- CD8_metadata[grep(Ags[i], CD8_metadata$Ag_range_10),]
  Ag_meta_O <- table(Ag_meta[grep("O",Ag_meta$Age),"TRA1_or_TRA2_TRB1_no_TRB2"])
  write.table(Ag_meta_O, paste(savedir,Ags[i],"/O_clone_freq.txt",sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
  Ag_meta_Y <- table(Ag_meta[grep("Y",Ag_meta$Age),"TRA1_or_TRA2_TRB1_no_TRB2"])
  write.table(Ag_meta_Y, paste(savedir,Ags[i],"/Y_clone_freq.txt",sep = ""), quote = F, col.names = T, row.names = F, sep = "\t")
}

#### Performing Differential for each Ag with other for K means clustering 
# memory_cells <- rownames(CD8@meta.data[grep("naïve",CD8@meta.data$celltypes,invert=TRUE),])
# CD8_mem <- subset(CD8, cells = memory_cells)
library(Seurat)
library(dplyr)
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/"
CD8 <- readRDS(paste(maindir, "saveRDS_obj/CD8_modality_integrate_cluster_0.6.RDS",sep = ""))
CD8@meta.data$TRB2 <- gsub(".*.-","",CD8@meta.data$TRA1_TRA2_TRB1_TRB2_cdrnt)
CD8@meta.data$TRB2 <- na_if(CD8@meta.data$TRB2, "NA")
TRB2_remove_cells <- rownames(CD8@meta.data[!is.na(CD8@meta.data$TRB2),])
memory_cells <- rownames(CD8@meta.data[grep("naïve",CD8@meta.data$celltypes,invert=TRUE),])
memory_cells_rem_TRB2 <- grep(paste(TRB2_remove_cells, collapse="|"), memory_cells, invert=TRUE, value=TRUE)
CD8_mem <- subset(CD8, cells = memory_cells_rem_TRB2)
CD8_mem@meta.data$Age_Ag_range_10 <- paste(CD8_mem@meta.data$Age, CD8_mem@meta.data$Ag_range_10,sep="_")
mem_cells_specific_Ag <- rownames(CD8_mem@meta.data[grep(",",CD8_mem@meta.data$Age_Ag_range_10,invert=TRUE),])
CD8_mem_spec_Ag <- subset(CD8_mem, cells = mem_cells_specific_Ag)
CD8_mem_spec_Ag@meta.data$Antigen2 <- gsub(".*._","",CD8_mem_spec_Ag@meta.data$Ag_range_10)
CD8_mem_spec_Ag@meta.data$AgeAntigen2 <- paste(CD8_mem_spec_Ag@meta.data$Age,CD8_mem_spec_Ag@meta.data$Antigen2,sep=".")

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)

Ag = unique(CD8_mem_spec_Ag@meta.data$Antigen2)
AgeAg <- unique(CD8_mem_spec_Ag@meta.data$AgeAntigen2)

for (i in 1:(length(AgeAg)-1)) {
  total_size = length(AgeAg)
  for (j in (i+1):total_size) {
  try({
    cluster = "all"
    group1 = AgeAg[i]
    group2 = AgeAg[j]
    remove_samples = c("O",grep(paste(AgeAg[i],"|",AgeAg[j],sep = ""),AgeAg,value = TRUE, invert = TRUE))
    remove_samples_name <- paste(remove_samples,collapse="_and_")
    DefaultAssay(CD8_mem_spec_Ag) <- "RNA"
    # savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run01_to_10_analysis/vaccine/S_vs_Z/"
    source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/COVID_pseudobulk_PCA_within_AJ_2.R")
    CD8_subset_dds <- COVID_pseudobulk_within_cluster_AJ(obj = CD8_mem_spec_Ag, savedir = savedir, group1 = group1, group2 = group2,
                                                         grouping_by = "AgeAntigen2", cluster = cluster,cell_freq = 20,remove_samples = remove_samples,
                                                         cluster_group = "seurat_clusters",sample_col = "orig.ident", batch_col = "Run",
                                                         gene_min_counts =  5, 
                                                         column_name_split = c("sampleid","virus","age","age_number","gender","Run"))
    
    cluster2 <- paste(cluster, sep="_", collapse="_")
    savedir2 <- paste(savedir,"pseudobulk/clus_",cluster2,"_",group1,"_vs_",group2,"/",sep = "")
    # savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run04/Analysis/pseudobulk/clus_13_removed__DMSO_vs_Sprotein/"
    source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/RNAseq/pipeline_function/RNAseq_limma_EdgeR.R")
    design0 <- model.matrix(~ 0 + AgeAntigen2 + Run, data = colData(CD8_subset_dds))
    colnames(design0) <- c(AgeAg[i],AgeAg[j],paste("Run",1:(ncol(design0)-2),sep = ""))
    differential <- paste(AgeAg[i],"_VS_",AgeAg[j],"=",AgeAg[i],"-",AgeAg[j],sep = "")
    print(differential)
    cm <- makeContrasts(differential,levels = design0)
    desl_clus <- LimmaEdgeR_differential(dds = CD8_subset_dds,
                                         design0 = design0,
                                         cm = cm,
                                         savedir = savedir2,
                                         logfc = 0.5,
                                         p_value_adj = 0.05)
  })
  }
}

#### Generating for the BMRF1

Ag = unique(CD8_mem_spec_Ag@meta.data$Antigen2)
AgeAg <- unique(CD8_mem_spec_Ag@meta.data$AgeAntigen2)

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Range_10_rem_TRB2_rem_batch_effect/"
cluster = "all"
group1 = "O.BMRF1"
group2 = "Y.BMRF1"
remove_samples = c(grep(paste(group1,"|",group2,sep = ""),AgeAg,value = TRUE, invert = TRUE))
remove_samples_name <- paste(remove_samples,collapse="_and_")
DefaultAssay(CD8_mem_spec_Ag) <- "RNA"
# savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run01_to_10_analysis/vaccine/S_vs_Z/"
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/COVID_pseudobulk_PCA_within_AJ_2.R")
CD8_subset_dds <- COVID_pseudobulk_within_cluster_AJ(obj = CD8_mem_spec_Ag, savedir = savedir, group1 = group1, group2 = group2,
                                                     grouping_by = "AgeAntigen2", cluster = cluster,cell_freq = 20,remove_samples = remove_samples,
                                                     cluster_group = "seurat_clusters",sample_col = "orig.ident", batch_col = "Run",
                                                     gene_min_counts =  5, 
                                                     column_name_split = c("sampleid","virus","age","age_number","gender","Run"), sample_num = 2)

cluster2 <- paste(cluster, sep="_", collapse="_")
savedir2 <- paste(savedir,"pseudobulk/clus_",cluster2,"_",group1,"_vs_",group2,"/",sep = "")
# savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run04/Analysis/pseudobulk/clus_13_removed__DMSO_vs_Sprotein/"
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/RNAseq/pipeline_function/RNAseq_limma_EdgeR.R")
design0 <- model.matrix(~ 0 + AgeAntigen2 + Run, data = colData(CD8_subset_dds))
colnames(design0) <- c("O.BMRF1","Y.BMRF1",paste("Run",1:(ncol(design0)-2),sep = ""))
differential <- paste("O.BMRF1","_VS_","Y.BMRF1","=","O.BMRF1","-","Y.BMRF1",sep = "")
print(differential)
cm <- makeContrasts(differential,levels = design0)
desl_clus <- LimmaEdgeR_differential(dds = CD8_subset_dds,
                                     design0 = design0,
                                     cm = cm,
                                     savedir = savedir2,
                                     logfc = 0.5,
                                     p_value_adj = 0.05,
                                     column_name_split = c("sampleid","virus","age","age_number","gender","Run"),
                                     grouping_by = "AgeAntigen2")


### Since some of the cases we need to remove those samples in which individuals are less than 2 in either of the group
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/pseudobulk/"
filenames <- list.files(maindir,recursive = TRUE, pattern = "CPM_counts.txt",full.names = FALSE)
setwd(maindir)

df2 <- as.data.frame(matrix(nrow = length(filenames), ncol = 2))
colnames(df2) <- c("grp1","grp2")
rownames(df2) <- filenames

for (i in 1:length(filenames)) {
  try({
  samplename <- colnames(read.table(filenames[i]))
  df2[i,] <- as.data.frame(table(gsub("*.*Exp0[0-9]_","",samplename)))[,2]
  })
}

low_number <- df2[df2$grp1 < 3 | df2$grp2 < 3,]
low_num_name <- gsub("/.*.","",rownames(low_number))
write.table(low_num_name,paste("sample_less_than3.txt"),row.names = F, col.names = F, sep = "\t", quote = F)


### K means clustering #####
### Now we have genes that are more than 3 individuals in each group and we have performed uniq for the genes that is required
### Since we have performed the same analysis with the ATACseq 
### Based on the analysis 90 samples have more than 20 cells and found to involve in differentiation so further taking only those samples
required_samples <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/pseudobulk/required_samples")[,1]
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/"
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/specific_samples_batch_effect_removed_BMRF1_LMP/"
diff_genes <- read.table(paste(maindir,"pseudobulk/final_genes_sig_uniq.txt",sep = ""))[,1]
Ind_Antigen <- unique(paste(CD8_mem_spec_Ag@meta.data$orig.ident,"_",CD8_mem_spec_Ag@meta.data$Age,".", gsub("EBV_|VZV_","",CD8_mem_spec_Ag@meta.data$Ag_range_10),sep=""))

### Since we need the groups from each of the individuals
remove_samples = grep(paste(required_samples,collapse = "|"), Ind_Antigen, invert = TRUE, value = TRUE)
remove_samples = c(remove_samples,"BMRF1","LMP","Y3_EBV_Y_28_F_Exp02_Y.BMLF1")
remove_samples_name <- paste(remove_samples,collapse="_and_")
DefaultAssay(CD8_mem_spec_Ag) <- "RNA"
# savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run01_to_10_analysis/vaccine/S_vs_Z/"
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/COVID_pseudobulk_PCA_within_AJ_2.R")
CD8_subset_rld <- COVID_pseudobulk_within_cluster_AJ(obj = CD8_mem_spec_Ag, savedir = savedir, group1 = "all", group2 = "all",
                                                     grouping_by = "AgeAntigen2", cluster = "all", cell_freq = 20,remove_samples = remove_samples,
                                                     cluster_group = "seurat_clusters",sample_col = "orig.ident", batch_col = "Run",
                                                     gene_min_counts =  5, 
                                                     column_name_split = c("sampleid","virus","age","age_number","gender","Run"))

saveRDS(CD8_subset_rld, paste(savedir,"CD8_subset_rld.RDS",sep = ""))
# sampleinfo$Group <- factor(sampleinfo$AgeAntigen2)

## Run this in the 
### Run this in r/4.2.2
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/specific_samples_batch_effect_removed_BMRF1_LMP/"
CD8_subset_rld <- readRDS(paste(savedir,"CD8_subset_rld.RDS",sep = ""))
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/"
diff_genes <- read.table(paste(maindir,"pseudobulk/final_genes_sig_uniq.txt",sep = ""))[,1]

sampleinfo <- colData(CD8_subset_rld)
dir.create(savedir, showWarnings = FALSE)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/ATACseq/k_mer_clustering.R")
gap_stat_2 <- kmeans_clustering_ATAC(differential_peaks = diff_genes, rlog_dds = CD8_subset_rld,
                       sampleinfo = sampleinfo, Group = "AgeAntigen2", saveDir = savedir)

#### We are performing without removing the batch effect also
# savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/no_effect_remove"
dir.create(savedir,showWarnings=FALSE)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/ATACseq/k_mer_clustering.R")
kmeans_clustering_ATAC(differential_peaks = diff_genes, rlog_dds = CD8_subset_rld,
                       sampleinfo = sampleinfo, Group = "AgeAntigen2", saveDir = savedir)

### Now plotting it
# gap_stat_2 <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/Differential/cluster_gap.RDS")
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/ATACseq/k_means_plots.R")
dir.create(paste(savedir,"Clustering",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"Clustering/cluster_number.pdf",sep = ""), width = 7, height = 7)
print(fviz_gap_stat(gap_stat_2))
dev.off()


df_scale <- read.table(paste(savedir,"Differential/median_counts_for_Group_significant_normalized.txt",sep = ""))
set.seed(123)
final_2 <- kmeans(df_scale, 5, nstart = 25) # based on the optimal clusters
print(final_2)

saveRDS(final_2,paste(savedir,"Clustering/final_2.Rds",sep = ""))


peaks_arranged <- as.data.frame(final_2$cluster[order(as.data.frame(final_2$cluster))])
colnames(peaks_arranged) <- "clusters"
peaks_arranged$clusters <- paste("C",peaks_arranged$clusters,sep="")
table(peaks_arranged$clusters)

# C1   C2   C3   C4   C5
# 1017 1511 1846 1571 1435


df_scale_ordered <- df_scale[match(rownames(peaks_arranged), rownames(df_scale)),]

library(ComplexHeatmap)
library(circlize)
h=Heatmap(df_scale_ordered,
          name = "scaled", #title of legend
          column_title = "Groups", row_title = "Genes",
          row_names_gp = gpar(fontsize = 0.0001), # Text size for row names
          column_names_gp = gpar(fontsize = 10),
          cluster_columns = TRUE,
          cluster_rows = FALSE,
          col=colorRamp2(c(-0.5, 0, 0.5), 
                         c("blue", "white", "red"))
          #rect_gp = gpar(col = "white", lwd = 1),
          # right_annotation = row_anno,
          # top_annotation = column_anno
          #column_split = factor(c(rep("A", 9), rep("B",8)), levels = c("A","B")),
)

pdf(paste(savedir,"Clustering/heatmap_5_clusters.pdf",sep = ""))
h
dev.off()

write.table(peaks_arranged, paste(savedir,"genes_order.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")


### Performing the Kmeans for O vs Y ####
required_samples <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/pseudobulk/required_samples")[,1]
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/"
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/O_vs_Y_all_Ag/"
dir.create(savedir, showWarnings = FALSE)
diff_genes <- read.table(paste(maindir,"pseudobulk/final_sig_genes_Y_vs_O.txt",sep = ""))[,1]
Ind_Antigen <- unique(paste(CD8_mem_spec_Ag@meta.data$orig.ident,"_",CD8_mem_spec_Ag@meta.data$Age,".", gsub("EBV_|VZV_","",CD8_mem_spec_Ag@meta.data$Ag_range_10),sep=""))

### Since we need the groups from each of the individuals
# remove_samples = grep(paste(required_samples,collapse = "|"), Ind_Antigen, invert = TRUE, value = TRUE)
# remove_samples = c(remove_samples)
# remove_samples_name <- paste(remove_samples,collapse="_and_")
remove_samples = c("BALF4")
DefaultAssay(CD8_mem_spec_Ag) <- "RNA"
# savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run01_to_10_analysis/vaccine/S_vs_Z/"
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/COVID_pseudobulk_PCA_within_AJ_2.R")
CD8_subset_rld <- COVID_pseudobulk_within_cluster_AJ(obj = CD8_mem_spec_Ag, savedir = savedir, group1 = "all_Ag", group2 = "all_Ag",
                                                     grouping_by = "AgeAntigen2", cluster = "all", cell_freq = 20,remove_samples = remove_samples,
                                                     cluster_group = "seurat_clusters",sample_col = "orig.ident", batch_col = "Run",
                                                     gene_min_counts =  5,
                                                     column_name_split = c("sampleid","virus","age","age_number","gender","Run"))

saveRDS(CD8_subset_rld, paste(savedir,"CD8_subset_rld.RDS",sep = ""))
# sampleinfo$Group <- factor(sampleinfo$AgeAntigen2)


### Please perform kmeans outside the conda environment
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/pseudobulk/final_required_comparison/"
CD8_subset_rld <- readRDS(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/O_vs_Y_all_Ag/CD8_subset_rld.RDS",sep = ""))
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/pseudobulk/final_required_comparison/"
diff_genes <- read.table(paste(maindir,"final_gene_list",sep = ""))[,1]

sampleinfo <- colData(CD8_subset_rld)
dir.create(savedir, showWarnings = FALSE)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/ATACseq/k_mer_clustering.R")
gap_stat_2 <- kmeans_clustering_ATAC(differential_peaks = diff_genes, rlog_dds = CD8_subset_rld,
                                     sampleinfo = sampleinfo, Group = "AgeAntigen2", saveDir = savedir)

#### We are performing without removing the batch effect also
# savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/no_effect_remove"
# dir.create(savedir,showWarnings=FALSE)
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/ATACseq/k_mer_clustering.R")
# kmeans_clustering_ATAC(differential_peaks = diff_genes, rlog_dds = CD8_subset_rld,
#                        sampleinfo = sampleinfo, Group = "AgeAntigen2", saveDir = savedir)

### Now plotting it
# gap_stat_2 <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/Differential/cluster_gap.RDS")
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/ATACseq/k_means_plots.R")
dir.create(paste(savedir,"Clustering",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"Clustering/cluster_number.pdf",sep = ""), width = 7, height = 7)
print(fviz_gap_stat(gap_stat_2))
dev.off()

clusters_num <- c(3,5,7)
for (i in 1:length(clusters_num)) {
  df_scale <- read.table(paste(savedir,"Differential/median_counts_for_Group_significant_normalized.txt",sep = ""))
  set.seed(123)
  final_2 <- kmeans(df_scale, clusters_num[i], nstart = 25) # based on the optimal clusters
  print(final_2)
  
  dir.create(paste(savedir,"Clustering",sep = ""))
  saveRDS(final_2,paste(savedir,"Clustering/final_clusters_",clusters_num[i],".Rds",sep = ""))

  peaks_arranged <- as.data.frame(final_2$cluster[order(as.data.frame(final_2$cluster))])
  colnames(peaks_arranged) <- "clusters"
  peaks_arranged$clusters <- paste("C",peaks_arranged$clusters,sep="")
  print(table(peaks_arranged$clusters))
  
  
  # C1  C2  C3  C4  C5
  # 322 246 269 289 293
  
  df_scale_ordered <- df_scale[match(rownames(peaks_arranged), rownames(df_scale)),]
  
  library(ComplexHeatmap)
  library(circlize)
  h=Heatmap(df_scale_ordered,
            name = "scaled", #title of legend
            column_title = "Groups", row_title = "Genes",
            row_names_gp = gpar(fontsize = 0.0001), # Text size for row names
            column_names_gp = gpar(fontsize = 10),
            cluster_columns = TRUE,
            cluster_rows = FALSE,
            col=colorRamp2(c(-0.5, 0, 0.5), 
                           c("blue", "white", "red"))
            #rect_gp = gpar(col = "white", lwd = 1),
            # right_annotation = row_anno,
            # top_annotation = column_anno
            #column_split = factor(c(rep("A", 9), rep("B",8)), levels = c("A","B")),
  )
  
  pdf(paste(savedir,"Clustering/heatmap_",clusters_num[i],"_clusters_change_scale_clustered.pdf",sep = ""))
  print(h)
  dev.off()
  
  ### Ordering the columns
  col_order <- c("Y.BMLF1", "O.BMLF1", "Y.BRLF1", "O.BRLF1", "Y.BMRF1", "O.BMRF1", "Y.EBNA1", "O.EBNA1", "Y.EBNA3C", "O.EBNA3C", "Y.LMP1", "O.LMP1", "Y.LMP2", "O.LMP2", "Y.IE62", "O.IE62")
  df_scale_ordered <- df_scale[match(rownames(peaks_arranged), rownames(df_scale)),match(col_order, colnames(df_scale))]
  
  library(ComplexHeatmap)
  library(circlize)
  h=Heatmap(df_scale_ordered,
            name = "scaled", #title of legend
            column_title = "Groups", row_title = "Genes",
            row_names_gp = gpar(fontsize = 0.0001), # Text size for row names
            column_names_gp = gpar(fontsize = 10),
            cluster_columns = FALSE,
            cluster_rows = FALSE,
            col=colorRamp2(c(-0.5, 0, 0.5), 
                           c("blue", "white", "red"))
            #rect_gp = gpar(col = "white", lwd = 1),
            # right_annotation = row_anno,
            # top_annotation = column_anno
            #column_split = factor(c(rep("A", 9), rep("B",8)), levels = c("A","B")),
  )
  
  pdf(paste(savedir,"Clustering/heatmap_",clusters_num[i],"_clusters_change_scale.pdf",sep = ""))
  print(h)
  dev.off()
  
  write.table(peaks_arranged, paste(savedir,"genes_order_",clusters_num[i],".txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
}


#### For visualization purposes we have changed the order of the clusters
df_scale <- read.table(paste(savedir,"Differential/median_counts_for_Group_significant_normalized.txt",sep = ""))

peaks_arranged <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/pseudobulk/final_required_comparison/gene_order_new_format.txt")

colnames(peaks_arranged) <- "clusters"
print(table(peaks_arranged$clusters))


# C1  C2  C3  C4  C5
# 322 246 269 289 293

df_scale_ordered <- df_scale[match(rownames(peaks_arranged), rownames(df_scale)),]

library(ComplexHeatmap)
library(circlize)
h=Heatmap(df_scale_ordered,
          name = "scaled", #title of legend
          column_title = "Groups", row_title = "Genes",
          row_names_gp = gpar(fontsize = 0.0001), # Text size for row names
          column_names_gp = gpar(fontsize = 10),
          cluster_columns = TRUE,
          cluster_rows = FALSE,
          col=colorRamp2(c(-0.5, 0, 0.5), 
                         c("blue", "white", "red"))
          #rect_gp = gpar(col = "white", lwd = 1),
          # right_annotation = row_anno,
          # top_annotation = column_anno
          #column_split = factor(c(rep("A", 9), rep("B",8)), levels = c("A","B")),
)

pdf(paste(savedir,"Clustering/heatmap_",clusters_num[i],"_clusters_change_scale_clustered_ordered.pdf",sep = ""))
print(h)
dev.off()

### Ordering the columns
col_order <- c("Y.BMLF1", "O.BMLF1", "Y.BRLF1", "O.BRLF1", "Y.BMRF1", "O.BMRF1", "Y.EBNA1", "O.EBNA1", "Y.EBNA3C", "O.EBNA3C", "Y.LMP1", "O.LMP1", "Y.LMP2", "O.LMP2", "Y.IE62", "O.IE62")
df_scale_ordered <- df_scale[match(rownames(peaks_arranged), rownames(df_scale)),match(col_order, colnames(df_scale))]

library(ComplexHeatmap)
library(circlize)
h=Heatmap(df_scale_ordered,
          name = "scaled", #title of legend
          column_title = "Groups", row_title = "Genes",
          row_names_gp = gpar(fontsize = 0.0001), # Text size for row names
          column_names_gp = gpar(fontsize = 10),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          col=colorRamp2(c(-0.5, 0, 0.5), 
                         c("blue", "white", "red"))
          #rect_gp = gpar(col = "white", lwd = 1),
          # right_annotation = row_anno,
          # top_annotation = column_anno
          #column_split = factor(c(rep("A", 9), rep("B",8)), levels = c("A","B")),
)

pdf(paste(savedir,"Clustering/heatmap_",clusters_num[i],"_clusters_change_scale.pdf",sep = ""))
print(h)
dev.off()

write.table(peaks_arranged, paste(savedir,"genes_order_",clusters_num[i],"_ordered.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")



### Performing the GSEA for the Cluster 5 for innate and adaptive
## Since it will be checked for each Antigen
### GSEA #####
library(fgsea)
library(data.table)
library(ggplot2)

### For the cluster 7 genes we have to take the Cluster 5, 2, 3 and 2+3 for innate and adaptive
### Starting with innateness
genes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/O_vs_Y_all_Ag/genes_order_7.txt")
clusters <- c(5,2,3)
genes$clusters2 <- genes$clusters

### combining cluster 2 and 3
genes$clusters2 <- gsub("2|3","2_3",genes$clusters2)
savedir = c("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/O_vs_Y_all_Ag/GSEA/")

dir.create(savedir, showWarnings = FALSE)
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Range_10_rem_TRB2_rem_batch_effect/pseudobulk/"
files <- grep("Y_vs_O",list.files(maindir, pattern = "_genes_no_filter.txt", full.names = TRUE, recursive = TRUE),invert = TRUE, value = TRUE)
Ag_name <- gsub(maindir,"",files) %>% gsub("/Differential/O_VS_Y_all_genes_no_filter.txt","",.) %>% gsub("/clus_all_O.","",.) %>% gsub("_vs_.*.","",.)

clusters <- "2_3"
for (i in 1:length(clusters)) {
  genes_cluster <- genes[grep(clusters[i],genes$clusters2),]
  dir.create(paste(savedir,"clus",clusters[i],sep = ""), showWarnings = FALSE)
  
  ### Genes to consider for GSEA
  req_genes <- rownames(genes_cluster)
  
  for (j in 1:length(Ag_name)) {
    Ag_diff <- read.table(files[j], header = TRUE, sep = "\t")
    
    ## extracting the genes required 
    Ag_diff_req <- Ag_diff[match(req_genes,rownames(Ag_diff), nomatch = 0),]
    Ag_diff_req_Entrez <-  bitr(rownames(Ag_diff_req), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")
    Ag_diff_req_Entrez[duplicated(Ag_diff_req_Entrez$SYMBOL)  # Identify the duplicates
                    | duplicated(Ag_diff_req_Entrez$SYMBOL, fromLast=TRUE),]
    # Ag_diff_req_Entrez_unique <- tops_all_Entrez[Ag_diff_req_Entrez$ENTREZID != 122394733 & Ag_diff_req_Entrez$ENTREZID != 100124696,]
    ## Adding the fold Change
    Ag_diff_req_Entrez$logFC <- Ag_diff_req[match(Ag_diff_req_Entrez$SYMBOL, rownames(Ag_diff_req)),"logFC"]
    
    ### Saving the file
    write.csv(Ag_diff_req_Entrez, file = paste(savedir,"clus",clusters[i],"/",Ag_name[j],"_diff_req_Entrez.csv",sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    library(clusterProfiler)
    geneList = Ag_diff_req_Entrez$logFC ## logfc
    names(geneList) = as.character(Ag_diff_req_Entrez$ENTREZID)
    geneList = sort(geneList, decreasing = TRUE)
    
    library(fgsea)
    library(data.table)
    library(ggplot2)
    
    T_cells <- c("adaptiveness","innateness")
    
    for (k in 1:length(T_cells)) {
      tfh_genes <- read.table(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/Resource/",T_cells[k],sep = ""))[,1]
      Tfh_entrez <-  bitr(tfh_genes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")
      Tfh_entrez[,2]
      pathways <- list()
      pathways[[T_cells[k]]] <- Tfh_entrez[,2]
      geneList
      
      fgseaRes <- fgsea(pathways = pathways, 
                        stats    = geneList,
                        minSize  = 8,
                        eps = 0.0,
                        maxSize  = 500)
      
      ### for this it 
      pdf(paste(savedir,"clus",clusters[i],"/",Ag_name[j],"_",T_cells[k],".pdf",sep = ""))
      print(plotEnrichment(pathways[[T_cells[k]]],geneList,gseaParam = 1) + labs(title=T_cells[k]))
      dev.off()
      
      pdf(paste(savedir,"clus",clusters[i],"/",Ag_name[j],"_",T_cells[k],"_table.pdf",sep = ""))
      print(plotGseaTable(pathways[T_cells[k]],geneList, fgseaRes,gseaParam=0.5))
      dev.off()
    }
  }
}

### GSEA for O vs Y ####
library(fgsea)
library(data.table)
library(ggplot2)
library(clusterProfiler)
library(dplyr)

### For the cluster 7 genes we have to take the Cluster 5, 2, 3 and 2+3 for innate and adaptive
### combining cluster 2 and 3
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/O_vs_Y_all_Ag/GSEA/all_genes/"
dir.create(savedir, showWarnings = FALSE)
maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Range_10_rem_TRB2_rem_batch_effect/pseudobulk/"
files <- grep("Y_vs_O",list.files(maindir, pattern = "_genes_no_filter.txt", full.names = TRUE, recursive = TRUE),invert = TRUE, value = TRUE)
Ag_name <- gsub(maindir,"",files) %>% gsub("/Differential/O_VS_Y_all_genes_no_filter.txt","",.) %>% gsub("/clus_all_O.","",.) %>% gsub("_vs_.*.","",.)

#j=4 
# add j= x number corresponding to number in Ag_name when only running a subset of comparisons, in that case skip the next line "for (....)" and keep running with "Ag_diff..."

for (j in 1:length(Ag_name)) {
  Ag_diff <- read.table(files[j], header = TRUE, sep = "\t")
  
  ## extracting the genes required 
  Ag_diff_req <- Ag_diff
  Ag_diff_req_Entrez <-  bitr(rownames(Ag_diff_req), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")
  Ag_diff_req_Entrez[duplicated(Ag_diff_req_Entrez$SYMBOL)  # Identify the duplicates
                     | duplicated(Ag_diff_req_Entrez$SYMBOL, fromLast=TRUE),]
  # Ag_diff_req_Entrez_unique <- tops_all_Entrez[Ag_diff_req_Entrez$ENTREZID != 122394733 & Ag_diff_req_Entrez$ENTREZID != 100124696,]
  ## Adding the fold Change
  Ag_diff_req_Entrez$logFC <- Ag_diff_req[match(Ag_diff_req_Entrez$SYMBOL, rownames(Ag_diff_req)),"logFC"]
  
  ### Saving the file
  
  write.csv(Ag_diff_req_Entrez, file = paste(savedir,Ag_name[j],"_all_O_vs_Y_diff_req_Entrez.csv",sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  library(clusterProfiler)
  geneList = Ag_diff_req_Entrez$logFC ## logfc
  names(geneList) = as.character(Ag_diff_req_Entrez$ENTREZID)
  geneList = sort(geneList, decreasing = TRUE)
  
  library(fgsea)
  library(data.table)
  library(ggplot2)
  
  T_cells <- c("GSEA_Exhausted-vs-Effector_CD8-LCMV_GSE9650.txt","GSEA_Exhausted-vs-Memory_CD8-LCMV_GSE41867.txt", "GSEA_Senescence_GOBP.txt","GSEA_Senescence_Reactome.txt", "GSEA_Senescence_SenMayo.txt", "GSEA_Senescence_SenSig_UP.txt")
  #T_cells <- c("adaptiveness","innateness")
  
  ## To distinguish different pathways please put the overall pathway name like for adaptiveness and innateness, adapt_innate
  filename_pathway <- "EX-SNC"
  
  # pathways_combined <- list()
  pathways <- list()
  for (k in 1:length(T_cells)) {
    tfh_genes <- read.table(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/Resource/",T_cells[k],sep = ""))[,1]
    Tfh_entrez <-  bitr(tfh_genes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")
    Tfh_entrez[,2]
    # 
    pathways[[T_cells[k]]] <- Tfh_entrez[,2]
    geneList
    
    fgseaRes <- fgsea(pathways = pathways,
                      stats    = geneList,
                      minSize  = 8,
                      eps = 0.0,
                      maxSize  = 1500)
    
    ## If you get the warning it is due to the  genes in your preranked list having the same numerical/ranking value. For example, if you rank genes by log2 fold change,
    # genes with the same log2 fold change value will be ordered randomly once processed. You can use rank() function and a ties.method parameter to control how genes 
    # with the same value are ranked:
    
    ### for this it 
    pdf(paste(savedir,Ag_name[j],"_",T_cells[k],"_all_O_vs_Y_2.pdf",sep = ""))
    print(plotEnrichment(pathways[[T_cells[k]]],geneList,gseaParam = 1) + labs(title=T_cells[k]))
    dev.off()
    
  }

  pdf(paste(savedir,Ag_name[j],"_",filename_pathway,"_all_O_vs_Y_table.pdf",sep = ""))
  print(plotGseaTable(pathways[T_cells],geneList, fgseaRes,gseaParam=0.5))
  dev.off()
  
  fgsea_df <- as.matrix(fgseaRes)
  write.table(fgsea_df, paste(savedir,Ag_name[j],"_",filename_pathway,"_all_O_vs_Y_table.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
}


### Multiple testing
T_cells <- c("GSEA_Exhausted-vs-Effector_CD8-LCMV_GSE9650.txt","GSEA_Exhausted-vs-Memory_CD8-LCMV_GSE41867.txt","GSEA_StemLikeExhausted-vs-TerminalExhausted_CD8-LCMV_GSE84105.txt","GSEA_TerminalExhausted-vs-StemLikeExhausted_CD8-LCMV_GSE84105.txt","GSEA_Senescence_GOBP.txt","GSEA_Senescence_Reactome.txt", "GSEA_Senescence_SenMayo.txt", "GSEA_Senescence_SenSig_UP.txt")
pathways <- list()
for (i in 1:length(T_cells)) {
  innate <- read.table(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/Resource/",T_cells[i],sep = ""))[,1]
  innate_entrez <-  bitr(innate, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")
  pathways[[T_cells[i]]] <- innate_entrez[,2]
}

fgseaRes <- fgsea(pathways = pathways, 
                  stats= geneList, 
                  minSize  = 8, 
                  eps = 0.0, maxSize  = 1500)

### for this it 
pdf(paste(savedir,Ag_name[j],"_",T_cells[k],"_all_O_vs_Y_2.pdf",sep = ""))
print(plotEnrichment(pathways[["adaptiveness"]],geneList) + labs(title=T_cells[k]))
dev.off()

pdf(paste(savedir,Ag_name[j],"_",T_cells[k],"_all_O_vs_Y_table_2.pdf",sep = ""))
print(plotGseaTable(pathways[T_cells[k]],geneList, fgseaRes,gseaParam=0.5))
dev.off()



### Y for all naive ####
required_samples <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/pseudobulk/required_samples")[,1]
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/"
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/Y_all_Ags_comp/Naive/"
dir.create(savedir, showWarnings = FALSE)
diff_genes <- read.table(paste(maindir,"pseudobulk/diff_genes_Y_all_Ags",sep = ""))[,1]
Ind_Antigen <- unique(paste(CD8_mem_spec_Ag@meta.data$orig.ident,"_",CD8_mem_spec_Ag@meta.data$Age,".", gsub("EBV_|VZV_","",CD8_mem_spec_Ag@meta.data$Ag_range_10),sep=""))

cellname_req_Ag <- rownames(CD8_all@meta.data[grep(",",CD8_all@meta.data$Age_Ag,invert=TRUE),])
CD8_all_req_Ag <- subset(CD8_all, cells=cellname_req_Ag)
CD8_all_req_Ag@meta.data$Age_Ag <- gsub("_VZV_|_EBV_",".",CD8_all_req_Ag@meta.data$Age_Ag)

### Since we need the groups from each of the individuals
# remove_samples = grep(paste(required_samples,collapse = "|"), Ind_Antigen, invert = TRUE, value = TRUE)
# remove_samples = c(remove_samples)
# remove_samples_name <- paste(remove_samples,collapse="_and_")
remove_samples = c("O","BALF4","Y7_EBV_Y_39_M_Exp03_Y.BMRF1")
DefaultAssay(CD8_all_req_Ag) <- "RNA"
# savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run01_to_10_analysis/vaccine/S_vs_Z/"
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/COVID_pseudobulk_PCA_within_AJ_2.R")
CD8_subset_rld <- COVID_pseudobulk_within_cluster_AJ(obj = CD8_all_req_Ag, savedir = savedir, group1 = "all_Ag", group2 = "all_Ag",
                                                     grouping_by = "Age_Ag", cluster = "all", cell_freq = 20,remove_samples = remove_samples,
                                                     cluster_group = "seurat_clusters",sample_col = "orig.ident", batch_col = "Run",
                                                     gene_min_counts =  5,
                                                     column_name_split = c("sampleid","virus","age","age_number","gender","Run"))

saveRDS(CD8_subset_rld, paste(savedir,"CD8_subset_rld.RDS",sep = ""))

#### Y comparison all Ags memory ####
required_samples <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/pseudobulk/required_samples")[,1]
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/"
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/Y_all_Ags_comp/"
dir.create(savedir, showWarnings = FALSE)
diff_genes <- read.table(paste(maindir,"pseudobulk/diff_genes_Y_all_Ags",sep = ""))[,1]
Ind_Antigen <- unique(paste(CD8_mem_spec_Ag@meta.data$orig.ident,"_",CD8_mem_spec_Ag@meta.data$Age,".", gsub("EBV_|VZV_","",CD8_mem_spec_Ag@meta.data$Ag_range_10),sep=""))

### Since we need the groups from each of the individuals
# remove_samples = grep(paste(required_samples,collapse = "|"), Ind_Antigen, invert = TRUE, value = TRUE)
# remove_samples = c(remove_samples)
# remove_samples_name <- paste(remove_samples,collapse="_and_")
remove_samples = c("O","BALF4")
DefaultAssay(CD8_mem_spec_Ag) <- "RNA"
# savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run01_to_10_analysis/vaccine/S_vs_Z/"
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/COVID_pseudobulk_PCA_within_AJ_2.R")
CD8_subset_rld <- COVID_pseudobulk_within_cluster_AJ(obj = CD8_mem_spec_Ag, savedir = savedir, group1 = "all_Ag", group2 = "all_Ag",
                                                     grouping_by = "AgeAntigen2", cluster = "all", cell_freq = 20,remove_samples = remove_samples,
                                                     cluster_group = "seurat_clusters",sample_col = "orig.ident", batch_col = "Run",
                                                     gene_min_counts =  5,
                                                     column_name_split = c("sampleid","virus","age","age_number","gender","Run"))

saveRDS(CD8_subset_rld, paste(savedir,"CD8_subset_rld.RDS",sep = ""))
# sampleinfo$Group <- factor(sampleinfo$AgeAntigen2)

### Please perform kmeans outside the conda environment
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/Y_all_Ags_comp/"
CD8_subset_rld <- readRDS(paste(savedir,"CD8_subset_rld.RDS",sep = ""))
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/"
diff_genes <- read.table(paste(maindir,"pseudobulk/diff_genes_Y_all_Ags",sep = ""))[,1]

sampleinfo <- colData(CD8_subset_rld)
dir.create(savedir, showWarnings = FALSE)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/ATACseq/k_mer_clustering.R")
gap_stat_2 <- kmeans_clustering_ATAC(differential_peaks = diff_genes, rlog_dds = CD8_subset_rld,
                                     sampleinfo = sampleinfo, Group = "AgeAntigen2", saveDir = savedir)

#### We are performing without removing the batch effect also
# savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/no_effect_remove"
# dir.create(savedir,showWarnings=FALSE)
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/ATACseq/k_mer_clustering.R")
# kmeans_clustering_ATAC(differential_peaks = diff_genes, rlog_dds = CD8_subset_rld,
#                        sampleinfo = sampleinfo, Group = "AgeAntigen2", saveDir = savedir)

### Now plotting it
# gap_stat_2 <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/k_means/Differential/cluster_gap.RDS")
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/ATACseq/k_means_plots.R")
dir.create(paste(savedir,"Clustering",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"Clustering/cluster_number.pdf",sep = ""), width = 7, height = 7)
print(fviz_gap_stat(gap_stat_2))
dev.off()

clusters_num <- c(3,5,7)
for (i in 1:length(clusters_num)) {
  df_scale <- read.table(paste(savedir,"Differential/median_counts_for_Group_significant_normalized.txt",sep = ""))
  set.seed(123)
  final_2 <- kmeans(df_scale, clusters_num[i], nstart = 25) # based on the optimal clusters
  print(final_2)
  
  saveRDS(final_2,paste(savedir,"Clustering/final_clusters_",clusters_num[i],".Rds",sep = ""))
  
  peaks_arranged <- as.data.frame(final_2$cluster[order(as.data.frame(final_2$cluster))])
  colnames(peaks_arranged) <- "clusters"
  peaks_arranged$clusters <- paste("C",peaks_arranged$clusters,sep="")
  print(table(peaks_arranged$clusters))
  
  
  # C1  C2  C3  C4  C5
  # 322 246 269 289 293
  
  df_scale_ordered <- df_scale[match(rownames(peaks_arranged), rownames(df_scale)),]
  
  library(ComplexHeatmap)
  library(circlize)
  h=Heatmap(df_scale_ordered,
            name = "scaled", #title of legend
            column_title = "Groups", row_title = "Genes",
            row_names_gp = gpar(fontsize = 0.0001), # Text size for row names
            column_names_gp = gpar(fontsize = 10),
            cluster_columns = TRUE,
            cluster_rows = FALSE,
            col=colorRamp2(c(-0.5, 0, 0.5), 
                           c("blue", "white", "red"))
            #rect_gp = gpar(col = "white", lwd = 1),
            # right_annotation = row_anno,
            # top_annotation = column_anno
            #column_split = factor(c(rep("A", 9), rep("B",8)), levels = c("A","B")),
  )
  
  pdf(paste(savedir,"Clustering/heatmap_",clusters_num[i],"_clusters_change_scale_col_clustered.pdf",sep = ""))
  print(h)
  dev.off()
  
  ### Ordering the columns
  # col_order <- c("Y.BMLF1", "O.BMLF1", "Y.BRLF1", "O.BRLF1", "Y.BMRF1", "O.BMRF1", "Y.EBNA1", "O.EBNA1", "Y.EBNA3C", "O.EBNA3C", "Y.LMP1", "O.LMP1", "Y.LMP2", "O.LMP2", "Y.IE62", "O.IE62")
  col_order <- grep("Y",col_order,value=TRUE)
  df_scale_ordered <- df_scale[match(rownames(peaks_arranged), rownames(df_scale)),match(col_order, colnames(df_scale))]
  
  library(ComplexHeatmap)
  library(circlize)
  h=Heatmap(df_scale_ordered,
            name = "scaled", #title of legend
            column_title = "Groups", row_title = "Genes",
            row_names_gp = gpar(fontsize = 0.0001), # Text size for row names
            column_names_gp = gpar(fontsize = 10),
            cluster_columns = FALSE,
            cluster_rows = FALSE,
            col=colorRamp2(c(-0.5, 0, 0.5), 
                           c("blue", "white", "red"))
            #rect_gp = gpar(col = "white", lwd = 1),
            # right_annotation = row_anno,
            # top_annotation = column_anno
            #column_split = factor(c(rep("A", 9), rep("B",8)), levels = c("A","B")),
  )
  
  pdf(paste(savedir,"Clustering/heatmap_",clusters_num[i],"_clusters_change_scale.pdf",sep = ""))
  print(h)
  dev.off()
  
  
  write.table(peaks_arranged, paste(savedir,"genes_order_",clusters_num[i],".txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
}

### Making an Upset plot ####
dir_name <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Range_10_rem_TRB2_rem_batch_effect/pseudobulk/"

diff_files <- list.files(dir_name, recursive = TRUE, full.names = TRUE, pattern = "O_VS_Y_sig_up.txt")
Ags <- gsub(paste(dir_name,"/clus_all_O.",sep = ""),"",diff_files) %>%
  gsub("/Differential/O_VS_Y_sig_up.txt","",.) %>% gsub("_vs_.*.","",.)

for (i in 1:length(diff_files)) {
  spec_Ag <- rownames(read.table(diff_files[i]))
  assign(Ags[i],spec_Ag)
}

library(data.table)
BMLF1 <- get(Ags[1])
BRLF1 <- get(Ags[2])
EBNA1 <- get(Ags[3])
EBNA3C <- get(Ags[4])
IE62 <- get(Ags[5])
LMP1 <- get(Ags[6])
LMP2 <- get(Ags[7])

n <- max(length(BMLF1), length(BRLF1), length(EBNA1), length(EBNA3C), length(IE62), length(LMP1), length(LMP2))

length(BMLF1) = n                    
length(BRLF1) = n
length(EBNA1) = n
length(EBNA3C) = n
length(IE62) = n
length(LMP1) = n
length(LMP2) = n 

df = as.data.frame(cbind(BMLF1, BRLF1, EBNA1, EBNA3C, IE62, LMP1, LMP2))

### Can be manually be provided the dataframe rather than running it. Should not be blank instead of that it should be NA

df$V1 <- NULL
fld <- fromList(as.list(df)) ## this is an intersect function

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Range_10_rem_TRB2_rem_batch_effect/pseudobulk/"
pdf(paste(savedir,"Upregulated_upset_plot.pdf",sep = ""), width = 10, height = 7)
upset(fld, nsets = 8, text.scale = 1.5)
dev.off()

### To find the interseected genes
list_filter <- list("BMLF1" = BMLF1, "BRLF1" = BRLF1, "EBNA1"= EBNA1, 
                    "EBNA3C" = EBNA3C, "IE62" = IE62, "LMP1" = LMP1, "LMP2" = LMP2)

df2 <- data.frame(gene=unique(unlist(list_filter)))

df1 <- lapply(list_filter,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

library(dplyr)
final <- df_int %>%
  group_by(int) %>%
  dplyr::summarise(n = n()) %>%
  arrange(desc(n))

write.table(df_int, paste(savedir,"Upregulated_gene_intersect.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir,"Upregulated_gene_number.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

### Downregulated Genes
diff_files <- list.files("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Range_10_rem_TRB2_rem_batch_effect/pseudobulk/",
                         recursive = TRUE, full.names = TRUE, pattern = "O_VS_Y_sig_down.txt")

Ags <- gsub("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Range_10_rem_TRB2_rem_batch_effect/pseudobulk//clus_all_O.","",diff_files) %>% gsub("/Differential/O_VS_Y_sig_down.txt","",.) %>% gsub("_vs_.*.","",.)

for (i in 1:length(diff_files)) {
  spec_Ag <- rownames(read.table(diff_files[i]))
  assign(Ags[i],spec_Ag)
}

library(data.table)
library(UpSetR)
BMLF1 <- get(Ags[1])
BRLF1 <- get(Ags[2])
EBNA1 <- get(Ags[3])
EBNA3C <- get(Ags[4])
IE62 <- get(Ags[5])
LMP1 <- get(Ags[6])
LMP2 <- get(Ags[7])

n <- max(length(BMLF1), length(BRLF1), length(EBNA1),length(EBNA3C), length(IE62), length(LMP1), length(LMP2))

length(BMLF1) = n
length(BRLF1) = n
length(EBNA1) = n
length(EBNA3C) = n
length(IE62) = n
length(LMP1) = n
length(LMP2) = n 

df = as.data.frame(cbind(BMLF1, BRLF1, EBNA1, EBNA3C, IE62, LMP1, LMP2))

df$V1 <- NULL

fld <- fromList(as.list(df))

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Range_10_rem_TRB2_rem_batch_effect/pseudobulk/"
pdf(paste(savedir,"downregulated_upset_plot.pdf",sep = ""), width = 10, height = 7)
upset(fld, nsets = 8, text.scale = 1.5)
dev.off()

### To find the interseected genes
list_filter <- list("BMLF1" = BMLF1, "BRLF1" = BRLF1, "EBNA1"= EBNA1, 
                    "EBNA3C" = EBNA3C, "IE62" = IE62, "LMP1" = LMP1, "LMP2" = LMP2)

df2 <- data.frame(gene=unique(unlist(list_filter)))

df1 <- lapply(list_filter,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")


df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

library(dplyr)
final <- df_int %>%
  group_by(int) %>%
  dplyr::summarise(n = n()) %>%
  arrange(desc(n))

write.table(df_int, paste(savedir,"downregulated_gene_intersect.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(final, paste(savedir,"downregulated_gene_number.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

### Absolute PCA and Lutian Analysis ####
## Mem
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/mem/"
dir.create(savedir, showWarnings = FALSE)
Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/mem/Ag_mem_table.txt", header = TRUE, sep = "\t")

rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Samples)
Ag_mem <- Ag_mem[,-1]

metadata <- data.frame(samples = rownames(Ag_mem))
rownames(metadata) <- metadata$samples
metadata$Ags <- gsub(".*.Exp0[1-3]_","",metadata$samples)
metadata$Age <- gsub("_[0-9][0-9]_.*.","",metadata$samples) %>% gsub(".*._EBV_","",.)
metadata$id <- gsub("_.*.","",metadata$samples)
metadata$Run <- gsub(".*._[M|F]_","",metadata$samples) %>% gsub("_.*.","",.)
metadata$gender <- gsub(".*._M_.*.","M",metadata$samples) %>% gsub(".*._F_.*.","F",.)

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
all(metadata$samples == rownames(Ag_mem)) # if TRUE move forward
# Th17_metadata$samplename <- gsub("gE_|DMSO_|_DMSO|_gE|B22015","",Th17_metadata$samples)
## Run will be treated as a covariate in the regression model,
Ag_mem_t <- t(Ag_mem)
pc <- prcomp(Ag_mem_t)
attributes(pc)
pc$center
summary(pc)

pc$rotation

rownames(metadata) <- metadata$samples
pc_with_metadata <- merge(pc$rotation,metadata,by = 'row.names', all = TRUE)

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/mem/"
dir.create(savedir,showWarnings = FALSE)
write.table(pc_with_metadata, file = paste(savedir,"PCA_coordinate_matched_samples_with_metadata.txt",sep = ""),
            quote = F, col.names = T, row.names = F, sep = "\t")

# pc_with_metadata <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/all_run_analysis_CR7/PCA/PCA/absolute_PCA_coordinate_with_metadata.txt",
#                                header = TRUE)
# pc_with_metadata$vaccine <- gsub("S$","Shingrix",pc_with_metadata$vaccine)
# pc_with_metadata$vaccine <- gsub("Z$","Zostavax",pc_with_metadata$vaccine)

pc_summary <- summary(pc)
pc_summary_imp <- as.data.frame(pc_summary$importance)
pc_with_metadata
pc_with_metadata$Ags <- gsub("EBV_|VZV_","",pc_with_metadata$Ags)
pc_with_metadata$id_Ags_gender_Age <- paste(pc_with_metadata$id,pc_with_metadata$Ags,pc_with_metadata$gender, pc_with_metadata$Age,sep = "_")
PC1 <-ggplot(pc_with_metadata, aes(PC1, PC2, label=id_Ags_gender_Age, color=Ags, shape=Age)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(pc_summary_imp[2,1])*100,"% variance")) +
  ylab(paste0("PC2: ",(pc_summary_imp[2,2])*100,"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir,"PCA_norm_samples.pdf",sep = ""),width = 9, height = 8)
PC1
dev.off()


library(umap)
set.seed(2021)
# vsd_cnts <- as.data.frame(vsd@assays@data[[1]])
cnts_forDESeq_umap <- umap(Ag_mem) # others are default
cnts_forDESeq_umap_layout <- as.data.frame(cnts_forDESeq_umap$layout)
colnames(cnts_forDESeq_umap_layout) <- c("UMAP1","UMAP2")
umap_with_metadata <- merge(cnts_forDESeq_umap_layout,metadata,by = 'row.names', all = TRUE)
umap_with_metadata$Ags <- gsub("EBV_|VZV_","",umap_with_metadata$Ags)

umap_with_metadata$id_Ags_gender_Age <- paste(umap_with_metadata$id,umap_with_metadata$Ags,umap_with_metadata$gender, umap_with_metadata$Age,sep = "_")
UMAP1 <-ggplot(umap_with_metadata, aes(UMAP1, UMAP2, label=id_Ags_gender_Age, color=Ags, shape=Age)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("UMAP1")) +
  ylab(paste0("UMAP2")) +
  ggtitle("UMAP") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir,"UMAP1_norm_samples.pdf",sep = ""),width = 9, height = 8)
UMAP1
dev.off()

write.table(umap_with_metadata, paste(savedir,"UMAP_norm_samples.txt",sep = ""), sep = "\t", row.names = T, col.names = T, quote = F)


Ag_mem_Y <- Ag_mem[grep("Y",rownames(Ag_mem)),]
cnts_forDESeq_umap <- umap(Ag_mem_Y) # others are default
cnts_forDESeq_umap_layout <- as.data.frame(cnts_forDESeq_umap$layout)
colnames(cnts_forDESeq_umap_layout) <- c("UMAP1","UMAP2")
umap_with_metadata <- merge(cnts_forDESeq_umap_layout,metadata,by = 'row.names', all = TRUE)
umap_with_metadata$Ags <- gsub("EBV_|VZV_","",umap_with_metadata$Ags)

umap_with_metadata$id_Ags_gender_Age <- paste(umap_with_metadata$id,umap_with_metadata$Ags,umap_with_metadata$gender, umap_with_metadata$Age,sep = "_")
UMAP1 <-ggplot(umap_with_metadata, aes(UMAP1, UMAP2, label=id_Ags_gender_Age, color=Ags, shape=Age)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("UMAP1")) +
  ylab(paste0("UMAP2")) +
  ggtitle("UMAP") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir,"UMAP1_norm_samples_Y.pdf",sep = ""),width = 9, height = 8)
UMAP1
dev.off()

write.table(umap_with_metadata, paste(savedir,"UMAP_norm_samples_Y.txt",sep = ""), sep = "\t", row.names = T, col.names = T, quote = F)

### Lutian Analysis
Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/mem/Ag_mem_table.txt", header = TRUE, sep = "\t")
rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Samples)
Ag_mem <- Ag_mem[,-1]
as.data.frame(print(rowSums(Ag_mem))) ### Sanity Check

metadata <- data.frame(samples = rownames(Ag_mem))
rownames(metadata) <- metadata$samples
metadata$Ags <- gsub(".*.Exp0[1-3]_","",metadata$samples)
metadata$Age <- gsub("_[0-9][0-9]_.*.","",metadata$samples) %>% gsub(".*._EBV_","",.)
metadata$id <- gsub("_.*.","",metadata$samples)
metadata$Run <- gsub(".*._[M|F]_","",metadata$samples) %>% gsub("_.*.","",.)
metadata$gender <- gsub(".*._M_.*.","M",metadata$samples) %>% gsub(".*._F_.*.","F",.)

Ags <- gsub(".*.Exp0[1-3]_","",metadata$samples)
Age <- gsub("_[0-9][0-9]_.*.","",metadata$samples) %>% gsub(".*._EBV_","",.)
id <- gsub("_.*.","",metadata$samples)
Run <- gsub(".*._[M|F]_","",metadata$samples) %>% gsub("_.*.","",.)
Sex <- gsub(".*._M_.*.","M",metadata$samples) %>% gsub(".*._F_.*.","F",.)
Age <- as.numeric(factor(Age))
# age=data[,time[j]]
female=1*(Sex=="F")
# infection=(data$Infection=="Yes")
datax = Ag_mem
print(dim(datax))
  
clus <- dim(datax)[2]
z=rep(NA,clus)
for(k in 1:clus)
  {
    fit=summary(lm(datax[,k]~Age+female))$coef
    z[k]=fit[2,3]
  }
  z.obs=sum(z^2)
  
  B=10000
  z.null=rep(NA, B)
  n1=sum(female==1)
  n0=sum(female==0)
  for(b in 1:B){
    Age0=Age
    Age0[female==1]=sample(Age[female==1], n1, F)
    Age0[female==0]=sample(Age[female==0], n0, F)
    z=rep(NA,clus)
    for(k in 1:clus)
    {
      fit=summary(lm(datax[,k]~Age0+female))$coef
      z[k]=fit[2,3]
    }
    z.null[b]=sum(z^2)
  }
  # print(paste("Removed this vaccine",vaccine[i],"for this time",time[j],sep = " "))
  # print(paste("Removed",vac[i],sep = " "))
print(mean(z.null>=z.obs))

pvalue <- mean(z.null>=z.obs)
write.table(pvalue, paste(savedir,"Lutian_pvalue.txt",sep = ""),quote = F, row.names = F, col.names = F, sep = "\t")

#### Performing for each Antigen
Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/mem/Ag_mem_table.txt", header = TRUE, sep = "\t")
# Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/mem/Ag_mem_table2.txt", header = TRUE, sep = "\t")
rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Samples)
Ag_mem <- Ag_mem[,-1]
as.data.frame(print(rowSums(Ag_mem))) ### Sanity Check

Ag_unique <- unique(gsub(".*.Exp0[0-9]_","",rownames(Ag_mem)))
df = as.data.frame(matrix(, nrow = length(Ag_unique), ncol = 1))
rownames(df) <- Ag_unique
colnames(df) <- "pvalue"

for (i in 1:length(Ag_unique)) {
  Ag_mem_subset <- Ag_mem[grep(Ag_unique[i],rownames(Ag_mem)),]
  metadata <- data.frame(samples = rownames(Ag_mem_subset))
  rownames(metadata) <- metadata$samples
  metadata$Ags <- gsub(".*.Exp0[1-3]_","",metadata$samples)
  metadata$Age <- gsub("_[0-9][0-9]_.*.","",metadata$samples) %>% gsub(".*._EBV_","",.)
  metadata$id <- gsub("_.*.","",metadata$samples)
  metadata$Run <- gsub(".*._[M|F]_","",metadata$samples) %>% gsub("_.*.","",.)
  metadata$gender <- gsub(".*._M_.*.","M",metadata$samples) %>% gsub(".*._F_.*.","F",.)
  
  Ags <- gsub(".*.Exp0[1-3]_","",metadata$samples)
  Age <- gsub("_[0-9][0-9]_.*.","",metadata$samples) %>% gsub(".*._EBV_","",.)
  id <- gsub("_.*.","",metadata$samples)
  Run <- gsub(".*._[M|F]_","",metadata$samples) %>% gsub("_.*.","",.)
  Sex <- gsub(".*._M_.*.","M",metadata$samples) %>% gsub(".*._F_.*.","F",.)
  Age <- as.numeric(factor(Age))
  # age=data[,time[j]]
  female=1*(Sex=="F")
  # infection=(data$Infection=="Yes")
  datax = Ag_mem_subset
  print(dim(datax))
  
  clus <- dim(datax)[2]
  z=rep(NA,clus)
  for(k in 1:clus)
  {
    fit=summary(lm(datax[,k]~Age+female))$coef
    z[k]=fit[2,3]
  }
  z <- z[!is.na(z)]
  z.obs=sum(z^2)
  
  B=25000
  z.null=rep(NA, B)
  n1=sum(female==1)
  n0=sum(female==0)
  for(b in 1:B){
    Age0=Age
    Age0[female==1]=sample(Age[female==1], n1, F)
    Age0[female==0]=sample(Age[female==0], n0, F)
    z=rep(NA,clus)
    for(k in 1:clus)
    {
      fit=summary(lm(datax[,k]~Age0+female))$coef
      z[k]=fit[2,3]
    }
    z <- z[!is.na(z)]
    z.null[b]=sum(z^2)
  }
  # print(paste("Removed this vaccine",vaccine[i],"for this time",time[j],sep = " "))
  # print(paste("Removed",vac[i],sep = " "))
  print(mean(z.null>=z.obs))
  
  df[i,] <- mean(z.null>=z.obs)
}

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/mem/"
write.table(df, paste(savedir,"Lutian_pvalue_no_clus17_25k_iterations.txt",sep = ""),quote = F, row.names = T, col.names = T, sep = "\t")

### Only Young
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/mem/"
dir.create(savedir, showWarnings = FALSE)
Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/mem/Ag_mem_table.txt", header = TRUE, sep = "\t")

rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Samples)
Ag_mem <- Ag_mem[,-1]

Ag_mem_Y <- Ag_mem[grep("^Y",rownames(Ag_mem)),]

metadata <- data.frame(samples = rownames(Ag_mem_Y))
rownames(metadata) <- metadata$samples
metadata$Ags <- gsub(".*.Exp0[1-3]_","",metadata$samples)
metadata$Age <- gsub("_[0-9][0-9]_.*.","",metadata$samples) %>% gsub(".*._EBV_","",.)
metadata$id <- gsub("_.*.","",metadata$samples)
metadata$Run <- gsub(".*._[M|F]_","",metadata$samples) %>% gsub("_.*.","",.)
metadata$gender <- gsub(".*._M_.*.","M",metadata$samples) %>% gsub(".*._F_.*.","F",.)

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

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/mem/"
dir.create(savedir,showWarnings = FALSE)
write.table(pc_with_metadata, file = paste(savedir,"PCA_coordinate_matched_samples_with_metadata_Young.txt",sep = ""),
            quote = F, col.names = T, row.names = F, sep = "\t")

# pc_with_metadata <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/all_run_analysis_CR7/PCA/PCA/absolute_PCA_coordinate_with_metadata.txt",
#                                header = TRUE)
# pc_with_metadata$vaccine <- gsub("S$","Shingrix",pc_with_metadata$vaccine)
# pc_with_metadata$vaccine <- gsub("Z$","Zostavax",pc_with_metadata$vaccine)

pc_summary <- summary(pc)
pc_summary_imp <- as.data.frame(pc_summary$importance)
pc_with_metadata
pc_with_metadata$Ags <- gsub("EBV_|VZV_","",pc_with_metadata$Ags)
pc_with_metadata$id_Ags_gender_Age <- paste(pc_with_metadata$id,pc_with_metadata$Ags,pc_with_metadata$gender, pc_with_metadata$Age,sep = "_")
PC1 <-ggplot(pc_with_metadata, aes(PC1, PC2, label=id_Ags_gender_Age, color=Ags, shape=Age)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(pc_summary_imp[2,1])*100,"% variance")) +
  ylab(paste0("PC2: ",(pc_summary_imp[2,2])*100,"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir,"PCA_Y_norm_samples.pdf",sep = ""),width = 9, height = 8)
PC1
dev.off()

### Antigen Lutian Analysis
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/mem/"
dir.create(savedir, showWarnings = FALSE)
Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/mem/Ag_mem_table.txt", header = TRUE, sep = "\t")

rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Samples)
Ag_mem <- Ag_mem[,-1]

Ag_mem_Y <- Ag_mem[grep("^Y",rownames(Ag_mem)),]
Ag <- unique(gsub(".*._EBV_|*.*_VZV_","",rownames(Ag_mem_Y)))

df <- data.frame(matrix(nrow = length(Ag), ncol = length(Ag)))
rownames(df) <- Ag
colnames(df) <- Ag

for (i in 1:(length(Ag)-1)) {
  Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/mem/Ag_mem_table.txt", header = TRUE, sep = "\t")
  rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Samples)
  Ag_mem <- Ag_mem[,-1]
  Ag_mem_Y <- Ag_mem[grep("^Y",rownames(Ag_mem)),]
  
  for (j in 1:length(Ag)) {
    print(paste("Performing Lutian analysis for ",Ag[i]," and ",Ag[j]))
    data <- Ag_mem_Y[grep(paste(Ag[i],"|",Ag[j],sep = ""),rownames(Ag_mem_Y)),]
    
    if (Ag[i] == Ag[j]) {
      print(paste(Ag[i]," is same as ",Ag[j]))
    } else if (Ag[i] != Ag[j]){
      Sex <- gsub(".*._M_.*","M",rownames(data)) %>% gsub(".*._F_.*.","F",.)
      Ag_var <- gsub(".*._EBV_|.*._VZV_","",rownames(data))
      Ag_var <- as.numeric(factor(Ag_var))
      female=1*(Sex=="F")
      
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

Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/mem/Ag_mem_table.txt", header = TRUE, sep = "\t")
rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Samples)
Ag_mem <- Ag_mem[,-1]
Ag_mem_Y <- Ag_mem[grep("^Y",rownames(Ag_mem)),]

data <- Ag_mem_Y[grep("VZV",rownames(Ag_mem_Y),invert = T),]
# print(rownames(data))
Sex <- gsub(".*._M_.*","M",rownames(data)) %>% gsub(".*._F_.*.","F",.)
Ag_var <- gsub(".*._BMLF1|.*._BMRF1|.*._BRLF1","Lyt",rownames(data)) %>% gsub(".*._EBNA1|.*._EBNA3C|.*._LMP1|.*._LMP2","Lat",.)
Ag_var <- as.numeric(factor(Ag_var))
print(Ag_var)
female=1*(Sex=="F")

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
pvalue <- mean(z.null>=z.obs)

write.table(pvalue, paste(savedir,"Lutian_Y_only_Ag_Lyt_vs_Lat.txt",sep = ""), sep = "\t", quote = F, row.names = T, col.names = T)




## All Clusters
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/all_clus/"
dir.create(savedir, showWarnings = FALSE)
Ag_mem <- read.table(paste(savedir,"Ag_table.txt",sep = ""), header = TRUE, sep = "\t")

rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Sample)
Ag_mem <- Ag_mem[,-1]

metadata <- data.frame(samples = rownames(Ag_mem))
rownames(metadata) <- metadata$samples
metadata$Ags <- gsub(".*.Exp0[1-3]_","",metadata$samples)
metadata$Age <- gsub("_[0-9][0-9]_.*.","",metadata$samples) %>% gsub(".*._EBV_","",.)
metadata$id <- gsub("_.*.","",metadata$samples)
metadata$Run <- gsub(".*._[M|F]_","",metadata$samples) %>% gsub("_.*.","",.)
metadata$gender <- gsub(".*._M_.*.","M",metadata$samples) %>% gsub(".*._F_.*.","F",.)

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
all(metadata$samples == rownames(Ag_mem)) # if TRUE move forward
# Th17_metadata$samplename <- gsub("gE_|DMSO_|_DMSO|_gE|B22015","",Th17_metadata$samples)
## Run will be treated as a covariate in the regression model,
Ag_mem_t <- t(Ag_mem)
pc <- prcomp(Ag_mem_t)
attributes(pc)
pc$center
summary(pc)

pc$rotation

rownames(metadata) <- metadata$samples
pc_with_metadata <- merge(pc$rotation,metadata,by = 'row.names', all = TRUE)

write.table(pc_with_metadata, file = paste(savedir,"PCA_coordinate_matched_samples_with_metadata.txt",sep = ""),
            quote = F, col.names = T, row.names = F, sep = "\t")

pc_summary <- summary(pc)
pc_summary_imp <- as.data.frame(pc_summary$importance)
pc_with_metadata
pc_with_metadata$Ags <- gsub("EBV_|VZV_","",pc_with_metadata$Ags)
pc_with_metadata$id_Ags_gender_Age <- paste(pc_with_metadata$id,pc_with_metadata$Ags,pc_with_metadata$gender, pc_with_metadata$Age,sep = "_")
PC1 <-ggplot(pc_with_metadata, aes(PC1, PC2, label=id_Ags_gender_Age, color=Ags, shape=Age)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(pc_summary_imp[2,1])*100,"% variance")) +
  ylab(paste0("PC2: ",(pc_summary_imp[2,2])*100,"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir,"PCA_norm_samples.pdf",sep = ""),width = 9, height = 8)
PC1
dev.off()


### UMAP
library(umap)
set.seed(2021)
# vsd_cnts <- as.data.frame(vsd@assays@data[[1]])
cnts_forDESeq_umap <- umap(Ag_mem) # others are default
cnts_forDESeq_umap_layout <- as.data.frame(cnts_forDESeq_umap$layout)
colnames(cnts_forDESeq_umap_layout) <- c("UMAP1","UMAP2")
umap_with_metadata <- merge(cnts_forDESeq_umap_layout,metadata,by = 'row.names', all = TRUE)
umap_with_metadata$Ags <- gsub("EBV_|VZV_","",umap_with_metadata$Ags)

umap_with_metadata$id_Ags_gender_Age <- paste(umap_with_metadata$id,umap_with_metadata$Ags,umap_with_metadata$gender, umap_with_metadata$Age,sep = "_")
UMAP1 <-ggplot(umap_with_metadata, aes(UMAP1, UMAP2, label=id_Ags_gender_Age, color=Ags, shape=Age)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("UMAP1")) +
  ylab(paste0("UMAP2")) +
  ggtitle("UMAP") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir,"UMAP1_norm_samples.pdf",sep = ""),width = 9, height = 8)
UMAP1
dev.off()

write.table(umap_with_metadata, paste(savedir,"UMAP_norm_samples.txt",sep = ""), sep = "\t", row.names = T, col.names = T, quote = F)


Ag_mem_Y <- Ag_mem[grep("Y",rownames(Ag_mem)),]
cnts_forDESeq_umap <- umap(Ag_mem_Y) # others are default
cnts_forDESeq_umap_layout <- as.data.frame(cnts_forDESeq_umap$layout)
colnames(cnts_forDESeq_umap_layout) <- c("UMAP1","UMAP2")
umap_with_metadata <- merge(cnts_forDESeq_umap_layout,metadata,by = 'row.names', all = TRUE)
umap_with_metadata$Ags <- gsub("EBV_|VZV_","",umap_with_metadata$Ags)

umap_with_metadata$id_Ags_gender_Age <- paste(umap_with_metadata$id,umap_with_metadata$Ags,umap_with_metadata$gender, umap_with_metadata$Age,sep = "_")
UMAP1 <-ggplot(umap_with_metadata, aes(UMAP1, UMAP2, label=id_Ags_gender_Age, color=Ags, shape=Age)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("UMAP1")) +
  ylab(paste0("UMAP2")) +
  ggtitle("UMAP") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir,"UMAP1_norm_samples_Y.pdf",sep = ""),width = 9, height = 8)
UMAP1
dev.off()

write.table(umap_with_metadata, paste(savedir,"UMAP_norm_samples_Y.txt",sep = ""), sep = "\t", row.names = T, col.names = T, quote = F)

### Lutian Analysis
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/all_clus/"
Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/all_clus/Ag_table.txt", header = TRUE, sep = "\t")
rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Sample)
Ag_mem <- Ag_mem[,-1]
as.data.frame(print(rowSums(Ag_mem))) ### Sanity Check

Ags <- gsub(".*.Exp0[1-3]_","",metadata$samples)
Age <- gsub("_[0-9][0-9]_.*.","",metadata$samples) %>% gsub(".*._EBV_","",.)
id <- gsub("_.*.","",metadata$samples)
Run <- gsub(".*._[M|F]_","",metadata$samples) %>% gsub("_.*.","",.)
Sex <- gsub(".*._M_.*.","M",metadata$samples) %>% gsub(".*._F_.*.","F",.)
Age <- as.numeric(factor(Age))
# age=data[,time[j]]
female=1*(Sex=="F")
# infection=(data$Infection=="Yes")
datax = Ag_mem
print(dim(datax))

clus <- dim(datax)[2]
z=rep(NA,clus)
for(k in 1:clus){
  fit=summary(lm(datax[,k]~Age+female))$coef
  z[k]=fit[2,3]
}
z.obs=sum(z^2)

B=10000
z.null=rep(NA, B)
n1=sum(female==1)
n0=sum(female==0)
for(b in 1:B){
  Age0=Age
  Age0[female==1]=sample(Age[female==1], n1, F)
  Age0[female==0]=sample(Age[female==0], n0, F)
  z=rep(NA,clus)
  for(k in 1:clus)
  {
    fit=summary(lm(datax[,k]~Age0+female))$coef
    z[k]=fit[2,3]
  }
  z.null[b]=sum(z^2)
}
# print(paste("Removed this vaccine",vaccine[i],"for this time",time[j],sep = " "))
# print(paste("Removed",vac[i],sep = " "))
print(mean(z.null>=z.obs))
pvalue <- mean(z.null>=z.obs)
write.table(pvalue, paste(savedir,"Lutian_pvalue.txt",sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")

### All Antigen ####
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/all_clus/"
Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/all_clus/Ag_table.txt", header = TRUE, sep = "\t")
rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Sample)
Ag_mem <- Ag_mem[,-1]
# as.data.frame(print(rowSums(Ag_mem))) ### Sanity Check

Ag_unique <- unique(gsub(".*.Exp0[0-9]_","",rownames(Ag_mem)))
df = as.data.frame(matrix(, nrow = length(Ag_unique), ncol = 1))
rownames(df) <- Ag_unique
colnames(df) <- "pvalue"

for (i in 1:length(Ag_unique)) {
  Ag_mem_subset <- Ag_mem[grep(Ag_unique[i],rownames(Ag_mem)),]
  metadata <- data.frame(samples = rownames(Ag_mem_subset))
  rownames(metadata) <- metadata$samples
  metadata$Ags <- gsub(".*.Exp0[1-3]_","",metadata$samples)
  metadata$Age <- gsub("_[0-9][0-9]_.*.","",metadata$samples) %>% gsub(".*._EBV_","",.)
  metadata$id <- gsub("_.*.","",metadata$samples)
  metadata$Run <- gsub(".*._[M|F]_","",metadata$samples) %>% gsub("_.*.","",.)
  metadata$gender <- gsub(".*._M_.*.","M",metadata$samples) %>% gsub(".*._F_.*.","F",.)
  
  Ags <- gsub(".*.Exp0[1-3]_","",metadata$samples)
  Age <- gsub("_[0-9][0-9]_.*.","",metadata$samples) %>% gsub(".*._EBV_","",.)
  id <- gsub("_.*.","",metadata$samples)
  Run <- gsub(".*._[M|F]_","",metadata$samples) %>% gsub("_.*.","",.)
  Sex <- gsub(".*._M_.*.","M",metadata$samples) %>% gsub(".*._F_.*.","F",.)
  Age <- as.numeric(factor(Age))
  # age=data[,time[j]]
  female=1*(Sex=="F")
  # infection=(data$Infection=="Yes")
  datax = Ag_mem_subset
  print(dim(datax))
  
  clus <- dim(datax)[2]
  z=rep(NA,clus)
  for(k in 1:clus)
  {
    fit=summary(lm(datax[,k]~Age+female))$coef
    z[k]=fit[2,3]
  }
  z <- z[!is.na(z)]
  z.obs=sum(z^2)
  
  B=10000
  z.null=rep(NA, B)
  n1=sum(female==1)
  n0=sum(female==0)
  for(b in 1:B){
    Age0=Age
    Age0[female==1]=sample(Age[female==1], n1, F)
    Age0[female==0]=sample(Age[female==0], n0, F)
    z=rep(NA,clus)
    for(k in 1:clus)
    {
      fit=summary(lm(datax[,k]~Age0+female))$coef
      z[k]=fit[2,3]
    }
    z <- z[!is.na(z)]
    z.null[b]=sum(z^2)
  }
  # print(paste("Removed this vaccine",vaccine[i],"for this time",time[j],sep = " "))
  # print(paste("Removed",vac[i],sep = " "))
  print(mean(z.null>=z.obs))
  
  df[i,] <- mean(z.null>=z.obs)
}

write.table(df, paste(savedir,"Lutian_pvalue.txt",sep = ""),quote = F, row.names = T, col.names = T, sep = "\t")

#### Only Young
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/all_clus/"
dir.create(savedir, showWarnings = FALSE)
Ag_mem <- read.table(paste(savedir,"Ag_table.txt",sep = ""), header = TRUE, sep = "\t")

rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Sample)
Ag_mem <- Ag_mem[,-1]

Ag_mem_Y <- Ag_mem[grep("^Y",rownames(Ag_mem)),]

metadata <- data.frame(samples = rownames(Ag_mem_Y))
rownames(metadata) <- metadata$samples
metadata$Ags <- gsub(".*.Exp0[1-3]_","",metadata$samples)
metadata$Age <- gsub("_[0-9][0-9]_.*.","",metadata$samples) %>% gsub(".*._EBV_","",.)
metadata$id <- gsub("_.*.","",metadata$samples)
metadata$Run <- gsub(".*._[M|F]_","",metadata$samples) %>% gsub("_.*.","",.)
metadata$gender <- gsub(".*._M_.*.","M",metadata$samples) %>% gsub(".*._F_.*.","F",.)

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

write.table(pc_with_metadata, file = paste(savedir,"PCA_coordinate_matched_samples_with_metadata_Young.txt",sep = ""),
            quote = F, col.names = T, row.names = F, sep = "\t")

pc_summary <- summary(pc)
pc_summary_imp <- as.data.frame(pc_summary$importance)
pc_with_metadata
pc_with_metadata$Ags <- gsub("EBV_|VZV_","",pc_with_metadata$Ags)
pc_with_metadata$id_Ags_gender_Age <- paste(pc_with_metadata$id,pc_with_metadata$Ags,pc_with_metadata$gender, pc_with_metadata$Age,sep = "_")
PC1 <-ggplot(pc_with_metadata, aes(PC1, PC2, label=id_Ags_gender_Age, color=Ags, shape=Age)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(pc_summary_imp[2,1])*100,"% variance")) +
  ylab(paste0("PC2: ",(pc_summary_imp[2,2])*100,"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir,"PCA_Young_norm_samples.pdf",sep = ""),width = 9, height = 8)
PC1
dev.off()

### Antigen Lutian Analysis
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/all_clus/"
dir.create(savedir, showWarnings = FALSE)
Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/all_clus/Ag_table.txt", header = TRUE, sep = "\t")

rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Sample)
Ag_mem <- Ag_mem[,-1]

Ag_mem_Y <- Ag_mem[grep("^Y",rownames(Ag_mem)),]
Ag <- unique(gsub(".*._EBV_|*.*_VZV_","",rownames(Ag_mem_Y)))

df <- data.frame(matrix(nrow = length(Ag), ncol = length(Ag)))
rownames(df) <- Ag
colnames(df) <- Ag

for (i in 1:(length(Ag)-1)) {
  Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/all_clus/Ag_table.txt", header = TRUE, sep = "\t")
  rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Sample)
  Ag_mem <- Ag_mem[,-1]
  Ag_mem_Y <- Ag_mem[grep("^Y",rownames(Ag_mem)),]
  
  for (j in 1:length(Ag)) {
    print(paste("Performing Lutian analysis for ",Ag[i]," and ",Ag[j]))
    data <- Ag_mem_Y[grep(paste(Ag[i],"|",Ag[j],sep = ""),rownames(Ag_mem_Y)),]
    
    if (Ag[i] == Ag[j]) {
      print(paste(Ag[i]," is same as ",Ag[j]))
    } else if (Ag[i] != Ag[j]){
      print(rownames(data))
      Sex <- gsub(".*._M_.*","M",rownames(data)) %>% gsub(".*._F_.*.","F",.)
      Ag_var <- gsub(".*._EBV_|.*._VZV_","",rownames(data))
      Ag_var <- as.numeric(factor(Ag_var))
      print(Ag_var)
      female=1*(Sex=="F")
      
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

write.table(df, paste(savedir,"Lutian_Y_only_Ag_all_cluster.txt",sep = ""), sep = "\t", quote = F, row.names = T, col.names = T)


### Lytic VS Latent
Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Lutian_analysis/all_clus/Ag_table.txt", header = TRUE, sep = "\t")
rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Sample)
Ag_mem <- Ag_mem[,-1]
Ag_mem_Y <- Ag_mem[grep("^Y",rownames(Ag_mem)),]

data <- Ag_mem_Y[grep("VZV",rownames(Ag_mem_Y),invert = T),]
# print(rownames(data))
Sex <- gsub(".*._M_.*","M",rownames(data)) %>% gsub(".*._F_.*.","F",.)
Ag_var <- gsub(".*._BMLF1|.*._BMRF1|.*._BRLF1","Lyt",rownames(data)) %>% gsub(".*._EBNA1|.*._EBNA3C|.*._LMP1|.*._LMP2","Lat",.)
Ag_var <- as.numeric(factor(Ag_var))
print(Ag_var)
female=1*(Sex=="F")

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
pvalue <- mean(z.null>=z.obs)




write.table(df, paste(savedir,"Lutian_Y_only_Ag_Lyt_vs_Lat.txt",sep = ""), sep = "\t", quote = F, row.names = T, col.names = T)


