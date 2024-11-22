### We have process the data /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/Ines_BEAM_T.R
library(Seurat)
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/"
CD8_bulk <- readRDS(paste(savedir,"saveRDS_obj/CD8_modality_integrate_cluster_0.6.RDS",sep = ""))

memory_cells <- rownames(CD8_bulk@meta.data[grep("naïve",CD8_bulk@meta.data$celltypes,invert=TRUE),])
CD8_mem <- subset(CD8_bulk, cells = memory_cells)

### Performing Differential
group1 = "O"
group2 = "Y"
remove_samples = NULL
remove_samples_name <- paste(remove_samples,collapse="_and_")
DefaultAssay(CD8_mem) <- "RNA"
CD8_mem
# savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run01_to_10_analysis/vaccine/S_vs_Z/"
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/COVID_pseudobulk_PCA_within_AJ_2.R")
CD8_subset_dds <- COVID_pseudobulk_within_cluster_AJ(obj = CD8_mem, savedir = savedir, group1 = group1, group2 = group2,
                                                     grouping_by = "Age", cluster = "all", cell_freq = 20, remove_samples = remove_samples,
                                                     cluster_group = "seurat_clusters", sample_col = "orig.ident", batch_col = "Run",
                                                     gene_min_counts =  5, column_name_split = c("sampleid","virus","age","age_number","gender","Run"))

cluster2 <- paste(cluster, sep="_", collapse="_")
savedir2 <- paste(savedir,"pseudobulk/clus_",cluster2,"_",group1,"_vs_",group2,"/",sep = "")
# savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run04/Analysis/pseudobulk/clus_13_removed__DMSO_vs_Sprotein/"
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/RNAseq/pipeline_function/RNAseq_limma_EdgeR.R")
design0 <- model.matrix(~ 0 + Age + Run, data = colData(CD8_subset_dds))
colnames(design0) <- c("O","Y",paste("Run",1:(ncol(design0)-2),sep = ""))
cm <- makeContrasts(O_VS_Y = O-Y,levels = design0)
desl_clus <- LimmaEdgeR_differential(dds = CD8_subset_dds,
                                     design0 = design0,
                                     cm = cm, 
                                     savedir = savedir2,
                                     logfc = 0.5,
                                     p_value_adj = 0.05)


### Absolute PCA and Lutian Analysis ####
## Mem
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/Lutian_analysis/"
dir.create(savedir, showWarnings = FALSE)
Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/Lutian_analysis/mem/bulk_memory", header = TRUE, sep = "\t")

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
library(Glimma)
library(ggplot2)
library(ggrepel)
# library(session)
# library(DEFormats)
library(EnhancedVolcano)
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

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/Lutian_analysis/mem/"
dir.create(savedir,showWarnings = FALSE)
write.table(pc_with_metadata, file = paste(savedir,"PCA_coordinate_matched_samples_with_metadata.txt",sep = ""),
            quote = F, col.names = T, row.names = F, sep = "\t")

pc_summary <- summary(pc)
pc_summary_imp <- as.data.frame(pc_summary$importance)
pc_with_metadata
pc_with_metadata$Ags <- gsub("EBV_|VZV_","",pc_with_metadata$Ags)
pc_with_metadata$id_gender_Age <- paste(pc_with_metadata$id,pc_with_metadata$gender, pc_with_metadata$Age,sep = "_")
PC1 <-ggplot(pc_with_metadata, aes(PC1, PC2, label=id_gender_Age, color=Age, shape=gender)) +
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

### Lutian Analysis
Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/Lutian_analysis/mem/bulk_memory", header = TRUE, sep = "\t")
rownames(Ag_mem) <- gsub(" ","_",Ag_mem$Sample)
Ag_mem <- Ag_mem[,-1]
as.data.frame(print(rowSums(Ag_mem))) ### Sanity Check

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

## All Clusters
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/Lutian_analysis/all_clus/"
dir.create(savedir, showWarnings = FALSE)
Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/Lutian_analysis/all_clus/bulk_all", header = TRUE, sep = "\t")

rownames(Ag_mem) <- gsub(" ","_",Ag_mem$X)
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
library(Glimma)
library(ggplot2)
library(ggrepel)
# library(session)
# library(DEFormats)
library(EnhancedVolcano)
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

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/Lutian_analysis/all_clus//"
dir.create(savedir,showWarnings = FALSE)
write.table(pc_with_metadata, file = paste(savedir,"PCA_coordinate_matched_samples_with_metadata.txt",sep = ""),
            quote = F, col.names = T, row.names = F, sep = "\t")

pc_summary <- summary(pc)
pc_summary_imp <- as.data.frame(pc_summary$importance)
pc_with_metadata
pc_with_metadata$Ags <- gsub("EBV_|VZV_","",pc_with_metadata$Ags)
pc_with_metadata$id_gender_Age <- paste(pc_with_metadata$id,pc_with_metadata$gender, pc_with_metadata$Age,sep = "_")
PC1 <-ggplot(pc_with_metadata, aes(PC1, PC2, label=id_gender_Age, color=Age, shape=gender)) +
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

### Lutian Analysis
Ag_mem <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/Lutian_analysis/all_clus//bulk_all", header = TRUE, sep = "\t")
rownames(Ag_mem) <- gsub(" ","_",Ag_mem$X)
Ag_mem <- Ag_mem[,-1]
as.data.frame(print(rowSums(Ag_mem))) ### Sanity Check

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
savedir = c("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/GSEA/")
dir.create(savedir, showWarnings = FALSE)

maindir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/pseudobulk/"
files <- grep("Y_vs_O",list.files(maindir, pattern = "_genes_no_filter.txt", full.names = TRUE, recursive = TRUE),invert = TRUE, value = TRUE)
Ag_name <- "Bulk"

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


