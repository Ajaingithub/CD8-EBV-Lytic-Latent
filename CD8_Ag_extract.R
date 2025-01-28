# https://omarwagih.github.io/ggseqlogo/
#### Extracting out the TCR and Ag extraction for uploading to the VDJ database. So creating TCR alpha amino acid, TCR beta amino acid
library(Seurat)
library(dplyr)
# CD8_Ag <- readRDS("/diazlab/data3/.abhinav/.immune/cancer_combined/project/resources/GSE275633_CD8_Antigen_BEAM_T.RDS")
CD8_required <- CD8_Ag@meta.data[, c("TRA1_TRA2_TRB1_TRB2_cdraa", "TRA1_TRA2_TRB1_TRB2_vdjc", "Ag_range_10")]
CD8_required_noNA <- na.omit(CD8_required)
# unique_clones <- unique(CD8_required_noNA$TRA1_TRA2_TRB1_TRB2_cdraa)

# CD8_required_noNA_2 <- CD8_required_noNA[match(unique_clones, CD8_required_noNA$TRA1_TRA2_TRB1_TRB2_cdraa),]

### Extracting out Alpha CDR3 amino acids and V and J alpha
CD8_required_noNA$cdr3.alpha <- gsub("-.*.", "", CD8_required_noNA$TRA1_TRA2_TRB1_TRB2_cdraa) %>%
    gsub("*.*_", "", .)
CD8_required_noNA$v.alpha <- gsub("-[A-Z].*.", "", CD8_required_noNA$TRA1_TRA2_TRB1_TRB2_vdjc) %>%
    gsub("__.*.", "", .)
CD8_required_noNA$j.alpha <- gsub("-[A-Z].*.", "", CD8_required_noNA$TRA1_TRA2_TRB1_TRB2_vdjc) %>%
    gsub(".*.__", "", .) %>%
    gsub("_TRAC", "", .)

# Beta CDR3 amino acid, V, D, J
CD8_required_noNA$cdr3.beta <- gsub(".*.--", "", CD8_required_noNA$TRA1_TRA2_TRB1_TRB2_cdraa) %>%
    gsub("-.*.", "", .) %>%
    gsub(".*._", "", .)
CD8_required_noNA$v.beta <- gsub(".*.--", "", CD8_required_noNA$TRA1_TRA2_TRB1_TRB2_vdjc) %>%
    gsub("_.*.", "", .)
CD8_required_noNA$d.beta <- gsub(".*.--", "", CD8_required_noNA$TRA1_TRA2_TRB1_TRB2_vdjc) %>%
    gsub("_TRBJ.*.", "", .) %>%
    gsub("*.*_", "", .) %>%
    gsub("NA-NA", "", .)
CD8_required_noNA$j.beta <- gsub(".*.--", "", CD8_required_noNA$TRA1_TRA2_TRB1_TRB2_vdjc) %>%
    gsub("_TRBC[1-2]*.*", "", .) %>%
    gsub(".*._", "", .)

#### Making a TCR and individual matrix
TCR_matrix <- table(CD8_Ag@meta.data$CDR3, CD8_Ag@meta.data$samplename)
TCR_matrix2 <- TCR_matrix
TCR_matrix2[TCR_matrix2 > 0] <- 1
shared_TCRs <- TCR_matrix[rowSums(TCR_matrix2) > 1, ]

savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/"
write.table(shared_TCRs, paste0(savedir, "Table/TCR_individual_matrix_shared_TCRs_CDR3_aa.txt"), sep = "\t", row.names = T, col.names = T, quote = F)

TCR_matrix <- table(CD8_Ag@meta.data$TRA1_or_TRA2_TRB1_no_TRB2, CD8_Ag@meta.data$samplename)
TCR_matrix2 <- TCR_matrix
TCR_matrix2[TCR_matrix2 > 0] <- 1
shared_TCRs <- TCR_matrix[rowSums(TCR_matrix2) > 1, ]

savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/"
write.table(shared_TCRs, paste0(savedir, "Table/TCR_individual_matrix_shared_TCRs_nucleotides.txt"), sep = "\t", row.names = T, col.names = T, quote = F)

#### Making the UMAP for the clones which are shared between individuals
top13 <- read.table("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/Ines_selected_top13.txt", header = FALSE)[, 1]

dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE)
rm(plot_list)
plot_list <- list()
as.data.frame(CD8_Ag@reductions$wnn.umap@cell.embeddings) -> wnn_umap_embed
for (i in 8:length(top13)) {
    CD8_Ag@meta.data$same_clone_diff_ind <- "zz"
    index <- grep(top13[i], CD8_Ag@meta.data$TRA1_or_TRA2_TRB1_no_TRB2)
    CD8_Ag@meta.data$same_clone_diff_ind[index] <- CD8_Ag@meta.data$orig.ident[index]
    stopifnot(rownames(CD8_Ag@meta.data) == rownames(wnn_umap_embed))
    wnn_umap_embed$same_clone_diff_ind <- CD8_Ag@meta.data$same_clone_diff_ind
    wnn_umap_embed$same_clone_diff_ind[grep("O5|Y6", wnn_umap_embed$same_clone_diff_ind)] <- "zz"
    wnn_umap_embed_ordered <- wnn_umap_embed[order(wnn_umap_embed$same_clone_diff_ind, decreasing = TRUE), ]

    p <- ggplot(wnn_umap_embed_ordered, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = same_clone_diff_ind, size = same_clone_diff_ind)) +
        geom_point() +
        scale_color_manual(values = c("deeppink", "dodgerblue3", "grey87")) +
        scale_size_manual(values = c(2, 2, 0.1)) +
        theme(
            plot.title = element_text(hjust = 0.5, size = 4),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3),
            panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.minor = element_line(colour = "white"),
            panel.grid.major = element_line(colour = "white")
        ) +
        # NoLegend() +
        ggtitle(top13[i])

    plot_list[[i]] <- p
    pdf(paste0(savedir, "UMAP/clone", i, ".pdf"))
    print(p)
    dev.off()
}

pdf(paste0(savedir, "UMAP/all_sharing_clones.pdf"))
plot_list
dev.off()

### Making a table to identify that each antigen has same epitope
for (i in 1:length(top13)) {
    index <- grep(top13[i], CD8_Ag@meta.data$TRA1_or_TRA2_TRB1_no_TRB2)
    clone_table <- CD8_Ag@meta.data[index, c(
        "orig.ident", "Ag_range_10", "specificity_score_range_10",
        "max_Ag", "max_specificity_score", "EBV_BALF4",
        "EBV_BMLF1", "EBV_BMRF1", "EBV_BRLF1", "EBV_EBNA1", "EBV_EBNA3C",
        "EBV_LMP2", "VZV_IE62", "VZV_IE63", "TRA1_or_TRA2_TRB1_no_TRB2"
    )]
    write.table(clone_table, paste0(savedir, "Table/clone_table", i, ".txt"), quote = F, row.names = T, col.names = T, sep = "\t")
}

#### Performing the analysis for amino acids
#### Making the UMAP for the clones which are shared between individuals
top13 <- read.table("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/Ines_selected_top13_aa.txt", header = FALSE)[, 1]
top13 <- paste0(gsub("NA-", "NA--", top13), "-NA")

dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE)
rm(plot_list)
plot_list <- list()
as.data.frame(CD8_Ag@reductions$wnn.umap@cell.embeddings) -> wnn_umap_embed
for (i in 12:13) {
    CD8_Ag@meta.data$same_clone_diff_ind <- "zz"
    index <- grep((top13[i]), CD8_Ag@meta.data$TRA1_TRA2_TRB1_TRB2_cdraa)
    CD8_Ag@meta.data$same_clone_diff_ind[index] <- CD8_Ag@meta.data$orig.ident[index]
    stopifnot(rownames(CD8_Ag@meta.data) == rownames(wnn_umap_embed))
    wnn_umap_embed$same_clone_diff_ind <- CD8_Ag@meta.data$same_clone_diff_ind
    wnn_umap_embed$same_clone_diff_ind[grep("O7", wnn_umap_embed$same_clone_diff_ind)] <- "zz"
    wnn_umap_embed_ordered <- wnn_umap_embed[order(wnn_umap_embed$same_clone_diff_ind, decreasing = TRUE), ]

    p <- ggplot(wnn_umap_embed_ordered, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = same_clone_diff_ind, size = same_clone_diff_ind)) +
        geom_point() +
        scale_color_manual(values = c("deeppink", "dodgerblue3", "grey87")) +
        scale_size_manual(values = c(2, 2, 0.1)) +
        theme(
            plot.title = element_text(hjust = 0.5, size = 4),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3),
            panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.minor = element_line(colour = "white"),
            panel.grid.major = element_line(colour = "white")
        ) +
        # NoLegend() +
        ggtitle(top13[i])

    plot_list[[i]] <- p
    pdf(paste0(savedir, "UMAP/clone_aa_", i, ".pdf"))
    print(p)
    dev.off()
}

pdf(paste0(savedir, "UMAP/all_sharing_clones.pdf"))
plot_list
dev.off()

### Making a table to identify that each antigen has same epitope
for (i in 1:length(top13)) {
    index <- grep(top13[i], CD8_Ag@meta.data$TRA1_TRA2_TRB1_TRB2_cdraa)
    clone_table <- CD8_Ag@meta.data[index, c(
        "orig.ident", "Ag_range_10", "specificity_score_range_10",
        "max_Ag", "max_specificity_score", "EBV_BALF4",
        "EBV_BMLF1", "EBV_BMRF1", "EBV_BRLF1", "EBV_EBNA1", "EBV_EBNA3C",
        "EBV_LMP2", "VZV_IE62", "VZV_IE63", "TRA1_or_TRA2_TRB1_no_TRB2", "celltypes", "seurat_clusters"
    )]
    write.table(clone_table, paste0(savedir, "Table/clone_table_aa_", i, ".txt"), quote = F, row.names = T, col.names = T, sep = "\t")
}

##### Performing single cell differential between the individuals has same TCR
# CD8_Ag_RNA <- DietSeurat(CD8_Ag, assays = "RNA")
savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/"
clone_files <- list.files("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/clone_sharing_table_aa/", pattern = "clone_table_aa")
for(i in 1:length(clone_files)){
    try({
        clones_common <- read.table(paste0(savedir,"Table/clone_sharing_table_aa/",clone_files[i]), header = TRUE, sep = "\t")
        CD8_Ag_RNA_clones <- subset(CD8_Ag_RNA,cells = rownames(clones_common))
        CD8_Ag_RNA_clones <- NormalizeData(CD8_Ag_RNA_clones)
        samplenames <- unique(CD8_Ag_RNA_clones@meta.data$samplename)
        markers <- FindMarkers(CD8_Ag_RNA_clones, ident.1 = samplenames[1], ident.2 = samplenames[2], group.by = "samplename")
        write.table(markers, paste0(savedir,"Table/clone_sharing_table_aa/scdiff/",samplenames[1],"_vs_",samplenames[2],"_",clone_files[i]),
        sep = "\t", row.names = T, col.names = T, quote = F)
    })
}

savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/"
clone_files <- list.files("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/clone_sharing_table_aa/", pattern = "clone_table_aa")
df$TRA1_TRA2_TRB1_aaa <- "NA"
df$clone_filename <- "NA"

for(i in 1:length(clone_files)){
    try({
        clones <- read.table(paste0(savedir,"Table/clone_sharing_table_aa/",clone_files[i]), header = TRUE, sep = "\t")
        df[i,"TRA1_TRA2_TRB1_aaa"] <- unique(CD8_Ag@meta.data[grep(clones$TRA1_or_TRA2_TRB1_no_TRB2[1],CD8_Ag@meta.data$TRA1_or_TRA2_TRB1_no_TRB2),"TRA1_TRA2_TRB1_TRB2_cdraa"])
        df[i,"clone_filename"] <- clone_files[i]
    })
}

write.table(df, "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/clone_sharing_table_aa/scdiff/clones_id_TCRaa.txt", quote = F, row.names = F, sep ="\t")

savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/"
clone_files <- list.files("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/clone_sharing_table_aa/scdiff/", pattern = "clone_table_aa")
clone_filename <- gsub(".*._Exp0[0-9]_","",clone_files)
rm(cell_list, gene_list)
cell_list <- list()
gene_list <- list()

for(i in 1:length(clone_files)){
    #### extract out the significant genes
    diff_genes <- read.table(paste0("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/clone_sharing_table_aa/scdiff/",clone_files[i]), header = TRUE)
    upreg_genes <- rownames(diff_genes[diff_genes$p_val < 0.01 & diff_genes$avg_log2FC > 0.5,])
    downreg_genes <- rownames(diff_genes[diff_genes$p_val < 0.01 & diff_genes$avg_log2FC < -0.5,])
    sig_genes <- c(upreg_genes, downreg_genes)
    cells <- rownames(read.table(paste0("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/clone_sharing_table_aa/",clone_filename[i]), header = TRUE, sep = "\t"))
    cell_list[[clone_filename[i]]] <- cells
    gene_list[[clone_filename[i]]] <- sig_genes
}

cells_unlisted <- as.data.frame(unlist(cell_list))[,"unlist(cell_list)"]
genes_unlisted <- as.data.frame(unlist(gene_list))[,"unlist(gene_list)"]

cell_list_md <- as.data.frame(unlist(cell_list))
cell_list_md$clones <- gsub(".txt.*.","",rownames(cell_list_md))
colnames(cell_list_md) <- c("cellnames","clones")
# cell_list_md_ordered <- cell_list_md[match(cell_list_md$cellnames, rownames(CD8_Ag_cells@meta.data)),]
cell_list_md_ordered <- cell_list_md[match(rownames(CD8_Ag_cells@meta.data), cell_list_md$cellnames),]
all(cell_list_md_ordered$cellnames == rownames(CD8_Ag_cells@meta.data))

CD8_Ag_cells@meta.data$clones <- cell_list_md_ordered$clones

CD8_Ag_cells <- subset(CD8_Ag, cells =cells_unlisted)
CD8_Ag_cells <- NormalizeData(CD8_Ag_cells)
CD8_Ag_cells <- ScaleData(CD8_Ag_cells)

CD8_Ag_cells@meta.data$clones2 <- gsub("table_aa_","",CD8_Ag_cells@meta.data$clones)
h <- DoHeatmap(CD8_Ag_cells, features = genes_unlisted, group.by = "clones2")

pdf(paste0("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/clone_sharing_table_aa/scdiff/diff_genes.pdf"))
h
dev.off()

CD8_Ag_cells@meta.data$sample_clones <- paste(CD8_Ag_cells@meta.data$clones2,CD8_Ag_cells@meta.data$sample_id)
h <- DoHeatmap(CD8_Ag_cells, features = genes_unlisted, group.by = "sample_clones")

pdf(paste0("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/clone_sharing_table_aa/scdiff/diff_genes_sample_clones.pdf"))
h
dev.off()

write.table(as.data.frame(unlist(gene_list)), "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/clone_sharing_table_aa/scdiff/gene_list.txt", row.names = T, col.names = T, sep = "\t", quote = F)

specific_cells <- rownames(CD8_Ag_cells@meta.data[grep("_1$|_2$|_3$|_12$",CD8_Ag_cells@meta.data$clones2),])
CD8_Ag_cells_specific <- subset(CD8_Ag_cells, cells = specific_cells)
CD8_Ag_cells_specific <- NormalizeData(CD8_Ag_cells_specific)
CD8_Ag_cells_specific <- ScaleData(CD8_Ag_cells_specific)

h <- DoHeatmap(CD8_Ag_cells_specific, features = genes_unlisted, group.by = "sample_clones")

pdf(paste0("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/clone_sharing_table_aa/scdiff/diff_genes_sample_specific_clones.pdf"))
h
dev.off()

#### Performing pseudobulk differential for each clones
savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/"
clone_files <- list.files("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/clone_sharing_table_aa/", pattern = "clone_table_aa")

for(i in 1:length(clone_files)){
    try({
        clones_common <- read.table(paste0(savedir,"Table/clone_sharing_table_aa/",clone_files[i]), header = TRUE, sep = "\t")
        CD8_Ag_RNA_clones <- subset(CD8_Ag_RNA,cells = rownames(clones_common))
        CD8_Ag_RNA_clones <- NormalizeData(CD8_Ag_RNA_clones)
        samplenames <- unique(CD8_Ag_RNA_clones@meta.data$samplename)
        markers <- FindMarkers(CD8_Ag_RNA_clones, ident.1 = samplenames[1], ident.2 = samplenames[2], group.by = "samplename")
        write.table(markers, paste0(savedir,"Table/clone_sharing_table_aa/scdiff/",samplenames[1],"_vs_",samplenames[2],"_",clone_files[i]),
        sep = "\t", row.names = T, col.names = T, quote = F)
    })
}


# Load required library

for(i in 1:length(clone_files)){
    library(edgeR)
    clones_common <- read.table(paste0(savedir,"Table/clone_sharing_table_aa/",clone_files[i]), header = TRUE, sep = "\t")
    CD8_Ag_RNA_clones <- subset(CD8_Ag,cells = rownames(clones_common))
    Treg_cts <- AggregateExpression(CD8_Ag_RNA_clones, group.by = "samplename", assays = 'RNA', slot = "counts", return.seurat = FALSE)
    Treg_cts_2 <- as.data.frame(Treg_cts$RNA)
    colnames(Treg_cts_2) <- gsub("-","_",colnames(Treg_cts_2))
    
    ### Both the samples should have at least some value
    Treg_cts_3 <- Treg_cts_2[Treg_cts_2[,1] > 0 & Treg_cts_2[,2] > 0,]
    # Normalize data using CPM (Counts Per Million)
    dge <- DGEList(counts = Treg_cts_3)
    cpm_data <- cpm(dge, log = TRUE)  # Log2-transformed CPM 
    cpm_data_df <- as.data.frame(cpm_data)
    cpm_data_df$logFC <- log2(cpm_data_df[,1] / cpm_data_df[,2])
    cpm_data_df_order <- cpm_data_df[order(-cpm_data_df$logFC),]
    deg <- cpm_data_df_order[abs(cpm_data_df_order$logFC) > 0.5, ]
    write.table(deg, paste0(savedir,"Table/clone_sharing_table_aa/pseudobulk/",colnames(deg)[1],"_vs_",colnames(deg)[2],"_",clone_files[i]),
                sep = "\t", row.names = T, col.names = T, quote = F)
}

### Now we have identified significant genes we can perform the herarical clustering for all the samples
### Now making a heatmap for all the clones
rm(clones)
clones <- list()
for(i in 1:length(clone_files)){
    clones[[i]] <- rownames(read.table(paste0(savedir,"Table/clone_sharing_table_aa/",clone_files[i]), header = TRUE, sep = "\t"))
}

extract_cells <- unlist(clones)
CD8_Ag_extracted = subset(CD8_Ag, cells = extract_cells)
Treg_cts <- AggregateExpression(CD8_Ag_extracted, group.by = "samplename", assays = 'RNA', slot = "counts", return.seurat = FALSE)
Treg_cts_2 <- as.data.frame(Treg_cts$RNA)
colnames(Treg_cts_2) <- gsub("-","_",colnames(Treg_cts_2))
diff_genes <- read.table("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/clone_sharing_table_aa/pseudobulk/differential_genes.txt", header = FALSE)[,1]
Treg_cts_3 <- Treg_cts_2[match(diff_genes,rownames(Treg_cts_2), nomatch=0),]

dge <- DGEList(counts = Treg_cts_3)
cpm_data <- cpm(dge, log = TRUE)  # Log2-transformed CPM 

metadata <- as.data.frame(colnames(cpm_data))
colnames(metadata) <- "orig_ident"

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

# Convert split elements to dataframe
sample_split <- do.call(rbind, lapply(metadata$orig_ident, split_to_dataframe))
colnames(sample_split) <- c("sampleID","Ag","Age","Age","Gender","Exp")
metadata2 <- cbind(metadata, sample_split)

library(ComplexHeatmap)
library(circlize)
column_anno <- HeatmapAnnotation(
    Age = metadata2$Age,
    Gender = metadata2$Gender,
    col = list(
        Age = c("Y" = "#8fa0ed", "O" = "#e19292"),
        Gender = c("F" = "#09a941", "M" = "#af0000")
    )
)

cpm_data_scaled <- t(scale(t(cpm_data)))
# rld_gene_scaled_ordered <- rld_gene_scaled[match(tot_genes, rownames(rld_gene_scaled)), ]
h <- Heatmap(cpm_data_scaled,
    name = "Z-score", # title of legend
    column_title = "Differential genes", row_title = "Samples",
    row_names_gp = gpar(fontsize = 12), # Text size for row names
    column_names_gp = gpar(fontsize = 8),
    col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    top_annotation = column_anno
)

pdf(paste("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/clone_sharing_table_aa/pseudobulk/heatmap_differential_genes.pdf", sep = ""), height = 7, width = 6)
print(h)
dev.off()



#### Checking epitopes of VDJdb matching with our epitopes
# VDJdb
# EBV-BMLF1	GLCTLVAML
# EBV-BRLF1	YVLDHLIVV
# EBV-LMP2	CLGGLLTMV
# EBV-BALF4 FLDKGTYTL
# EBV-EBNA3C (EBNA3B)	LLDFVRFMGV
# EBV-BMRF1	TLDYKPLSV
# EBV-LMP1	YLQQNWWTL
# EBV-EBNA1	FMVFLQTHI
# VZV-IE62	ALWALPHAA
# VZV-IE63 RLVEDINRV

### Ines is identifying the common TCR between the samples
CD8_Ag_VDJC_sample_Ag <- CD8_Ag@meta.data[c("TRA1_TRA2_TRB1_TRB2_cdrnt", "TRA1_TRA2_TRB1_TRB2_vdjc", "sample_Ag")]
CD8_Ag_VDJC_sample_Ag_noNA <- CD8_Ag_VDJC_sample_Ag[!is.na(CD8_Ag_VDJC_sample_Ag$TRA1_TRA2_TRB1_TRB2_cdrnt), ] ### removing empty clones
CD8_Ag_VDJC_sample_Ag_noNA_2 <- CD8_Ag_VDJC_sample_Ag_noNA[grep("^NA", CD8_Ag_VDJC_sample_Ag_noNA$TRA1_TRA2_TRB1_TRB2_vdjc, invert = TRUE), ] ### Remoing empty TRAV
CD8_Ag_VDJC_sample_Ag_noNA_2_single <- CD8_Ag_VDJC_sample_Ag_noNA_2[grep(",", CD8_Ag_VDJC_sample_Ag_noNA_2$sample_Ag, invert = TRUE), ]

#### Taking only unique clones
unique_clones <- unique(CD8_Ag_VDJC_sample_Ag_noNA_2_single$TRA1_TRA2_TRB1_TRB2_cdrnt)
CD8_Ag_VDJC_sample_Ag_noNA_2_single_unique <- CD8_Ag_VDJC_sample_Ag_noNA_2_single[match(unique_clones, CD8_Ag_VDJC_sample_Ag_noNA_2_single$TRA1_TRA2_TRB1_TRB2_cdrnt), ]

CD8_Ag_VDJC_sample_Ag_noNA_2_single_unique$TRAV <- gsub("__.*.", "", CD8_Ag_VDJC_sample_Ag_noNA_2_single_unique$TRA1_TRA2_TRB1_TRB2_vdjc)
TRAV_samples <- table(CD8_Ag_VDJC_sample_Ag_noNA_2_single_unique$TRAV, CD8_Ag_VDJC_sample_Ag_noNA_2_single_unique$sample_Ag)

column_sum <- colSums(TRAV_samples)

TRAV_samples_norm <- TRAV_samples
for (i in 1:length(column_sum)) {
    TRAV_samples_norm[, i] <- TRAV_samples[, i] / column_sum[i]
}

savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/"
write.table(TRAV_samples_norm, paste0(savedir, "v_gene_usage/TRAV_sampl_norm_unique_clones.txt"), quote = F, sep = "\t")
write.table(TRAV_samples, paste0(savedir, "v_gene_usage/TRAV_samples_unique_clones.txt"), quote = F, sep = "\t")

TRAV_samples_norm <- read.table(paste0(savedir, "v_gene_usage/TRAV_sampl_norm_unique_clones.txt"), header = T)
TRAV_samples_norm_logged= log1p(TRAV_samples_norm)
library(stats)
pca_result <- prcomp(TRAV_samples_norm_logged, scale. = TRUE)

library(ggplot2)
pca_df <- data.frame(pca_result$rotation)
pca_df$Ag <- gsub(".*._EBV_|.*._VZV_", "", rownames(pca_df))

p <- ggplot(pca_df, aes(PC1, PC2, color = Ag, fill = Ag)) +
    geom_point(stroke = 2) + 
    geom_encircle(expand=0, alpha=0.1) +
    xlab(paste0("PC1: ",round(pca_result$sdev[1],2),"% variance")) +
    ylab(paste0("PC2: ",round(pca_result$sdev[2],2),"% variance")) +
    geom_text(aes(label = Ag)) +
    scale_fill_manual(values=c("yellow", "limegreen","steelblue1","red3","gray2","orange","slateblue","azure1","beige","cadetblue3")) +
    scale_color_manual(values=c("yellow3", "green4","steelblue3","red4","gray2","orange3","slateblue4","azure1","beige","cadetblue3")) +
    theme_minimal() +
    ggtitle("TRAV gene usage norm unique clones")

pdf(paste0(paste0(savedir, "v_gene_usage/TRAV_sampl_norm_unique_clones_logged.pdf")), width = 10, height = 10)
p
dev.off()

TRAV_samples <- read.table(paste0(savedir, "v_gene_usage/TRAV_samples_unique_clones.txt"), header = T)
TRAV_samples_norm_logged= log1p(TRAV_samples_norm)
library(stats)
pca_result <- prcomp(TRAV_samples, scale. = TRUE)

library(ggalt)
library(ggplot2)
pca_df <- data.frame(pca_result$rotation)
pca_df$Ag <- gsub(".*._EBV_|.*._VZV_", "", rownames(pca_df))

p <- ggplot(pca_df, aes(PC1, PC2, color = Ag, fill = Ag)) +
    geom_point(stroke = 2) + 
    geom_encircle(expand=0, alpha=0.1) +
    geom_text(aes(label = Ag)) +
    scale_fill_manual(values=c("yellow", "limegreen","steelblue1","red3","lightcoral","orange","slateblue","azure1","beige","cadetblue3")) +
    scale_color_manual(values=c("yellow3", "green4","steelblue3","red4","coral4","orange3","slateblue4","azure1","beige","cadetblue3")) +
    theme_minimal() +
    ggtitle("TRAV gene usage")

pdf(paste0(paste0(savedir, "v_gene_usage/TRAV_sample_clone_unique.pdf")), width = 10, height = 10)
p
dev.off()

p <- ggplot(pca_df, aes(PC1, PC2, color = Ag, fill = Ag)) +
    geom_point(stroke = 2) + 
    geom_encircle(expand=0, alpha=0.1) +
    # geom_text(aes(label = Ag)) +
    scale_fill_manual(values=c("yellow", "limegreen","steelblue1","red3","lightcoral","orange","slateblue","azure1","beige","cadetblue3")) +
    scale_color_manual(values=c("yellow3", "green4","steelblue3","red4","coral4","orange3","slateblue4","azure1","beige","cadetblue3")) +
    theme_minimal() +
    ggtitle("TRAV gene usage")

pdf(paste0(paste0(savedir, "v_gene_usage/TRAV_sample_clone_unique_nolabels.pdf")), width = 10, height = 10)
p
dev.off()

write.table(pca_df, paste0(savedir, "v_gene_usage/TRAV_pca.txt"), sep = "\t", quote = F, row.names = T, col.names = T)

#### Planning to perform by combining all the samples and removing the BALF4 and IE63
TRAV_samples <- read.table(paste0(savedir, "v_gene_usage/TRAV_samples_unique_clones.txt"), header = T)
Ag <- gsub(".*._EBV_|.*._VZV_","",colnames(TRAV_samples)) %>% unique()
req_Ag <- grep("BALF4$|IE63$",Ag, invert = T, value = T)
# TRAV_samples_req <- TRAV_samples[grep("BALF4$|IE63$", colnames(TRAV_samples), invert = T, value = F),]

combined_TRAV <- as.data.frame(matrix(, nrow = nrow(TRAV_samples), ncol = length(req_Ag)))
rownames(combined_TRAV) <- rownames(TRAV_samples)
colnames(combined_TRAV) <- req_Ag

#### Combining all the samples into one for each V gene
for(i in 1:length(req_Ag)){
    combined_TRAV[,req_Ag[i]] <- rowSums(TRAV_samples[,grep(req_Ag[i], colnames(TRAV_samples))]) %>% as.vector()
}

# TRAV_samples_norm_logged= log1p(TRAV_samples_norm)
library(stats)
pca_result <- prcomp(combined_TRAV, scale. = TRUE)

library(ggalt)
library(ggplot2)
pca_df <- data.frame(pca_result$rotation)
pca_df$Ag <- gsub(".*._EBV_|.*._VZV_", "", rownames(pca_df))

p <- ggplot(pca_df, aes(PC1, PC2, color = Ag, fill = Ag)) +
    geom_point(stroke = 2) + 
    # geom_encircle(expand=0, alpha=0.1) +
    geom_text(aes(label = Ag)) +
        xlab(paste0("PC1: ",round(pca_result$sdev[1],2),"% variance")) +
    ylab(paste0("PC2: ",round(pca_result$sdev[2],2),"% variance")) +
    scale_fill_manual(values=c("yellow", "limegreen","steelblue1","red3","lightcoral","orange","slateblue","azure1","beige","cadetblue3")) +
    # scale_color_manual(values=c("yellow3", "green4","steelblue3","red4","coral4","orange3","slateblue4","azure1","beige","cadetblue3")) +
    theme_minimal() +
    ggtitle("TRAV gene usage individual combined")

pdf(paste0(paste0(savedir, "v_gene_usage/TRAV_sample_clone_unique_ind_combined.pdf")))
p
dev.off()

write.table(pca_df, paste0(savedir, "v_gene_usage/TRAV_sample_clone_unique_ind_combined.txt"), sep = "\t", quote = F, row.names = T, col.names = T)

# TRAV_samples_norm_logged= log1p(TRAV_samples_norm)
Ag_sum <- colSums(combined_TRAV)
combined_TRAV_norm <- combined_TRAV

for(i in 1:length(Ag_sum)){
    combined_TRAV_norm[,i] <- combined_TRAV[,i]/Ag_sum[i]
}

colSums(combined_TRAV_norm) #### Sanity Check

combined_TRAV_norm_logged <- log1p(combined_TRAV_norm)

library(stats)
pca_result <- prcomp(combined_TRAV_norm, scale. = TRUE)
# pca_result <- prcomp(combined_TRAV_norm_logged, scale. = TRUE)

library(ggalt)
library(ggplot2)
pca_df <- data.frame(pca_result$rotation)
pca_df$Ag <- gsub(".*._EBV_|.*._VZV_", "", rownames(pca_df))

p <- ggplot(pca_df, aes(PC1, PC2, color = Ag, fill = Ag)) +
    geom_point(stroke = 2) + 
    # geom_encircle(expand=0, alpha=0.1) +
    geom_text(aes(label = Ag)) +
        xlab(paste0("PC1: ",round(pca_result$sdev[1],2),"% variance")) +
    ylab(paste0("PC2: ",round(pca_result$sdev[2],2),"% variance")) +
    scale_fill_manual(values=c("yellow", "limegreen","steelblue1","red3","lightcoral","orange","slateblue","azure1","beige","cadetblue3")) +
    # scale_color_manual(values=c("yellow3", "green4","steelblue3","red4","coral4","orange3","slateblue4","azure1","beige","cadetblue3")) +
    theme_minimal() +
    ggtitle("TRAV normalized gene usage individual combined")

pdf(paste0(paste0(savedir, "v_gene_usage/TRAV_norm_clone_unique_ind_combined_log.pdf")))
p
dev.off()
write.table(pca_df, paste0(savedir, "v_gene_usage/TRAV_sample_clone_unique_ind_combined.txt"), sep = "\t", quote = F, row.names = T, col.names = T)

#### Old and Young
TRAV_samples <- read.table(paste0(savedir, "v_gene_usage/TRAV_samples_unique_clones.txt"), header = T)
Ag <- gsub(".*._EBV_|.*._VZV_","",colnames(TRAV_samples)) %>% unique()
req_Ag <- grep("BALF4$|IE63$",Ag, invert = T, value = T)
# TRAV_samples_req <- TRAV_samples[grep("BALF4$|IE63$", colnames(TRAV_samples), invert = T, value = F),]
### Old and Young
req_Ag_Y = paste0("Y_",req_Ag)
req_Ag_O = paste0("O_",req_Ag)
total_Ag <- c(req_Ag_Y, req_Ag_O)

# TRAV_samples_Y <- TRAV_samples[,grep("^O",colnames(TRAV_samples))]
combined_TRAV <- as.data.frame(matrix(, nrow = nrow(TRAV_samples), ncol = length(total_Ag)))
rownames(combined_TRAV) <- rownames(TRAV_samples)
colnames(combined_TRAV) <- total_Ag

#### Combining all the samples into one for each V gene
for(i in 1:length(req_Ag)){
    TRAV_samples_subset <- TRAV_samples[,grep(req_Ag[i], colnames(TRAV_samples))]
    
    ### Performing for old
    TRAV_samples_subset_O <- TRAV_samples_subset[,grep(paste0("^O"),colnames(TRAV_samples_subset))]
    combined_TRAV[,paste0("O_",req_Ag[i])] = rowSums(TRAV_samples_subset_O)
    
    ### Performing for young
    TRAV_samples_subset_Y <- TRAV_samples_subset[,grep(paste0("^Y"),colnames(TRAV_samples_subset))]
    combined_TRAV[,paste0("Y_",req_Ag[i])] = rowSums(TRAV_samples_subset_Y)
}

# TRAV_samples_norm_logged= log1p(TRAV_samples_norm)
library(stats)
pca_result <- prcomp(combined_TRAV, scale. = TRUE)

library(ggalt)
library(ggplot2)
library(ggrepel)
pca_df <- data.frame(pca_result$rotation)
pca_df$Ag <- gsub("Y_|O_", "", rownames(pca_df))
pca_df$Age <- gsub("_.*.|_.*.", "", rownames(pca_df))

p <- ggplot(pca_df, aes(PC1, PC2, color = Ag, fill = Ag, shape = Age)) +
    geom_point(stroke = 2) + 
    # geom_encircle(expand=0, alpha=0.1) +
    geom_text_repel(aes(label = Ag)) +
        xlab(paste0("PC1: ",round(pca_result$sdev[1],2),"% variance")) +
    ylab(paste0("PC2: ",round(pca_result$sdev[2],2),"% variance")) +
    scale_fill_manual(values=c("yellow", "limegreen","steelblue1","red3","lightcoral","orange","slateblue","azure1","beige","cadetblue3")) +
    # scale_color_manual(values=c("yellow3", "green4","steelblue3","red4","coral4","orange3","slateblue4","azure1","beige","cadetblue3")) +
    theme_minimal() +
    ggtitle("TRAV individual combined")

pdf(paste0(paste0(savedir, "v_gene_usage/TRAV_clone_unique_ind_combined_Age_withlabels.pdf")))
p
dev.off()

write.table(pca_df, paste0(savedir, "v_gene_usage/TRAV_clone_unique_ind_combined_Age.txt"), sep = "\t", quote = F, row.names = T, col.names = T)


# grep("--TRBV", CD8_Ag_VDJC_sample_Ag_noNA$TRA1_TRA2_TRB1_TRB2_vdjc, value = TRUE) %>%
#     gsub(".*.--", "", .) %>%
#     gsub("_.*.", "", .)

CD8_Ag_VDJC_sample_Ag_noNA_2 <- CD8_Ag_VDJC_sample_Ag_noNA[grep("--TRBV", CD8_Ag_VDJC_sample_Ag_noNA$TRA1_TRA2_TRB1_TRB2_vdjc, value = FALSE), ]
CD8_Ag_VDJC_sample_Ag_noNA_2_single <- CD8_Ag_VDJC_sample_Ag_noNA_2[grep(",", CD8_Ag_VDJC_sample_Ag_noNA_2$sample_Ag, invert = TRUE), ]
CD8_Ag_VDJC_sample_Ag_noNA_2_single_unique <- CD8_Ag_VDJC_sample_Ag_noNA_2_single[match(unique_clones, CD8_Ag_VDJC_sample_Ag_noNA_2_single$TRA1_TRA2_TRB1_TRB2_cdrnt), ]

CD8_Ag_VDJC_sample_Ag_noNA_2_single_unique$TRBV <- gsub(".*.--", "", CD8_Ag_VDJC_sample_Ag_noNA_2_single_unique$TRA1_TRA2_TRB1_TRB2_vdjc) %>% gsub("_.*.", "", .)
TRBV_samples <- table(CD8_Ag_VDJC_sample_Ag_noNA_2_single_unique$TRBV, CD8_Ag_VDJC_sample_Ag_noNA_2_single_unique$sample_Ag)

column_sum <- colSums(TRBV_samples)

TRBV_samples_norm <- TRBV_samples
for (i in 1:length(column_sum)) {
    TRBV_samples_norm[, i] <- TRBV_samples[, i] / column_sum[i]
}

write.table(TRBV_samples_norm, paste0(savedir, "v_gene_usage/TRBV_sampl_norm_unique_clones.txt"), quote = F, sep = "\t")
write.table(TRBV_samples, paste0(savedir, "v_gene_usage/TRBV_samples_unique_clones.txt"), quote = F, sep = "\t")

savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/"
TRAV_samples_norm <- read.table(paste0(savedir, "v_gene_usage/TRBV_sampl_norm_unique_clones.txt"), header = T)
TRAV_samples_norm_logged= log1p(TRAV_samples_norm)
library(stats)
pca_result <- prcomp(TRAV_samples_norm_logged, scale. = TRUE)

library(ggplot2)
pca_df <- data.frame(pca_result$rotation)
pca_df$Ag <- gsub(".*._EBV_|.*._VZV_", "", rownames(pca_df))

p <- ggplot(pca_df, aes(PC1, PC2, color = Ag, fill = Ag)) +
    geom_point(stroke = 2) + 
    geom_encircle(expand=0, alpha=0.1) +
    xlab(paste0("PC1: ",round(pca_result$sdev[1],2),"% variance")) +
    ylab(paste0("PC2: ",round(pca_result$sdev[2],2),"% variance")) +
    # geom_text(aes(label = Ag)) +
    scale_fill_manual(values=c("yellow", "limegreen","steelblue1","red3","gray2","orange","slateblue","azure1","beige","cadetblue3")) +
    scale_color_manual(values=c("yellow3", "green4","steelblue3","red4","gray2","orange3","slateblue4","azure1","beige","cadetblue3")) +
    theme_minimal() +
    ggtitle("TRAV gene usage norm unique clones")

pdf(paste0(paste0(savedir, "v_gene_usage/TRBV_sampl_norm_unique_clones_logged_nolabel.pdf")), width = 10, height = 10)
p
dev.off()

write.table(pca_df, paste0(savedir, "v_gene_usage/TRBV_pca.txt"), sep = "\t", quote = F, row.names = T, col.names = T)

#### Planning to perform by combining all the samples and removing the BALF4 and IE63
TRAV_samples <- read.table(paste0(savedir, "v_gene_usage/TRBV_samples_unique_clones.txt"), header = T)
Ag <- gsub(".*._EBV_|.*._VZV_","",colnames(TRAV_samples)) %>% unique()
req_Ag <- grep("BALF4$|IE63$",Ag, invert = T, value = T)
# TRAV_samples_req <- TRAV_samples[grep("BALF4$|IE63$", colnames(TRAV_samples), invert = T, value = F),]

combined_TRAV <- as.data.frame(matrix(, nrow = nrow(TRAV_samples), ncol = length(req_Ag)))
rownames(combined_TRAV) <- rownames(TRAV_samples)
colnames(combined_TRAV) <- req_Ag

#### Combining all the samples into one for each V gene
for(i in 1:length(req_Ag)){
    combined_TRAV[,req_Ag[i]] <- rowSums(TRAV_samples[,grep(req_Ag[i], colnames(TRAV_samples))]) %>% as.vector()
}

# TRAV_samples_norm_logged= log1p(TRAV_samples_norm)
library(stats)
pca_result <- prcomp(combined_TRAV, scale. = TRUE)

library(ggalt)
library(ggplot2)
pca_df <- data.frame(pca_result$rotation)
pca_df$Ag <- gsub(".*._EBV_|.*._VZV_", "", rownames(pca_df))

p <- ggplot(pca_df, aes(PC1, PC2, color = Ag, fill = Ag)) +
    geom_point(stroke = 2) + 
    # geom_encircle(expand=0, alpha=0.1) +
    geom_text(aes(label = Ag)) +
        xlab(paste0("PC1: ",round(pca_result$sdev[1],2),"% variance")) +
    ylab(paste0("PC2: ",round(pca_result$sdev[2],2),"% variance")) +
    scale_fill_manual(values=c("yellow", "limegreen","steelblue1","red3","lightcoral","orange","slateblue","azure1","beige","cadetblue3")) +
    # scale_color_manual(values=c("yellow3", "green4","steelblue3","red4","coral4","orange3","slateblue4","azure1","beige","cadetblue3")) +
    theme_minimal() +
    ggtitle("TRBV gene usage individual combined")

pdf(paste0(paste0(savedir, "v_gene_usage/TRBV_sample_clone_unique_ind_combined.pdf")))
p
dev.off()

write.table(pca_df, paste0(savedir, "v_gene_usage/TRBV_sample_clone_unique_ind_combined.txt"), sep = "\t", quote = F, row.names = T, col.names = T)

#### Planning to perform by combining all the samples and removing the BALF4 and IE63
TRAV_samples <- read.table(paste0(savedir, "v_gene_usage/TRBV_samples_unique_clones.txt"), header = T)
Ag <- gsub(".*._EBV_|.*._VZV_","",colnames(TRAV_samples)) %>% unique()
req_Ag <- grep("BALF4$|IE63$",Ag, invert = T, value = T)
# TRAV_samples_req <- TRAV_samples[grep("BALF4$|IE63$", colnames(TRAV_samples), invert = T, value = F),]
TRAV_samples_Y <- TRAV_samples[,grep("^Y",colnames(TRAV_samples))]

combined_TRAV <- as.data.frame(matrix(, nrow = nrow(TRAV_samples_Y), ncol = length(req_Ag)))
rownames(combined_TRAV) <- rownames(TRAV_samples_Y)
colnames(combined_TRAV) <- req_Ag

#### Combining all the samples into one for each V gene
### perform for Young individuals

for(i in 1:length(req_Ag)){
    combined_TRAV[,req_Ag[i]] <- rowSums(TRAV_samples_Y[,grep((req_Ag[i]), colnames(TRAV_samples_Y))]) %>% as.vector()
}

# TRAV_samples_norm_logged= log1p(TRAV_samples_norm)
library(stats)
pca_result <- prcomp(combined_TRAV, scale. = TRUE)

library(ggalt)
library(ggplot2)
pca_df <- data.frame(pca_result$rotation)
pca_df$Ag <- gsub(".*._EBV_|.*._VZV_", "", rownames(pca_df))

p <- ggplot(pca_df, aes(PC1, PC2, color = Ag, fill = Ag)) +
    geom_point(stroke = 2) + 
    # geom_encircle(expand=0, alpha=0.1) +
    geom_text(aes(label = Ag)) +
        xlab(paste0("PC1: ",round(pca_result$sdev[1],2),"% variance")) +
    ylab(paste0("PC2: ",round(pca_result$sdev[2],2),"% variance")) +
    scale_fill_manual(values=c("yellow", "limegreen","steelblue1","red3","lightcoral","orange","slateblue","azure1","beige","cadetblue3")) +
    # scale_color_manual(values=c("yellow3", "green4","steelblue3","red4","coral4","orange3","slateblue4","azure1","beige","cadetblue3")) +
    theme_minimal() +
    ggtitle("TRBV young gene usage individual combined")

pdf(paste0(paste0(savedir, "v_gene_usage/TRBV_young_sample_clone_unique_ind_combined.pdf")))
p
dev.off()

write.table(pca_df, paste0(savedir, "v_gene_usage/TRBV_young_sample_clone_unique_ind_combined.txt"), sep = "\t", quote = F, row.names = T, col.names = T)


### Old
TRAV_samples_Y <- TRAV_samples[,grep("^O",colnames(TRAV_samples))]

combined_TRAV <- as.data.frame(matrix(, nrow = nrow(TRAV_samples_Y), ncol = length(req_Ag)))
rownames(combined_TRAV) <- rownames(TRAV_samples_Y)
colnames(combined_TRAV) <- req_Ag

#### Combining all the samples into one for each V gene
### perform for Young individuals

for(i in 1:length(req_Ag)){
    combined_TRAV[,req_Ag[i]] <- rowSums(TRAV_samples_Y[,grep((req_Ag[i]), colnames(TRAV_samples_Y))]) %>% as.vector()
}

# TRAV_samples_norm_logged= log1p(TRAV_samples_norm)
library(stats)
pca_result <- prcomp(combined_TRAV, scale. = TRUE)

library(ggalt)
library(ggplot2)
pca_df <- data.frame(pca_result$rotation)
pca_df$Ag <- gsub(".*._EBV_|.*._VZV_", "", rownames(pca_df))

p <- ggplot(pca_df, aes(PC1, PC2, color = Ag, fill = Ag)) +
    geom_point(stroke = 2) + 
    # geom_encircle(expand=0, alpha=0.1) +
    geom_text(aes(label = Ag)) +
        xlab(paste0("PC1: ",round(pca_result$sdev[1],2),"% variance")) +
    ylab(paste0("PC2: ",round(pca_result$sdev[2],2),"% variance")) +
    scale_fill_manual(values=c("yellow", "limegreen","steelblue1","red3","lightcoral","orange","slateblue","azure1","beige","cadetblue3")) +
    # scale_color_manual(values=c("yellow3", "green4","steelblue3","red4","coral4","orange3","slateblue4","azure1","beige","cadetblue3")) +
    theme_minimal() +
    ggtitle("TRBV old gene usage individual combined")

pdf(paste0(paste0(savedir, "v_gene_usage/TRBV_old_sample_clone_unique_ind_combined.pdf")))
p
dev.off()

write.table(pca_df, paste0(savedir, "v_gene_usage/TRBV_old_sample_clone_unique_ind_combined.txt"), sep = "\t", quote = F, row.names = T, col.names = T)

#### Old and Young
TRAV_samples <- read.table(paste0(savedir, "v_gene_usage/TRBV_samples_unique_clones.txt"), header = T)
Ag <- gsub(".*._EBV_|.*._VZV_","",colnames(TRAV_samples)) %>% unique()
req_Ag <- grep("BALF4$|IE63$",Ag, invert = T, value = T)
# TRAV_samples_req <- TRAV_samples[grep("BALF4$|IE63$", colnames(TRAV_samples), invert = T, value = F),]
### Old and Young
req_Ag_Y = paste0("Y_",req_Ag)
req_Ag_O = paste0("O_",req_Ag)
total_Ag <- c(req_Ag_Y, req_Ag_O)

# TRAV_samples_Y <- TRAV_samples[,grep("^O",colnames(TRAV_samples))]
combined_TRAV <- as.data.frame(matrix(, nrow = nrow(TRAV_samples), ncol = length(total_Ag)))
rownames(combined_TRAV) <- rownames(TRAV_samples)
colnames(combined_TRAV) <- total_Ag

#### Combining all the samples into one for each V gene
for(i in 1:length(req_Ag)){
    TRAV_samples_subset <- TRAV_samples[,grep(req_Ag[i], colnames(TRAV_samples))]
    
    ### Performing for old
    TRAV_samples_subset_O <- TRAV_samples_subset[,grep(paste0("^O"),colnames(TRAV_samples_subset))]
    combined_TRAV[,paste0("O_",req_Ag[i])] = rowSums(TRAV_samples_subset_O)
    
    ### Performing for young
    TRAV_samples_subset_Y <- TRAV_samples_subset[,grep(paste0("^Y"),colnames(TRAV_samples_subset))]
    combined_TRAV[,paste0("Y_",req_Ag[i])] = rowSums(TRAV_samples_subset_Y)
}

# TRAV_samples_norm_logged= log1p(TRAV_samples_norm)
library(stats)
pca_result <- prcomp(combined_TRAV, scale. = TRUE)

library(ggalt)
library(ggplot2)
library(ggrepel)
pca_df <- data.frame(pca_result$rotation)
pca_df$Ag <- gsub("Y_|O_", "", rownames(pca_df))
pca_df$Age <- gsub("_.*.|_.*.", "", rownames(pca_df))

p <- ggplot(pca_df, aes(PC1, PC2, color = Ag, fill = Ag, shape = Age)) +
    geom_point(stroke = 2) + 
    # geom_encircle(expand=0, alpha=0.1) +
    geom_text_repel(aes(label = Ag)) +
        xlab(paste0("PC1: ",round(pca_result$sdev[1],2),"% variance")) +
    ylab(paste0("PC2: ",round(pca_result$sdev[2],2),"% variance")) +
    scale_fill_manual(values=c("yellow", "limegreen","steelblue1","red3","lightcoral","orange","slateblue","azure1","beige","cadetblue3")) +
    # scale_color_manual(values=c("yellow3", "green4","steelblue3","red4","coral4","orange3","slateblue4","azure1","beige","cadetblue3")) +
    theme_minimal() +
    ggtitle("TRBV individual combined")

pdf(paste0(paste0(savedir, "v_gene_usage/TRBV_clone_unique_ind_combined_Age_withlabels.pdf")))
p
dev.off()

write.table(pca_df, paste0(savedir, "v_gene_usage/TRBV_clone_unique_ind_combined_Age.txt"), sep = "\t", quote = F, row.names = T, col.names = T)

#### Performing the k-mer motif sequence analysis
# 1. Position Frequency Matrix (PFM)

#     What it is: The PFM shows the frequency of each amino acid (or nucleotide) at each position in the motif. It counts the number of occurrences for each residue at each position across a set of sequences.
#     Best for:
#         Exploratory analysis to see general residue distributions at each position.
#         Visualizing conserved and variable positions in a motif.
#         Initial motif discovery or analyzing unweighted motif occurrence data.
#     Limitations: The PFM doesn't account for sequence length bias or sequence background. It only shows raw frequencies and doesn’t normalize for the likelihood of observing each residue at a position.
#     When to use: If you're primarily interested in raw frequencies and distribution of residues at specific positions, or you're working with a dataset of sequence counts.

# 2. Position Probability Matrix (PPM)

#     What it is: The PPM normalizes the PFM by dividing each residue count at a position by the total number of sequences. This turns the frequencies into probabilities, representing the likelihood of observing each residue at a position.
#     Best for:
#         Comparing sequences with different lengths or different sample sizes, as it normalizes for sequence counts.
#         Visualizing the probability distribution of amino acids at each position.
#     Limitations: While PPM normalizes the data, it doesn’t take into account a background model (i.e., how likely each residue is to appear by chance), which can affect interpretations of sequence conservation.
#     When to use: If you want to compare the relative probabilities of different residues at each position and you don't need a background correction (for example, in motif discovery where you care about relative frequency).

# 3. Position Weight Matrix (PWM)

#     What it is: The PWM is a log-likelihood score that compares the observed probabilities of each amino acid at each position to the expected frequencies based on a background model (usually a random distribution or a biological background model).
#         The log-odds scores give the relative importance of each residue at each position.
#         PWM elements are calculated as:
#         PWM[i,j]=log⁡2(P(observed residue at position j)P(background residue at position j))
#         PWM[i,j]=log2​(P(background residue at position j)P(observed residue at position j)​)
#     Best for:
#         Motif scoring or sequence alignment where you need to quantify how much more likely an amino acid is at a given position compared to random or background distribution.
#         Biological interpretation, particularly when you're interested in understanding conservation and functional residues in a motif.
#         Machine learning applications, like predicting binding affinities or modeling sequence motifs.
#     Limitations: It requires background probabilities (which may not always be available or relevant) and may be affected by small sample sizes or sequence bias.
#     When to use: If you're working on motif scoring, functional analysis, or building predictive models, and you need to account for the likelihood of specific residues at each position compared to a random background.

# 4. Self-Information Matrix

#     What it is: The self-information matrix quantifies how informative each position in the motif is. It calculates the log2 of the inverse probability of each residue's occurrence at each position.
#     Best for:
#         Identifying key positions in a motif that carry the most information about the sequence. Positions with high self-information are typically conserved and play a central role in motif recognition.
#         Understanding which positions are more critical or functional in a sequence, based on their informativeness rather than just frequency or probability.
#     Limitations: It focuses on individual positions and doesn't provide a comprehensive view of the motif as a whole. Self-information does not directly indicate motif affinity or binding strength.
#     When to use: If you are analyzing motif conservation and want to understand how much information each residue provides about the motif, or if you're trying to evaluate which residues are the most significant in the recognition of the motif.

# Which One is Best?
# Use PWM (Position Weight Matrix) if:

#     You are interested in functional analysis, such as determining the importance of each position for recognizing a specific epitope or motif.
#     You need a log-likelihood based comparison that adjusts for background frequencies.
#     You are performing motif scoring or using machine learning to model sequence patterns.

# Use PPM (Position Probability Matrix) if:

#     You are working with large datasets and want to understand the probability distribution of residues at each position.
#     You want to perform comparative analysis across multiple datasets or motifs.
#     You don’t need a background correction (i.e., you are only interested in the relative likelihood of residue occurrence).

# Use PFM (Position Frequency Matrix) if:

#     You want a raw count or a simple summary of how often residues occur at each position.
#     You are at the initial stage of motif discovery and need to see the overall residue frequencies.

# Use Self-Information if:

#     You want to assess how informative or conserved each position is in the motif.
#     You are trying to identify key positions that are essential for the recognition or function of the motif.

### LMP2
LMP2_metadata <- CD8_Ag@meta.data[grep("^EBV_LMP2$", CD8_Ag@meta.data$Ag_range_10), ]
TRA_CDR3 <- gsub("-.*.", "", LMP2_metadata$TRA1_TRA2_TRB1_TRB2_cdraa) %>% gsub(".*._", "", .)
TRA_CDR3_noNA <- grep("^NA$", TRA_CDR3[!is.na(TRA_CDR3)], value = TRUE, invert = TRUE)
TRA_CDR3_noNA_df <- as.data.frame(TRA_CDR3_noNA)
colnames(TRA_CDR3_noNA_df) <- "CDR3.aa"
TRA_CDR3_noNA_df_tib <- as_tibble(TRA_CDR3_noNA_df)

TRA_CDR3_noNA_df_kmers <- getKmers(TRA_CDR3_noNA_df_tib, 4)
TRA_CDR3_noNA_df_kmersprofile <- kmer_profile(TRA_CDR3_noNA_df_kmers, "wei")

p1 <- vis(TRA_CDR3_noNA_df_kmersprofile)
p2 <- vis(TRA_CDR3_noNA_df_kmersprofile, .plot = "seq")

pdf(paste0(paste0(savedir, "kmer/EBV_LMP2_4_PWM.pdf")), width = 12, height = 6)
print(p1 + p2)
dev.off()

#### Performing the analysis with different antigen
Ag <- grep(",", unique(CD8_Ag@meta.data$Ag_range_10), invert = TRUE, value = TRUE)
methods <- c("freq", " ", "wei", "self")
ks <- c(3, 4, 5)

for (i in 1:length(Ag)) {
    LMP2_metadata <- CD8_Ag@meta.data[grep(paste0("^", Ag[i], "$"), CD8_Ag@meta.data$Ag_range_10), ]
    TRA_CDR3 <- gsub("-.*.", "", LMP2_metadata$TRA1_TRA2_TRB1_TRB2_cdraa) %>% gsub(".*._", "", .)
    TRA_CDR3_noNA <- grep("^NA$", TRA_CDR3[!is.na(TRA_CDR3)], value = TRUE, invert = TRUE)
    TRA_CDR3_noNA_df <- as.data.frame(TRA_CDR3_noNA)
    colnames(TRA_CDR3_noNA_df) <- "CDR3.aa"
    TRA_CDR3_noNA_df_tib <- as_tibble(TRA_CDR3_noNA_df)

    dir.create(paste0(savedir, "kmer/TRA/", Ag[i]), showWarnings = FALSE)

    for (k in 1:length(ks)) {
        TRA_CDR3_noNA_df_kmers <- getKmers(TRA_CDR3_noNA_df_tib, ks[k])
        for (j in 1:length(methods)) {
            TRA_CDR3_noNA_df_kmersprofile <- kmer_profile(TRA_CDR3_noNA_df_kmers, methods[j])
            p1 <- vis(TRA_CDR3_noNA_df_kmersprofile)
            p2 <- vis(TRA_CDR3_noNA_df_kmersprofile, .plot = "seq")

            pdf(paste0(paste0(savedir, "kmer/TRA/", Ag[i], "/", Ag[i], "_", ks[k], "_", methods[j], ".pdf")), width = 12, height = 6)
            print(p1 + p2)
            dev.off()
        }
    }
}


#### Running GLIPH2
Ag <- grep(",", unique(CD8_Ag@meta.data$Ag_range_10), invert = TRUE, value = TRUE)
#### Preparing the input
### Extracting out the CDR3 and VDJ for alpha and beta
CD8_Ag@meta.data$cdr3.alpha <- gsub("-.*.", "", CD8_Ag@meta.data$TRA1_TRA2_TRB1_TRB2_cdraa) %>%
    gsub("*.*_", "", .)
CD8_Ag@meta.data$v.alpha <- gsub("-[A-Z].*.", "", CD8_Ag@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc) %>%
    gsub("__.*.", "", .)
CD8_Ag@meta.data$j.alpha <- gsub("-[A-Z].*.", "", CD8_Ag@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc) %>%
    gsub(".*.__", "", .) %>%
    gsub("_TRAC", "", .)

# Beta CDR3 amino acid, V, D, J
CD8_Ag@meta.data$cdr3.beta <- gsub(".*.--", "", CD8_Ag@meta.data$TRA1_TRA2_TRB1_TRB2_cdraa) %>%
    gsub("-.*.", "", .) %>%
    gsub(".*._", "", .)
CD8_Ag@meta.data$v.beta <- gsub(".*.--", "", CD8_Ag@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc) %>%
    gsub("_.*.", "", .)
CD8_Ag@meta.data$d.beta <- gsub(".*.--", "", CD8_Ag@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc) %>%
    gsub("_TRBJ.*.", "", .) %>%
    gsub("*.*_", "", .) %>%
    gsub("NA-NA", "", .)
CD8_Ag@meta.data$j.beta <- gsub(".*.--", "", CD8_Ag@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc) %>%
    gsub("_TRBC[1-2]*.*", "", .) %>%
    gsub(".*._", "", .)

CD8_Ag@meta.data$samplename_Ag <- paste0(CD8_Ag@meta.data$samplename, ":", CD8_Ag@meta.data$Antigen)

gliph_req <- CD8_Ag@meta.data[c("cdr3.beta", "v.beta", "j.beta", "cdr3.alpha", "samplename_Ag")]
gliph_req_noNA <- gliph_req[!is.na(gliph_req$cdr3.beta), ]
gliph_req_noNA2 <- gliph_req_noNA[grep("^NA$", gliph_req_noNA$cdr3.beta, invert = TRUE), ]
savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/gliph2/"
for (i in 1:length(Ag)) {
    gliph_req_noNA2_Ag <- gliph_req_noNA2[grep(paste0(":", Ag[i], "$"), gliph_req_noNA2$samplename_Ag), ]
    gliph_req_noNA2_Ag$j.beta <- gsub("-NA","NA",gliph_req_noNA2_Ag$j.beta)
    gliph_req_noNA2_Ag$v.beta <- gsub("-NA","NA",gliph_req_noNA2_Ag$v.beta)
    # colnames(gliph_req_noNA2) <- c("#CDR3b","TRBV","TRBJ","CDR3a","subject:condition","count")
    write.table(gliph_req_noNA2_Ag, paste0(savedir, Ag[i], "_gliph2.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
}

sort EBV_BMLF1_gliph2.txt | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1}'
ls *gliph2.txt | sed 's/.txt//g' | awk '{print "sort "$1".txt | uniq -c | awk \x27{print $2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$1}\x27 > "$1"_2.txt"}'


#### TRB
Ag <- grep(",", unique(CD8_Ag@meta.data$Ag_range_10), invert = TRUE, value = TRUE)
methods <- c("freq", "prob", "wei", "self")
ks <- c(3, 4, 5)

for (i in 1:length(Ag)) {
    LMP2_metadata <- CD8_Ag@meta.data[grep(paste0("^", Ag[i], "$"), CD8_Ag@meta.data$Ag_range_10), ]
    TRA_CDR3 <- gsub(".*.--", "", LMP2_metadata$TRA1_TRA2_TRB1_TRB2_cdraa) %>%
        gsub(".*._", "", .) %>%
        gsub("-.*.", "", .)
    TRA_CDR3_noNA <- grep("^NA$", TRA_CDR3[!is.na(TRA_CDR3)], value = TRUE, invert = TRUE)
    TRA_CDR3_noNA_df <- as.data.frame(TRA_CDR3_noNA)
    colnames(TRA_CDR3_noNA_df) <- "CDR3.aa"
    TRA_CDR3_noNA_df_tib <- as_tibble(TRA_CDR3_noNA_df)

    dir.create(paste0(savedir, "kmer/TRB/", Ag[i]), showWarnings = FALSE)

    for (k in 1:length(ks)) {
        TRA_CDR3_noNA_df_kmers <- getKmers(TRA_CDR3_noNA_df_tib, ks[k])
        for (j in 1:length(methods)) {
            TRA_CDR3_noNA_df_kmersprofile <- kmer_profile(TRA_CDR3_noNA_df_kmers, methods[j])
            p1 <- vis(TRA_CDR3_noNA_df_kmersprofile)
            p2 <- vis(TRA_CDR3_noNA_df_kmersprofile, .plot = "seq")

            pdf(paste0(paste0(savedir, "kmer/TRB/", Ag[i], "/", Ag[i], "_", ks[k], "_", methods[j], ".pdf")), width = 12, height = 6)
            print(p1 + p2)
            dev.off()
        }
    }
}

#### Generating the k-means or PCA plot for the 3 amino acids
Ag <- grep(",", unique(CD8_Ag@meta.data$Ag_range_10), invert = TRUE, value = TRUE)
CD8_Ag@meta.data$TRB_CDR3 <- gsub(".*.--", "", CD8_Ag@meta.data$TRA1_TRA2_TRB1_TRB2_cdraa) %>%
    gsub(".*._", "", .) %>%
    gsub("-.*.", "", .)
CD8_Ag@meta.data$TRA_CDR3 <- gsub("-.*.", "", CD8_Ag@meta.data$TRA1_TRA2_TRB1_TRB2_cdraa) %>% gsub(".*._", "", .)

required_table <- CD8_Ag@meta.data[, c("Ag_range_10", "TRA_CDR3", "TRB_CDR3")]


for (i in 1:length(Ag)) {
    required_table_Ag <- required_table[grep(paste0("^", Ag[i], "$"), required_table$Ag_range_10), ]

    ## performing for TRA
    Ag_TRA_CDR3 <- grep("^NA$", required_table_Ag[!is.na(required_table_Ag[, "TRA_CDR3"]), "TRA_CDR3"], value = TRUE, invert = TRUE)
    Ag_TRA_CDR3_df <- as_tibble(Ag_TRA_CDR3)
    colnames(Ag_TRA_CDR3_df) <- "CDR3.aa"
    assign(paste0(Ag[i], "_TRA"), Ag_TRA_CDR3_df)

    ## performing for TRB
    Ag_TRB_CDR3 <- grep("^NA$", required_table_Ag[!is.na(required_table_Ag[, "TRB_CDR3"]), "TRB_CDR3"], value = TRUE, invert = TRUE)
    Ag_TRB_CDR3_df <- as_tibble(Ag_TRB_CDR3)
    colnames(Ag_TRB_CDR3_df) <- "CDR3.aa"
    assign(paste0(Ag[i], "_TRB"), Ag_TRB_CDR3_df)
}

#### TRA
## Since IE63 and BALF4 is very low
VZV_IE63_TRA_low <- VZV_IE63_TRA
EBV_BALF4_TRA_low <- EBV_BALF4_TRA
rm(VZV_IE63_TRA)
rm(EBV_BALF4_TRA)


TRA_samples <- ls(pattern = "_TRA$")
rm(TRA_list)
TRA_list <- list()
for (i in 1:length(TRA_samples)) {
    TRA_list[[TRA_samples[i]]] <- get(TRA_samples[i])
}

#### Performing the k-mers on the TRA
TRA_kmers <- getKmers(TRA_list, 4)
TRA_kmers[is.na(TRA_kmers)] <- 0

TRA_kmers_df <- as.data.frame(TRA_kmers)
rownames(TRA_kmers_df) <- TRA_kmers_df$Kmer
TRA_kmers_df <- TRA_kmers_df[, -1]

### Since there is a lot of dropout convert frequencies to log
TRA_kmers_df_log <- log1p(TRA_kmers_df)

library(stats)
pca_result <- prcomp(TRA_kmers_df_log, scale. = TRUE)

library(ggplot2)
pca_df <- data.frame(pca_result$rotation)
pca_df$Ag <- gsub("_TRA", "", rownames(pca_df))

p <- ggplot(pca_df, aes(PC1, PC2, color = Ag)) +
    geom_point() +
    geom_text(aes(label = Ag)) +
    theme_minimal() +
    ggtitle("PCA of CDR3 TRA k-mer Data")

pdf(paste0(paste0(savedir, "v_gene_usage/kmer_4_PCA_noBALF4_IE63.pdf")), width = 8, height = 6)
p
dev.off()

write.table(pca_df, paste0(savedir, "kmer/TRA/4_mers/kmer_4_PCA_noBALF4_IE63.txt"), sep = "\t", quote = F, row.names = T, col.names = T)

### Heatmap clustering
TRA_kmers_df_log_scale <- t(scale(t(TRA_kmers_df_log)))

library(ComplexHeatmap)

p <- Heatmap(
    as.matrix(TRA_kmers_df_log_scale),
    # name = "signature_contri",
    col = colorRampPalette(c("blue", "white", "red"))(50),
    column_names_gp = gpar(fontsize = 10),
    row_names_gp = gpar(fontsize = 12),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    heatmap_legend_param = list(title_position = "topcenter", legend_direction = "horizontal"),
    # top_annotation = col_anno # Use the column annotation as top annotation
)

pdf("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/kmer/TRA/4_mers/TRA_kmer_clustering_2.pdf")
draw(p, heatmap_legend_side = "top")
dev.off()

write.table(TRA_kmers_df_log_scale,
    "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/kmer/TRA/4_mers/TRA_4mer_clustering_2.txt",
    quote = F,
    row.names = T,
    col.names = T,
    sep = "\t"
)

### Since there are lot of drop out that is why we are missing VZV-IE63, LMP2, and BALF4 are separated out
## TRB
## Since IE63 and BALF4 is very low
VZV_IE63_TRB_low <- VZV_IE63_TRB
EBV_BALF4_TRB_low <- EBV_BALF4_TRB
rm(VZV_IE63_TRB)
rm(EBV_BALF4_TRB)

TRB_samples <- ls(pattern = "_TRB$")
rm(TRB_list)
TRB_list <- list()
for (i in 1:length(TRB_samples)) {
    TRB_list[[TRB_samples[i]]] <- get(TRB_samples[i])
}

#### Performing the k-mers on the TRB
TRB_kmers <- getKmers(TRB_list, 3)
TRB_kmers[is.na(TRB_kmers)] <- 0

TRB_kmers_df <- as.data.frame(TRB_kmers)
rownames(TRB_kmers_df) <- TRB_kmers_df$Kmer
TRB_kmers_df <- TRB_kmers_df[, -1]

### Since there is a lot of dropout convert frequencies to log
TRB_kmers_df_log <- log1p(TRB_kmers_df)

library(stats)
pca_result <- prcomp(TRB_kmers_df_log, scale. = TRUE)

library(ggplot2)
pca_df <- data.frame(pca_result$rotation)
pca_df$Ag <- gsub("_TRB", "", rownames(pca_df))

p <- ggplot(pca_df, aes(PC1, PC2, label = Ag, color = Ag)) +
    geom_point() +
    geom_text(aes(label = Ag)) +
    theme_minimal() +
    ggtitle("PCA of CDR3 TRB k-mer Data")

pdf(paste0(paste0(savedir, "kmer/TRB/3_mers/kmer_3_PCA_noBALF4_IE63.pdf")), width = 8, height = 6)
p
dev.off()
write.table(pca_df, paste0(savedir, "kmer/TRB/3_mers/kmer_3_PCA_noBALF4_IE63.txt"), sep = "\t", quote = F, row.names = T, col.names = T)

### Heatmap clustering
TRB_kmers_df_log_scale <- t(scale(t(TRB_kmers_df_log)))

library(ComplexHeatmap)

p <- Heatmap(
    as.matrix(TRB_kmers_df_log_scale),
    # name = "signature_contri",
    col = colorRampPalette(c("blue", "white", "red"))(50),
    column_names_gp = gpar(fontsize = 10),
    row_names_gp = gpar(fontsize = 12),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    heatmap_legend_param = list(title_position = "topcenter", legend_direction = "horizontal"),
    # top_annotation = col_anno # Use the column annotation as top annotation
)

pdf("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/kmer/TRB/3_mers/TRB_3mer_clustering_2.pdf")
draw(p, heatmap_legend_side = "top")
dev.off()

write.table(TRB_kmers_df_log_scale,
    "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/kmer/TRB/3_mers/TRB_3mer_clustering_2.txt",
    quote = F,
    row.names = T,
    col.names = T,
    sep = "\t"
)


### After running GLIPH2, performing the downstream analysis
# Load the required packages
library(msa)
library(ggseqlogo)
library(Biostrings)  # Needed to handle sequences

CDR3_BMLF1 = read.table("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/gliph2/Gliph2_output/downstream/CDR3_BMLF1_uniq.txt", header = FALSE)[,1]
seqs <- AAStringSet(CDR3_BMLF1)
alignment <- msa(seqs, method = "ClustalW")
aligned_seqs <- as.character(alignment)

savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/gliph2/Gliph2_output/downstream/"
pdf(paste0(savedir, "BMLF1_sequencelogo_probability.pdf"), width = 10, height = 6)
ggseqlogo(aligned_seqs, seq_type='aa', method = "prob")
dev.off()

savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/gliph2/Gliph2_output/downstream/"
pdf(paste0(savedir, "BMLF1_sequencelogo_bits.pdf"), width = 10, height = 6)
ggseqlogo(aligned_seqs, seq_type='aa', method = "bits")
dev.off()


#### Including all the TCR beta no unique
CDR3_BMLF1 = read.table("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/gliph2/EBV_BMLF1_gliph2.txt", header = FALSE)[,1]
seqs <- AAStringSet(CDR3_BMLF1)
alignment <- msa(seqs, method = "ClustalW")
aligned_seqs <- as.character(alignment)

savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/gliph2/"
pdf(paste0(savedir, "BMLF1_sequencelogo_allTCRs_probability.pdf"), width = 10, height = 6)
ggseqlogo(aligned_seqs, seq_type='aa', method = "prob")
dev.off()

pdf(paste0(savedir, "BMLF1_sequencelogo_allTCRs_bits.pdf"), width = 10, height = 6)
ggseqlogo(aligned_seqs, seq_type='aa', method = "bits")
dev.off()


#### Including all the TCR beta unique unique
CDR3_BMLF1 = read.table("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/gliph2/EBV_BMLF1_unique.txt", header = FALSE)[,1]
seqs <- AAStringSet(CDR3_BMLF1)
alignment <- msa(seqs, method = "Muscle")
aligned_seqs <- as.character(alignment)

savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/gliph2/"
pdf(paste0(savedir, "BMLF1_sequencelogo_unique_TCRs_probability.pdf"), width = 10, height = 6)
ggseqlogo(aligned_seqs, seq_type='aa', method = "prob")
dev.off()

pdf(paste0(savedir, "BMLF1_sequencelogo_unique_TCRs_bits.pdf"), width = 10, height = 6)
ggseqlogo(aligned_seqs, seq_type='aa', method = "bits")
dev.off()


#### Running for all the samples
# ls *interchange.csv | sed 's/.csv//g' | awk '{print "awk -F \",\" \x27{if ($3 < 0.05) print $14}\x27 "$1".csv | sed \x27/^$/d\x27 | sort | uniq > "$1"_CDR3_aa_sig.txt"}' > extracting_sig_TCR.sh
CDR3_sig <- list.files("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/gliph2/Gliph2_output/", pattern = "_interchange_CDR3_aa_sig.txt", full.names = TRUE)
savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/gliph2/Gliph2_output/downstream/"
Ag <- gsub("_GLIPH2.*.","",basename(CDR3_sig))

for(i in 1:length(CDR3_sig)){
    seqs <- AAStringSet(read.table(CDR3_sig[i], header = FALSE)[,1])
    alignment <- msa(seqs, method = "ClustalW")
    aligned_seqs <- as.character(alignment)
    
    pdf(paste0(savedir,Ag[i],"_sequencelogo_gliph2_CDR3_sig_TCRs_probability.pdf"), width = 10, height = 6)
    print(ggseqlogo(aligned_seqs, seq_type='aa', method = "prob"))
    dev.off()
    
    pdf(paste0(savedir, Ag[i],"_sequencelogo_gliph2_CDR3_sig_TCRs_bits.pdf"), width = 10, height = 6)
    print(ggseqlogo(aligned_seqs, seq_type='aa', method = "bits"))
    dev.off()
}


#### Running for all the samples
# ls *interchange.csv | sed 's/.csv//g' | awk '{print "awk -F \",\" \x27{if ($3 < 0.05) print $14}\x27 "$1".csv | sed \x27/^$/d\x27 | sort | uniq > "$1"_CDR3_aa_sig.txt"}' > extracting_sig_TCR.sh
CDR3_sig <- list.files("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/gliph2/Gliph2_output/old/", pattern = "_old_CDR3_aa_sig.txt", full.names = TRUE)
savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/gliph2/Gliph2_output/old/downstream/"
dir.create(savedir, showWarnings = FALSE)
Ag <- gsub("_old.*.","",basename(CDR3_sig))

for(i in 1:length(CDR3_sig)){
    try({
    seqs <- AAStringSet(read.table(CDR3_sig[i], header = FALSE)[,1])
    alignment <- msa(seqs, method = "ClustalW")
    aligned_seqs <- as.character(alignment)
    
    pdf(paste0(savedir,Ag[i],"_old_sequencelogo_gliph2_CDR3_sig_TCRs_probability.pdf"), width = 10, height = 6)
    print(ggseqlogo(aligned_seqs, seq_type='aa', method = "prob"))
    dev.off()
    
    pdf(paste0(savedir, Ag[i],"_old_sequencelogo_gliph2_CDR3_sig_TCRs_bits.pdf"), width = 10, height = 6)
    print(ggseqlogo(aligned_seqs, seq_type='aa', method = "bits"))
    dev.off()
    })
}

#### Running for young the samples
# ls *.csv | sed 's/.csv//g' | awk '{print "awk -F \",\" \x27{if ($3 < 0.05) print $14}\x27 "$1".csv | sed \x27/^$/d\x27 | sort | uniq > "$1"_CDR3_aa_sig.txt"}' > extracting_sig_TCR.sh
CDR3_sig <- list.files("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/gliph2/Gliph2_output/young/", pattern = "_young_CDR3_aa_sig.txt", full.names = TRUE)
savedir <- "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/gliph2/Gliph2_output/young/downstream/"
dir.create(savedir, showWarnings = FALSE)
Ag <- gsub("_young.*.","",basename(CDR3_sig))

for(i in 1:length(CDR3_sig)){
    try({
    seqs <- AAStringSet(read.table(CDR3_sig[i], header = FALSE)[,1])
    alignment <- msa(seqs, method = "ClustalW")
    aligned_seqs <- as.character(alignment)
    
    pdf(paste0(savedir,Ag[i],"_young_sequencelogo_gliph2_CDR3_sig_TCRs_probability.pdf"), width = 10, height = 6)
    print(ggseqlogo(aligned_seqs, seq_type='aa', method = "prob"))
    dev.off()
    
    pdf(paste0(savedir, Ag[i],"_young_sequencelogo_gliph2_CDR3_sig_TCRs_bits.pdf"), width = 10, height = 6)
    print(ggseqlogo(aligned_seqs, seq_type='aa', method = "bits"))
    dev.off()
    })
}

#### For affinity Proxy estimation
RNA_norm <- t(CD8_Ag@assays$RNA@data[c("CD3E","CD4"),])
colnames(RNA_norm) <- c("CD3E_RNA","CD4_RNA")
RNA_imputed <- t(CD8_Ag@assays$MAGIC_RNA@data[c("CD3E","CD4"),])
colnames(RNA_imputed) <- c("CD3E_imputed","CD4_imputed")

all(rownames(RNA_norm) == rownames(RNA_imputed))
RNA_imputed_df <- as.data.frame(RNA_imputed)
RNA_norm_df <- as.data.frame(RNA_norm)

RNA_norm_df$CD3E_imputed <- RNA_imputed_df$CD3E_imputed
RNA_norm_df_2 <- RNA_norm_df[,-2]
write.table()

##### Performing UCell for the Innate genes that Ines has provided
library(ArchR)
library(UCell)
library(ggplot2)
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
DefaultAssay(CD8_Ag) <- "MAGIC_RNA"
rm(markers)
markers <- list()
markers[["innate"]] <- read.table("/diazlab/data3/.abhinav/resources/gene_list/innateness.txt", header = FALSE, sep= "\t")[,1] %>% unique() %>% toupper()
CD8_Ag <- AddModuleScore(CD8_Ag, features = markers, slot="data")
colnames(CD8_Ag@meta.data)[grep("Cluster",colnames(CD8_Ag@meta.data))] <- "innate"
p <- FeaturePlot(CD8_Ag, features = "innate", reduction = "wnn.umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)

pdf("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/featureplot/innate_genes.pdf", width =5.5, height = 5.5)
p
dev.off()

innate_Ucell <- CD8_Ag@meta.data[,c("innate","CDR3")]
write.table(innate_Ucell, "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/innate_Ucell.txt", col.names = T, row.names = T, sep = "\t", quote = F)

DefaultAssay(CD8_Ag) <- "RNA"
rm(markers)
markers <- list()
markers[["innate"]] <- read.table("/diazlab/data3/.abhinav/resources/gene_list/innateness.txt", header = FALSE, sep= "\t")[,1] %>% unique() %>% toupper()
CD8_Ag <- AddModuleScore(CD8_Ag, features = markers, slot="data")
colnames(CD8_Ag@meta.data)[grep("Cluster",colnames(CD8_Ag@meta.data))] <- "innate_RNA"
p <- FeaturePlot(CD8_Ag, features = "innate_RNA", reduction = "wnn.umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)

pdf("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/featureplot/innate_genes_RNA.pdf", width =5.5, height = 5.5)
p
dev.off()

innate_Ucell <- CD8_Ag@meta.data[,c("innate","innate_RNA")]
colnames(innate_Ucell) <- c("innate_imputed","innate_unimputed")
write.table(innate_Ucell, "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/Table/innate_Ucell_imputed_unimputed.txt", col.names = T, row.names = T, sep = "\t", quote = F)

# TRSS redo: could you run this analysis again to make sure we have the correct data? please use all clones (defined as nucleotide CDR1+2+3 TRA1, A2, B1) using only 
# memory cells and only expanded cells, combining all clones of an Ag in Y vs O (not individual specific) , the make TRSS across cell types (TEMRA, CM, EM1, …)
library(Seurat)
library(dplyr)

CD8_metadata <- CD8_Ag@meta.data[!is.na(CD8_Ag@meta.data$TRA1_or_TRA2_TRB1_no_TRB2),] #### removing empty cells
CD8_metadata2 <- CD8_metadata[grep("naïve",CD8_metadata$celltypes,invert=TRUE),] #### removing the naive cells
CD8_metadata3 <- CD8_metadata2[grep(",",CD8_metadata2$Ag_range_10,invert=TRUE),] ### removing TCR cross presentation

CD8_metadata3$Ag <- gsub(".*._","",CD8_metadata3$Ag_range_10)
Ag = unique(CD8_metadata3$Ag)
CD8_metadata3$Age_Ag <- paste0(CD8_metadata3$Age,"_", CD8_metadata3$Ag)
AgeAg <- unique(CD8_metadata3$Age_Ag)

req_cells <- rownames(CD8_metadata3)
CD8_Ag_req = subset(CD8_Ag, cells = req_cells)

stopifnot(all(rownames(CD8_metadata3) == rownames(CD8_Ag_req@meta.data)))
CD8_Ag_req@meta.data <- CD8_metadata3

#### Cell types
savedir = "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/revision/TCR/"
dir.create(savedir, showWarnings = FALSE)
source("/diazlab/data3/.abhinav/resources/all_scripts/R/VZV_clonal_sharing_gt_1_celltypes.R")
CD8_Ag_req$Age_Ag <- paste(CD8_Ag_req@meta.data$Age ,CD8_Ag_req@meta.data$Ag, sep="_")
Age_Ag = grep(",",unique(CD8_Ag_req$Age_Ag),value=TRUE,invert=TRUE)
clus_length = length(unique(CD8_Ag_req@meta.data$celltypes))
clus_factor = c("CM","EM1","EM2","EM3","TEMRA","MAIT","HELIOS_high")

for (i in 1:length(Age_Ag)) {
  try({TCR_AJ_VZV_gt_1_CT(object = CD8_Ag_req, savedir = savedir, 
                          clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2", 
                          group_col = "Age_Ag", group_val = Age_Ag[i], 
                          split_col = "orig.ident", column_name = c("id","Ag","Age","Age_number","gender","Run"),
                          total_clusters = clus_length, clus_fact = clus_factor)
  })
}
