#### Extracting out the TCR and Ag extraction for uploading to the VDJ database. So creating TCR alpha amino acid, TCR beta amino acid
library(Seurat)
library(dplyr)
# CD8_Ag <- readRDS("/diazlab/data3/.abhinav/.immune/cancer_combined/project/resources/GSE275633_CD8_Antigen_BEAM_T.RDS")
CD8_required <- CD8_Ag@meta.data[, c("TRA1_TRA2_TRB1_TRB2_cdraa", "TRA1_TRA2_TRB1_TRB2_vdjc", "Ag_range_10")]
CD8_required_noNA <- na.omit(CD8_required)

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
