#### We are developing our model to identify TCR and epitope prediction 
### Extracting out the data from seurat object
# CD8_Ag_metadata <- CD8_Ag@meta.data[,c("TRA1_TRA2_TRB1_TRB2_cdraa","TRA1_or_TRA2_TRB1_no_TRB2","TRA1_TRA2_TRB1_TRB2_fwr_aa","EBV_BALF4","EBV_BMLF1","EBV_BMRF1","EBV_BRLF1","EBV_EBNA1","EBV_EBNA3C","EBV_LMP2","VZV_IE62","VZV_IE63")]
savedir = "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/"
write.table(CD8_Ag_metadata, paste0(savedir,"rawdata/TCR_EBV_peptide.txt"), sep = "\t", quote = F, col.names = T, row.names = F)