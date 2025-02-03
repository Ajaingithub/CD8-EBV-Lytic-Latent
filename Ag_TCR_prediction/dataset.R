### Buiding our model for EBV using TITAN architecture
library(Seurat)
library(dplyr)
# CD8_Ag <- readRDS("/diazlab/data3/.abhinav/.immune/cancer_combined/project/resources/GSE275633_CD8_Antigen_BEAM_T.RDS")
CD8_required <- CD8_Ag@meta.data[, c("TRA1_TRA2_TRB1_TRB2_cdraa", "TRA1_TRA2_TRB1_TRB2_vdjc", "Ag_range_10")]
CD8_required_noNA <- na.omit(CD8_required)

### Removing TRB2
CD8_required_noNA_TRB2 <- CD8_required_noNA[grep("-NA$", CD8_required_noNA$TRA1_TRA2_TRB1_TRB2_cdraa), ]

#### Let's start with only CDR3 of TRB2
CD8_required_noNA_TRB2$CDR3_TRB2 <- gsub("-NA", "", CD8_required_noNA_TRB2$TRA1_TRA2_TRB1_TRB2_cdraa) %>% gsub(".*._", "", .)
CD8_required_noNA_TRB2_Ag_spec <- CD8_required_noNA_TRB2[grep(",", CD8_required_noNA_TRB2$Ag_range_10, invert = TRUE), ]

### Removing any antigen that has less than 1000 TCRs for model to learn better
Ag_TCR_number <- table(CD8_required_noNA_TRB2_Ag_spec$Ag_range_10) %>% as.data.frame()
Ag_req <- Ag_TCR_number[Ag_TCR_number$Freq > 1000, "Var1"] %>% as.vector()

CD8_required_noNA_TRB2_Ag_spec[grep(paste0("^", Ag_req, "$", collapse = "|"), CD8_required_noNA_TRB2_Ag_spec$Ag_range_10), ]
CD8_required_noNA_TRB2_Ag_spec_req <- CD8_required_noNA_TRB2_Ag_spec[grep(paste0("^", Ag_req, "$", collapse = "|"), CD8_required_noNA_TRB2_Ag_spec$Ag_range_10), ]

patterns <- c("EBV_BMLF1", "EBV_BRLF1", "EBV_BMRF1", "EBV_BALF4", "EBV_EBNA1", "EBV_EBNA3C", "EBV_LMP1", "EBV_LMP2", "VZV_IE62")
replacements <- c("GLCTLVAML", "YVLDHLIVV", "TLDYKPLSV", "FLDKGTYTL", "FMVFLQTHI", "LLDFVRFMGV", "YLQQNWWTL", "CLGGLLTMV", "ALWALPHAA")

library(stringr)
CD8_required_noNA_TRB2_Ag_spec_req$eptiope <- CD8_required_noNA_TRB2_Ag_spec_req$Ag_range_10
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    CD8_required_noNA_TRB2_Ag_spec_req$eptiope <- str_replace_all(CD8_required_noNA_TRB2_Ag_spec_req$eptiope, pattern, replacements[i])
}

TCR_epitope <- paste(CD8_required_noNA_TRB2_Ag_spec_req$CDR3_TRB2, CD8_required_noNA_TRB2_Ag_spec_req$eptiope, sep = "_") %>% unique()
write.table(as.data.frame(TCR_epitope), "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/TCR_epitope.txt",
    sep = "\t", row.names = F, col.names = F, quote = F
)

TCR_epitope <- read.table("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/TCR_epitope.txt", header = FALSE)
TCR_epitope_df <- as.data.frame(TCR_epitope)
colnames(TCR_epitope_df) <- "TCR_epitope"
TCR_epitope_df$TCRs <- gsub("_.*.", "", TCR_epitope_df$TCR_epitope)
TCR_epitope_df$eptiope <- gsub(".*._", "", TCR_epitope_df$TCR_epitope)

patterns <- c("ALWALPHAA", "FMVFLQTHI", "GLCTLVAML", "LLDFVRFMGV", "YVLDHLIVV")
replacements <- as.character(c(1, 2, 3, 4, 5))

library(stringr)
TCR_epitope_df$eptiope_num <- TCR_epitope_df$eptiope
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    TCR_epitope_df$eptiope_num <- str_replace_all(TCR_epitope_df$eptiope_num, pattern, replacements[i])
}

patterns <- c("ALWALPHAA", "FMVFLQTHI", "GLCTLVAML", "LLDFVRFMGV", "YVLDHLIVV")
replacements <- as.character(c(1, 1, 1, 1, 1))

library(stringr)
TCR_epitope_df$label <- TCR_epitope_df$eptiope
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    TCR_epitope_df$label <- str_replace_all(TCR_epitope_df$label, pattern, replacements[i])
}

TCR_epitope_df$TCR_num <- c(1:nrow(TCR_epitope_df))
TCR_epitope_df$index <- c(0:(nrow(TCR_epitope_df) - 1))

#### Shifflfing the data
req_TCR <- TCR_epitope_df[, c("TCRs", "eptiope")]
TITAN_input <- TCR_epitope_df[, c("eptiope_num", "TCR_num", "label")]

write.table(TITAN_input, "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/full_titan_dataset.txt", quote = F, row.names = F, col.names = F, sep = "\t")

TITAN_input <- read.table("/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/full_titan_shuffled.txt", header = FALSE, sep = "\t")
set.seed(123)
test_indices <- sample(nrow(TITAN_input), size = 0.2 * nrow(TITAN_input))

# Split the dataframe
test_df <- TITAN_input[test_indices, ] # 20% test set
train_df <- TITAN_input[-test_indices, ] # Remaining 80% training set

test_df$index <- c(0:(nrow(test_df) - 1))
test_df2 <- test_df[, c(4, 1, 2, 3)]
train_df$index <- c(0:(nrow(train_df) - 1))
train_df2 <- train_df[, c(4, 1, 2, 3)]

write.table(test_df2, "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/Ag_test.csv", sep = ",", row.names = F, col.names = F, quote = F)
write.table(train_df2, "/diazlab/data3/.abhinav/.immune/CD8-EBV-Lytic-Latent/Ag_TCR_prediction/data/TITAN/Ag_train.csv", sep = ",", row.names = F, col.names = F, quote = F)
