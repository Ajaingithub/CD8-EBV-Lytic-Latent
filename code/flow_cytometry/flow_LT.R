### Only Young
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Lutian_analysis/"
dir.create(savedir, showWarnings = FALSE)
Ag_mem_Y <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Lutian_analysis/young/young_freq", header = TRUE, sep = "\t", row.names = 1)

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

write.table(pc_with_metadata, file = paste(savedir,"young/PCA_coordinate_matched_samples_with_metadata_Young.txt",sep = ""),
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
plot_list[[1]] <-ggplot(pc_with_metadata, aes(PC1, PC2, label=samples, color=Ag, shape=Ag)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(pc_summary_imp[2,1])*100,"% variance")) +
  ylab(paste0("PC2: ",(pc_summary_imp[2,2])*100,"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

plot_list[[2]] <-ggplot(pc_with_metadata, aes(PC3, PC4, label=samples, color=Ag, shape=Ag)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(pc_summary_imp[2,3])*100,"% variance")) +
  ylab(paste0("PC2: ",(pc_summary_imp[2,4])*100,"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

plot_list[[3]] <-ggplot(pc_with_metadata, aes(PC5, PC6, label=samples, color=Ag, shape=Ag)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(pc_summary_imp[2,5])*100,"% variance")) +
  ylab(paste0("PC2: ",(pc_summary_imp[2,6])*100,"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))


pdf(paste(savedir,"young/PCA_Y_norm_samples.pdf",sep = ""),width = 9, height = 8)
plot_list
dev.off()

### Antigen Lutian Analysis
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Lutian_analysis/"
dir.create(savedir, showWarnings = FALSE)
Ag_mem_Y <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Lutian_analysis/young/young_freq_more_info", header = TRUE, sep = "\t", row.names = 1)

### combining the samplename for the metadata
rownames(Ag_mem_Y) <- paste(rownames(Ag_mem_Y),Ag_mem_Y$Age, Ag_mem_Y$Gender, sep = "_")
### Only the Samplename
Ag_mem_Y <- Ag_mem_Y[,4:ncol(Ag_mem_Y)]

metadata <- data.frame(samples = rownames(Ag_mem_Y))
rownames(metadata) <- metadata$samples

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

# Convert split elements to dataframe
metadata_2 <- do.call(rbind, lapply(metadata$samples, split_to_dataframe))
colnames(metadata_2) <- c("samplename","Ag","Age","Sex")
metadata_2$samples <- metadata$samples


Ag <- unique(metadata_2$Ag)

df <- data.frame(matrix(nrow = length(Ag), ncol = length(Ag)))
rownames(df) <- Ag
colnames(df) <- Ag

for (i in 1:(length(Ag))) {
  for (j in 1:length(Ag)) {
    print(paste("Performing Lutian analysis for ",Ag[i]," and ",Ag[j]))
    data <- Ag_mem_Y[grep(paste(Ag[i],"|",Ag[j],sep = ""),rownames(Ag_mem_Y)),]
    
    if (Ag[i] == Ag[j]) {
      print(paste(Ag[i]," is same as ",Ag[j]))
    } else if (Ag[i] != Ag[j]){
      metadata_3 <- do.call(rbind, lapply(rownames(data), split_to_dataframe))
      colnames(metadata_3) <- c("samplename","Ag","Age","Sex")
      Sex <- metadata_3$Sex
      Ag_var <- metadata_3$Ag
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

write.table(df, paste(savedir,"young/Lutian_Y_only_Ag.txt",sep = ""), sep = "\t", quote = F, row.names = T, col.names = T)


#### Young and Old #####
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Lutian_analysis/"
dir.create(savedir, showWarnings = FALSE)
Ag_mem_Y <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Lutian_analysis/young_old/YO_freq", header = TRUE, sep = "\t", row.names = 1)

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

write.table(pc_with_metadata, file = paste(savedir,"young_old/PCA_coordinate_matched_samples_with_metadata_Young.txt",sep = ""),
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
plot_list[[1]] <-ggplot(pc_with_metadata, aes(PC1, PC2, label=samples, color=Ag, shape=Ag)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(pc_summary_imp[2,1])*100,"% variance")) +
  ylab(paste0("PC2: ",(pc_summary_imp[2,2])*100,"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

plot_list[[2]] <-ggplot(pc_with_metadata, aes(PC3, PC4, label=samples, color=Ag, shape=Ag)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(pc_summary_imp[2,3])*100,"% variance")) +
  ylab(paste0("PC2: ",(pc_summary_imp[2,4])*100,"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

plot_list[[3]] <-ggplot(pc_with_metadata, aes(PC5, PC6, label=samples, color=Ag, shape=Ag)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(pc_summary_imp[2,5])*100,"% variance")) +
  ylab(paste0("PC2: ",(pc_summary_imp[2,6])*100,"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))


pdf(paste(savedir,"young_old/PCA_YO_norm_samples.pdf",sep = ""),width = 9, height = 8)
plot_list
dev.off()

### Antigen Lutian Analysis
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Lutian_analysis/"
dir.create(savedir, showWarnings = FALSE)
Ag_mem_Y <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/flow_data/Lutian_analysis/young_old/YO_freq_more_info", header = TRUE, sep = "\t", row.names = 1)

### combining the samplename for the metadata
rownames(Ag_mem_Y) <- paste(rownames(Ag_mem_Y),Ag_mem_Y$Age,Ag_mem_Y$Age_num, Ag_mem_Y$Gender, sep = "_")
### Only the Samplename
Ag_mem_Y <- Ag_mem_Y[,5:ncol(Ag_mem_Y)]

metadata <- data.frame(samples = rownames(Ag_mem_Y))
rownames(metadata) <- metadata$samples

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

# Convert split elements to dataframe
metadata_2 <- do.call(rbind, lapply(metadata$samples, split_to_dataframe))
colnames(metadata_2) <- c("samplename","Ag","Age","Age_num","Sex")
metadata_2$samples <- metadata$samples
Ag <- unique(metadata_2$Ag)

df <- data.frame(matrix(nrow = length(Ag), ncol = 1))
rownames(df) <- Ag
colnames(df) <- "pvalue"

Ag_unique <- Ag

for (i in 1:length(Ag_unique)) {
  print(paste("performing Lutian analysis for ",Ag_unique[i]))
  Ag_mem_subset <- Ag_mem_Y[grep(Ag_unique[i],rownames(Ag_mem_Y)),]
  metadata_3 <- do.call(rbind, lapply(rownames(Ag_mem_subset), split_to_dataframe))
  colnames(metadata_3) <- c("samplename","Ag","Age","Age_num","Sex")
  
  Ags <- metadata_3$Ag
  Age <- metadata_3$Age_num
  id <-  metadata_3$samplename
  Sex <- metadata_3$Sex
  Age <- as.numeric(Age)
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

write.table(df, paste(savedir,"young_old/Lutian_pvalue.txt",sep = ""),quote = F, row.names = T, col.names = T, sep = "\t")
