#  https://www.10xgenomics.com/support/software/cell-ranger/tutorials/cr-tutorial-multi-beam-t 
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/refdata-gex-GRCh38-2020-A/
# [gene-expression]
# ref,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/refdata-gex-GRCh38-2020-A/
# 
# [feature]
# ref,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/training/5k_BEAM-T_Human_A0201_B0702_PBMC_5pv2_Multiplex_count_feature_reference.csv
# 
# [vdj]
# ref,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0
# 
# [libraries]
# fastq_id,fastqs,lanes,feature_types
# beamt_human_A0201_B0702_pbmc_ag,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/training/5k_BEAM-T_Human_A0201_B0702_PBMC_5pv2_fastqs/antigen_capture,1|2,Antigen Capture
# beamt_human_A0201_B0702_pbmc_vdj,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/training/5k_BEAM-T_Human_A0201_B0702_PBMC_5pv2_fastqs/vdj,1|2,VDJ-T
# beamt_human_A0201_B0702_pbmc_gex,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/training/5k_BEAM-T_Human_A0201_B0702_PBMC_5pv2_fastqs/gex,1|2,Gene Expression
# 
# [antigen-specificity]
# control_id,mhc_allele
# negative_control_A0201,HLA-A*02:01
# negative_control_B0702,HLA-B*07:02
#
## Made a bash file and saved in beam_run.sh
# source ~/.bashrc
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/bin/cellranger multi --id=beam_t_run --csv=../5k_BEAM-T_Human_A0201_B0702_PBMC_5pv2_Multiplex_config.csv
#
# Provide the permission to run chmod 777 beam_run.sh
# 
# Running it in slurm
# sbatch --mail-user=jain.abhinav@mayo.edu --mail-type=END,FAIL,TIME_LIMIT_50 --job-name=BEAM_t_run --partition=cpu-short --nodes=1 --tasks=5 --time=3-00:00:00 --mem=50G --chdir=./ --output=logs/%x.%N.%j.stdout --error=logs/%x.%N.%j.stderr <<EOF
# #!/bin/bash
# ./beam_run.sh
# EOF

### After running lets load the data
library("Matrix")
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/CR7_scCITESeq_QC.R")

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run06/Downstream/"
dir_run06 = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run06/post_processing/"
dir_run06_files <- list.files(dir_run06, pattern = "BB|IMAG")
dir_run06_files_2 <- dir_run06_files[grep(".csv",dir_run06_files, ignore.case = T, invert = T)]
dir_run06_2 <- vector()
for (i in 1:length(dir_run06_files_2)) {
  dir_run06_2[i] <- paste(dir_run06,dir_run06_files_2[i],"/outs/per_sample_outs/",dir_run06_files_2[i],"/",sep = "")
}

id =  basename(dir_run06_files_2)
samplepath = dir_run06_2[which(basename(dir_run06_files_2) %in% id)]
Sample = basename(samplepath)
table <- data.frame(Filter = integer(), Doublets = integer()) # Making an empty dataframe to put the count before and after filtering

i=1
for (i in 1:length(samplepath)) {
  table[i,1] <- paste(scCITE_QC(samplepath[i],Sample[i], saveDir = savedir, min_cells=3, min_genes=200, max_genes=5000, 
                                mitopercent=10, ADT_UMI=100000),collapse = ",")
}

rownames(table) <- Sample
dir.create(paste(savedir,"counts/",sep = ""), showWarnings = FALSE)
write.table(table, paste(savedir,"counts/cellfiltered.txt",sep = ""), sep = "\t", quote = F, row.names = T, col.names = T)

### Identify the resolution of the object for each sample with the cluster will be used for the doublet finder ####
## Before running this please change the source file, try to save the RDS files
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/NN_clustertree.R")
object <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "_RNA.RDS",full.names = TRUE)
object <- object[grep("clustertree",object,invert = TRUE)]
dims =20
samplename <- gsub("_RNA.RDS","",basename(object))
Assay = "RNA"
process = "nn_clustertree"

for (i in 1:length(object)) {
  obj = readRDS(object[i])
  obj = nearest_neigbour(Obj =  obj, Dims = dims, saveDir = savedir, samplename = samplename[i], Assay = Assay, process = process)
  assign(paste(samplename[i],"_nn",sep = ""), obj)
}

### Running and Filtering the Doublet Finder from the RNA and ADT ####
# https://github.com/satijalab/seurat/issues/1565 resolution
# Keeping minimum resolution to 0.4 and max 1.2 based on the Seurat essence, there is no correct clustering parameter,
# either you will over or under cluster your data. To compensate for what makes biological sense in the context of your experiment, you can merge certain clusters togethe
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITEseq_Doublet_Finder.R")
dims <- 20
res <- 0.4
# obj <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "_RNA.RDS",full.names = TRUE)
# obj <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "_nn_clustertree_RNA.RDS",full.names = TRUE)
obj = ls(pattern="_nn")
obj2 <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "ADT.RDS",full.names = TRUE)

samplename <- gsub("_nn","",basename(obj))
table2 <- read.table(paste(savedir,"counts/cellfiltered.txt",sep = ""),sep = "\t",  header = TRUE)

i=1
for (i in 1:length(obj)) {
  table2[i,2]<- paste(doublet_scCITEseq(Obj = get(obj[i]), dims = dims,res = res, saveDir = savedir,Obj2 = obj2[i],
                                        samplename = samplename[i], process = "doublet_finder", Assay = "RNA"), collapse =",")
}

dir.create(paste(savedir,"counts/",sep = ""), showWarnings = FALSE)
write.table(table2, paste(savedir,"counts/Doublets_Removed.txt",sep = ""), sep = "\t", quote = F, row.names = T, col.names = T)





