## Testing to run on Slurm mode
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/bin/cellranger count --id test --transcriptome /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/external/cellranger_tiny_ref --fastqs /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/external/cellranger_tiny_fastq --sample tinygex --jobmode=slurm --maxjobs=200 
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/bin/cellranger multi --id=Exp01_First_4 --csv=./Required_files/BEAM_T_multiconfig_Exp01_First_4.csv --jobmode=slurm --maxjobs=200
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/bin/cellranger multi --id=Exp02_First_4 --csv=./Required_files/BEAM_T_multiconfig_Exp02_First_4.csv --jobmode=slurm --maxjobs=200
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/bin/cellranger multi --id=Exp03_First_4 --csv=./Required_files/BEAM_T_multiconfig_Exp03_First_4.csv --jobmode=slurm --maxjobs=200
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/bin/cellranger multi --id=Exp03_second_4 --csv=./Required_files/BEAM_T_multiconfig_Exp03_second_4.csv --jobmode=slurm --maxjobs=200
## Preprocessing at this location
#
# Following and automating this to demultiplex 5 prime library 
# https://www.10xgenomics.com/cn/resources/analysis-guides/demultiplexing-and-analyzing-5%E2%80%99-immune-profiling-libraries-pooled-with-hashtags
# 
# source ~/.bashrc
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/bin/cellranger multi --id=Run04_CR7 --csv=/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run04/Analysis/preprocessing/Run04_CMO_reference.csv
# 
# [gene-expression]
# reference,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/refdata-gex-GRCh38-2020-A/
# cmo-set,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run04/Analysis/preprocessing/Run04_CMO_reference.csv
# 
# [libraries]
# fastq_id,fastqs,feature_types
# CoV2-04-GEX,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run04/rawdata,Gene Expression
# CoV2-04-CITE,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run04/rawdata,Multiplexing Capture
# 
# [samples]
# sample_id,cmo_ids
# DMSO_O_BB22037_Run04,HTO1
# DMSO_Y_IMAG015_Run04,HTO2
# DMSO_Y_IMAG031_Run04,HTO3
# DMSO_Y_IMAG017_Run04,HTO4
# Sprotein_O_BB22037_Run04,HTO5
# Sprotein_Y_IMAG015_Run04,HTO6
# Sprotein_Y_IMAG031_Run04,HTO7
# Sprotein_Y_IMAG017_Run04,HTO8
# 
# id,name,read,pattern,sequence,feature_type
# HTO1,HTO1,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,GTCAACTCTTTAGCG,Multiplexing Capture
# HTO2,HTO2,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,TGATGGCCTATTGGG,Multiplexing Capture
# HTO3,HTO3,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,TTCCGCCTCTCTTTG,Multiplexing Capture
# HTO4,HTO4,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,AGTAAGTTCAGCGTA,Multiplexing Capture
# HTO5,HTO5,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,AAGTATCGTTTCGCA,Multiplexing Capture
# HTO6,HTO6,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,GGTTGCCAGATGTCA,Multiplexing Capture
# HTO7,HTO7,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,TGTCTTTCCTGCCAG,Multiplexing Capture
# HTO8,HTO8,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,CTCCTCTGCAATTAC,Multiplexing Capture
#
## When we have Slurm to submit the jobs
## For fastqc
# #!/bin/bash
# 
# #SBATCH --mail-user=jain.abhinav@mayo.edu
# #SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
# #SBATCH --job-name=Run06
# #SBATCH --partition=cpu-short
# #SBATCH --nodes=6
# #SBATCH --tasks=5
# #SBATCH --time=3-00:00:00
# #SBATCH --mem=50G
# #SBATCH --chdir ./
# #SBATCH --output logs/%x.%N.%j.stdout
# #SBATCH --error logs/%x.%N.%j.stderr
# 
# #export OMP_NUM_THREADS=32
# echo "slurm job id $SLURM_JOB_ID"
# echo "slurm submit dir $SLURM_SUBMIT_ DIR"
# echo "Number of CPUs $SLURM_CPUS_ON_NODE"
# echo "slurm job nodelist $SLURM_JOB_NODELIST"
# 
# ./Run06_CR7.sh
# Now we have to make the bam2fastq files.
# Running this command
# #!/bin/bash
# input the file path
# input_path=$1
# output_dir=$2
# 
# # Adding the subdirectory
# mkdir -p ${input_path}"shell_scripts"
# input_path2=${input_path}"outs/per_sample_outs/"
# 
# # extracting the number of reads from any one of the summary metrics file
# dir_name=$(ls $input_path2 | head -1)
# metric_file=${input_path2}${dir_name}"/metrics_summary.csv"
# echo $metric_file
# cd ${input_path}"shell_scripts"
# lib_reads=$(grep "Number of reads in the library," $metric_file | awk -F "library," '{print $NF}' | sed 's/"//g' | sed 's/,//g')
# echo $lib_reads
# 
# # Creating the output directory
# mkdir -p $output_dir"bam_to_fastq"
# 
# # Now you can use the input_path variable in your commands
# # Make a shell script for generating the bam2fastq files
# ls $input_path2 | awk -v input_path="$input_path2" -v output_dir="$output_dir" -v lib_reads=$lib_reads '{print "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/cellranger-7.0.0/lib/bin/bamtofastq --reads-per-fastq=" lib_reads" "input_path $1"/count/sample_alignments.bam " output_dir"bam_to_fastq/"$1}' > bam2fastq.sh
# 
# # Moving to the directory
# # cd ${output_dir}"bam_to_fastq"
# 
# # Generating bash script for each sample
# awk '{print "head -"NR" bam2fastq.sh | tail -1 > sample_"NR"_b2f.sh"}' bam2fastq.sh > generating_bash.sh
# sh generating_bash.sh
# 
# # adding header
# echo "source ~/.bashrc" > header
# ls *_b2f.sh | sed 's/_/\t/g' | awk '{print "cat header "$1"_"$2"_b2f.sh > "$1"_"$2"_b2f_2.sh"}'  > adding_header.sh
# sh adding_header.sh
# 
# # Final job running command
# # ls *_2.sh | sed 's/.sh//g' | awk '{print "qsub -l h_vmem=32G -pe threaded 20 -N "$1" -q 1-day -o "$1".log -e "$1".err  -m ae -M jain.abhinav@mayo.edu -cwd "$1".sh"}' > final_running.sh
# # sh final_running.sh
# 
# ### bash executing_commands.sh /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run05/preprocessing/Run05_CR7_nosec/ /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run05/00_fastq/
# #!/bin/bash
# 
# #SBATCH --mail-user=jain.abhinav@mayo.edu
# #SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
# #SBATCH --job-name=Exp01_b2f_slurm_array
# #SBATCH --partition=cpu-short
# #SBATCH --time=02:00:00
# #SBATCH --mem=100G
# #SBATCH --chdir ./
# #SBATCH --output logs/%x.%A.%j.%a.stdout
# #SBATCH --error  logs/%x.%A.%j.%a.stderr
# 
# #SBATCH --array=1-4                     	# Array range and number of simultanous jobs
# file=$(ls *b2f_2.sh | sed -n ${SLURM_ARRAY_TASK_ID}p)
# 
# sh $file
# write this file as sample_cells
# Y1_31_M_Exp01	7576
# Y2_27_M_Exp01	6038
# O1_70_F_Exp01	7732
# O2_71_M_Exp01	5627

# Y3_28_F_Exp02	2607
# Y4_34_F_Exp02	9032
# Y5_37_M_Exp02	4079
# O3_78_F_Exp02	4923
# O4_68_F_Exp02	4929

# Y6_27_F_Exp03	4280
# Y7_39_M_Exp03	4616
# O5_65_M_Exp03	3092
# O6_73_F_Exp03	5360
# O7_73_F_Exp03	8649
# 
# this is Y6_EBV_27_F_Exp03_multiconfig.csv
# [gene-expression]
# reference,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/refdata-gex-GRCh38-2020-A/
# force-cells,4280
# check-library-compatibility,FALSE
# 
# [feature]
# reference,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/Exp03_feature_reference.csv
# 
# [vdj]
# reference,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0
# 
# [libraries]
# fastq_id,fastqs,lanes,feature_types
# bamtofastq,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/raw_data/fastq/igc1.salk.edu/NSC_Ines_231218_10X/bam_to_fastq/Y6_EBV_27_F_Exp03/Exp03_First_4_0_1_HKJNKDSX7,any,Gene Expression
# Exp03A_CITE,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/raw_data/fastq/igc1.salk.edu/NSC_Ines_231218_10X,any,Antibody Capture
# Exp03B_CITE,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/raw_data/fastq/igc1.salk.edu/NSC_Ines_231218_10X,any,Antibody Capture
# Exp03A_TCR,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/raw_data/fastq/igc1.salk.edu/NSC_Ines_231218_10X,any,VDJ-T
# Exp03B_TCR,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/raw_data/fastq/igc1.salk.edu/NSC_Ines_231218_10X,any,VDJ-T
# Exp03A_BEAM,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/raw_data/fastq/igc1.salk.edu/NSC_Ines_231218_10X,any,Antigen Capture
# Exp03B_BEAM,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/raw_data/fastq/igc1.salk.edu/NSC_Ines_231218_10X,any,Antigen Capture
# 
# [antigen-specificity]
# control_id,mhc_allele
# Neg,HLA-A*02:01
# # Change the number of from 166 to whatever you have put
# #!/bin/bash
# # input the file path
# input_path=$1
# # Generating multiple config files
# awk -v input_path=$input_path '{print "sed \x27s/4280/"$2"/g\x27 " input_path" | sed \x27s/Y6_EBV_27_F_Exp03/"$1"/g\x27 > "$1"_final_multiplex_config.csv"}' sample_cells > replace.sh
# sh replace.sh
# 
# ## adding a header and generating the command to submit jobs
# echo "source ~/.bashrc" > header
# ls *_final_multiplex_config.csv | sed 's/_final_/\t/g' | awk '{print "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/bin/cellranger multi --id="$1" --csv="$1"_final_multiplex_config.csv --jobmode=slurm --maxjobs=200"}' > command.sh
# awk '{print "head -"NR" command.sh | tail -1 > Sample_"NR".sh"}' command.sh > making_bash.sh
# sh making_bash.sh
# ls Sample_*.sh | sed 's/.sh//g' | awk '{print "cat header "$1".sh > "$1"_2.sh"}' > adding_header.sh
# sh adding_header.sh
## ls *_[0-9]_2.sh | sed 's/.sh//g' | awk '{print "qsub -l h_vmem=32G -pe threaded 20 -N "$1"  -o "$1".log -e "$1".err -m ae -M jain.abhinav@mayo.edu -cwd "$1".sh"}' > final_running.sh
## sh final_running.sh
# Run it like this
# bash executing_command_2.sh ./DMSO_O_BB22037_Run04_config.csv
# 
# source ~/.bashrc
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/bin/cellranger multi --id=O1_EBV_70_F_Exp01 --csv=O1_EBV_70_F_Exp01_final_multiplex_config.csv --jobmode=slurm --maxjobs=200
# Since Exp03 failed. We are currently exploring to increase the memory and thread of the cellranger using the json file
# {
#   "SC_MULTI_CS.SC_MULTI_CORE.MULTI_REPORTER.VDJ_T_REPORTER.BEAM_ANALYZER.CALCULATE_ANTIGEN_SPECIFICITY": {
#     "chunk.mem_gb": 10,
#     "chunk.threads": 5
#   }
# }
# 
# We ran the cellranger count for each experiment like Exp01A, Exp01B, Exp02A, Exp02B, Exp03A, Exp03B
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/bin/cellranger count --id=Exp01A --transcriptome=/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/refdata-gex-GRCh38-2020-A/ --sample=Exp01A_GEX --fastqs=/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/raw_data/fastq/igc1.salk.edu/NSC_Ines_231218_10X/ --jobmode=slurm --maxjobs=200
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/bin/cellranger count --id=Exp02A --transcriptome=/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/refdata-gex-GRCh38-2020-A/ --sample=Exp02A_GEX --fastqs=/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/raw_data/fastq/igc1.salk.edu/NSC_Ines_231218_10X/ --jobmode=slurm --maxjobs=200
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/bin/cellranger count --id=Exp02B --transcriptome=/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/refdata-gex-GRCh38-2020-A/ --sample=Exp02B_GEX --fastqs=/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/raw_data/fastq/igc1.salk.edu/NSC_Ines_231218_10X/ --jobmode=slurm --maxjobs=200
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/bin/cellranger count --id=Exp03A --transcriptome=/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/refdata-gex-GRCh38-2020-A/ --sample=Exp03A_GEX --fastqs=/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/raw_data/fastq/igc1.salk.edu/NSC_Ines_231218_10X/ --jobmode=slurm --maxjobs=200
# /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/cellranger-7.1.0/bin/cellranger count --id=Exp03B --transcriptome=/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/refdata-gex-GRCh38-2020-A/ --sample=Exp03B_GEX --fastqs=/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/raw_data/fastq/igc1.salk.edu/NSC_Ines_231218_10X/ --jobmode=slurm --maxjobs=200
# sbatch --mail-user=jain.abhinav@mayo.edu --mail-type=BEGIN,END,FAIL,TIME_LIMIT_50 --job-name=E1A_VCF --partition=cpu-short --nodes=1 --tasks=1 --time=4-00:00:00 --cpus-per-task=8 --mem=100G --chdir=./ --output=logs/%x.%N.%j.stdout --error=logs/%x.%N.%j.stderr <<EOF
# #!/bin/bash
# ./vcf_generation.sh
# EOF
# sbatch --mail-user=jain.abhinav@mayo.edu --mail-type=BEGIN,END,FAIL,TIME_LIMIT_50 --job-name=E1B_VCF --partition=cpu-short --nodes=1 --tasks=1 --time=4-00:00:00 --cpus-per-task=8 --mem=100G --chdir=./ --output=logs/%x.%N.%j.stdout --error=logs/%x.%N.%j.stderr <<EOF
# #!/bin/bash
# ./vcf_generation.sh
# EOF

### Got good recovery

### Downstream ####
## First identifying the number of cells present for each individuals and based on that in each antigen
## We have three experiments starting with Exp01
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01_First_4/post_processing/"
exp <- "Exp01"
files_list <- list.files(maindir)
samplename <- files_list[grep(paste("_",exp,"$",sep = ""),files_list)]
req_files <- paste(maindir,samplename,"/outs/per_sample_outs/",samplename,"/antigen_analysis/antigen_specificity_scores.csv",sep = "")
Ag <- unique(read.csv(req_files[1], header = TRUE)[,2])
rm(df)
df <- as.data.frame(matrix(ncol= length(Ag)))
colnames(df) <- Ag
# rownames(df) <- Ag

rm(df_list)
df_list <- list()
for (i in 1:length(req_files)) {
  Ag_file <- read.csv(req_files[i], header = TRUE)
  Ag_high_score <- Ag_file[Ag_file$antigen_specificity_score > 50,]
  Ag_high_score_df <- as.data.frame(table(Ag_high_score$antigen))
  rownames(Ag_high_score_df) <- Ag_high_score_df$Var1
  df_list[[samplename[i]]] <- Ag_high_score_df
}

df_list_bind <- bind_cols(df_list)
df_list_bind_2 <- df_list_bind[,-grep("Var1",colnames(df_list_bind))]
colnames(df_list_bind_2) <- names(df_list)
Exp01_df <- df_list_bind_2


maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02_First_4/post_processing/"
exp <- "Exp02"
files_list <- list.files(maindir)
samplename <- files_list[grep(paste("_",exp,"$",sep = ""),files_list)]
req_files <- paste(maindir,samplename,"/outs/per_sample_outs/",samplename,"/antigen_analysis/antigen_specificity_scores.csv",sep = "")
Ag <- unique(read.csv(req_files[1], header = TRUE)[,2])
rm(df)
df <- as.data.frame(matrix(ncol= length(Ag)))
colnames(df) <- Ag
# rownames(df) <- Ag

rm(df_list)
df_list <- list()
for (i in 1:length(req_files)) {
  Ag_file <- read.csv(req_files[i], header = TRUE)
  Ag_high_score <- Ag_file[Ag_file$antigen_specificity_score > 50,]
  Ag_high_score_df <- as.data.frame(table(Ag_high_score$antigen))
  rownames(Ag_high_score_df) <- Ag_high_score_df$Var1
  df_list[[samplename[i]]] <- Ag_high_score_df
}

df_list_bind <- bind_cols(df_list)
df_list_bind_2 <- df_list_bind[,-grep("Var1",colnames(df_list_bind))]
colnames(df_list_bind_2) <- names(df_list)

Exp02_df <- df_list_bind_2

maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/"
exp <- "Exp03"
files_list <- list.files(maindir)
samplename <- files_list[grep(paste("_",exp,"$",sep = ""),files_list)]
req_files <- paste(maindir,samplename,"/outs/per_sample_outs/",samplename,"/antigen_analysis/antigen_specificity_scores.csv",sep = "")
Ag <- unique(read.csv(req_files[1], header = TRUE)[,2])
rm(df)
df <- as.data.frame(matrix(ncol= length(Ag)))
colnames(df) <- Ag
# rownames(df) <- Ag

rm(df_list)
df_list <- list()
for (i in 1:length(req_files)) {
  Ag_file <- read.csv(req_files[i], header = TRUE)
  Ag_high_score <- Ag_file[Ag_file$antigen_specificity_score > 50,]
  Ag_high_score_df <- as.data.frame(table(Ag_high_score$antigen))
  rownames(Ag_high_score_df) <- Ag_high_score_df$Var1
  df_list[[samplename[i]]] <- Ag_high_score_df
}

df_list_bind <- bind_cols(df_list)
df_list_bind_2 <- df_list_bind[,-grep("Var1",colnames(df_list_bind))]
colnames(df_list_bind_2) <- names(df_list)

Exp03_df <- df_list_bind_2

Exp01_02_df <- merge(Exp01_df,Exp02_df,by="row.names")
rownames(Exp01_02_df) <- Exp01_02_df$Row.names
Exp01_02_df_2 <- Exp01_02_df[,-1]
Exp01_02_03_df <- merge(Exp01_02_df_2,Exp03_df, by = "row.names")

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/"
dir.create(savedir, showWarnings = FALSE)
dir.create(paste(savedir,"Table/",sep = ""), showWarnings = FALSE)
write.table(Exp01_02_03_df, paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/Table/Exp01_02_03_df.txt"),
            quote = F, row.names = F, col.names = T, sep = "\t")

Exp01_02_df_melted <- melt(Exp01_02_03_df)
colnames(Exp01_02_df_melted) <- c("Ag","ID","Ag_pos_cells")

dir.create(savedir, showWarnings = FALSE)
Exp01_02_df_melted$Ag <- factor(Exp01_02_df_melted$Ag)
p <- ggplot(Exp01_02_df_melted, aes(x=Ag, y=Ag_pos_cells, color = Ag)) + geom_boxplot() +
  geom_point() + ggtitle("Exp01 02 and 03 Ag Specificity score > 50") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"Table/Exp01_02_03_df.pdf",sep = ""), width = 8.5, height = 7)
p
dev.off()

#### Finding the cells belongs to which clusters
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01_First_4/post_processing/"
exp <- "Exp01"
files_list <- list.files(maindir)
samplename <- files_list[grep(paste("_",exp,"$",sep = ""),files_list)]
req_files <- paste(maindir,samplename,"/outs/per_sample_outs/",samplename,"/antigen_analysis/antigen_specificity_scores.csv",sep = "")
Ag <- unique(read.csv(req_files[1], header = TRUE)[,2])
# rownames(df) <- Ag

rm(df_list)
df_list <- list()
for (i in 1:length(req_files)) {
  Ag_file <- read.csv(req_files[i], header = TRUE)
  Ag_high_score <- Ag_file[Ag_file$antigen_specificity_score > 50,]
  Ag_high_score$barcode_antigen <- paste(Ag_high_score$barcode, Ag_high_score$antigen, sep = "_")
  df_list[[samplename[i]]] <- Ag_high_score$barcode_antigen
}

df_Exp01 <- as.data.frame(unlist(df_list))

maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02_First_4/post_processing/"
exp <- "Exp02"
files_list <- list.files(maindir)
samplename <- files_list[grep(paste("_",exp,"$",sep = ""),files_list)]
req_files <- paste(maindir,samplename,"/outs/per_sample_outs/",samplename,"/antigen_analysis/antigen_specificity_scores.csv",sep = "")
Ag <- unique(read.csv(req_files[1], header = TRUE)[,2])
# rownames(df) <- Ag

rm(df_list)
df_list <- list()
for (i in 1:length(req_files)) {
  Ag_file <- read.csv(req_files[i], header = TRUE)
  Ag_high_score <- Ag_file[Ag_file$antigen_specificity_score > 50,]
  Ag_high_score$barcode_antigen <- paste(Ag_high_score$barcode, Ag_high_score$antigen, sep = "_")
  df_list[[samplename[i]]] <- Ag_high_score$barcode_antigen
}

df_Exp02 <- as.data.frame(unlist(df_list))

maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/"
exp <- "Exp03"
files_list <- list.files(maindir)
samplename <- files_list[grep(paste("_",exp,"$",sep = ""),files_list)]
req_files <- paste(maindir,samplename,"/outs/per_sample_outs/",samplename,"/antigen_analysis/antigen_specificity_scores.csv",sep = "")
Ag <- unique(read.csv(req_files[1], header = TRUE)[,2])
# rownames(df) <- Ag

rm(df_list)
df_list <- list()
for (i in 1:length(req_files)) {
  Ag_file <- read.csv(req_files[i], header = TRUE)
  Ag_high_score <- Ag_file[Ag_file$antigen_specificity_score > 50,]
  Ag_high_score$barcode_antigen <- paste(Ag_high_score$barcode, Ag_high_score$antigen, sep = "_")
  df_list[[samplename[i]]] <- Ag_high_score$barcode_antigen
}

df_Exp03 <- as.data.frame(unlist(df_list))

df_bind <- rbind(df_Exp01,df_Exp02,df_Exp03)
df_bind$sample_id <- gsub("Exp01.*.","Exp01",rownames(df_bind)) %>% gsub("Exp02.*.","Exp02",.) %>% gsub("Exp03.*.", "Exp03",.)
colnames(df_bind) <- c("barcode_Ag","sample_id")

df_bind$barcode <- gsub("_EBV_.*.|_VZV_.*.","",df_bind$barcode_Ag)
df_bind$sample_id_barcode <- paste(df_bind$sample_id, df_bind$barcode, sep = "_")
write.table(df_bind, paste(savedir,"Table/Antigen_positive_cells.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

#### Clones for which Antigen has been identified
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01_First_4/post_processing/"
exp <- "Exp01"
files_list <- list.files(maindir)
samplename <- files_list[grep(paste("_",exp,"$",sep = ""),files_list)]
req_files <- paste(maindir,samplename,"/outs/per_sample_outs/",samplename,"/antigen_analysis/antigen_specificity_scores.csv",sep = "")
Ag <- unique(read.csv(req_files[1], header = TRUE)[,2])
# rownames(df) <- Ag

rm(df_list)
df_list <- list()
for (i in 1:length(req_files)) {
  Ag_file <- read.csv(req_files[i], header = TRUE)
  Ag_high_score <- Ag_file[Ag_file$antigen_specificity_score > 50,]
  df_list[[samplename[i]]] <- Ag_high_score[grep("None",Ag_high_score$raw_clonotype_id, invert = TRUE),"raw_clonotype_id"]
}

df_Exp01 <- as.data.frame(unlist(df_list))

### Exp02
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02_First_4/post_processing/"
exp <- "Exp02"
files_list <- list.files(maindir)
samplename <- files_list[grep(paste("_",exp,"$",sep = ""),files_list)]
req_files <- paste(maindir,samplename,"/outs/per_sample_outs/",samplename,"/antigen_analysis/antigen_specificity_scores.csv",sep = "")
Ag <- unique(read.csv(req_files[1], header = TRUE)[,2])
# rownames(df) <- Ag

rm(df_list)
df_list <- list()
for (i in 1:length(req_files)) {
  Ag_file <- read.csv(req_files[i], header = TRUE)
  Ag_high_score <- Ag_file[Ag_file$antigen_specificity_score > 50,]
  df_list[[samplename[i]]] <- Ag_high_score[grep("None",Ag_high_score$raw_clonotype_id, invert = TRUE),"raw_clonotype_id"]
}

df_Exp02 <- as.data.frame(unlist(df_list))

## Exp03
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/"
exp <- "Exp03"
files_list <- list.files(maindir)
samplename <- files_list[grep(paste("_",exp,"$",sep = ""),files_list)]
req_files <- paste(maindir,samplename,"/outs/per_sample_outs/",samplename,"/antigen_analysis/antigen_specificity_scores.csv",sep = "")
Ag <- unique(read.csv(req_files[1], header = TRUE)[,2])
# rownames(df) <- Ag

rm(df_list)
df_list <- list()
for (i in 1:length(req_files)) {
  Ag_file <- read.csv(req_files[i], header = TRUE)
  Ag_high_score <- Ag_file[Ag_file$antigen_specificity_score > 50,]
  df_list[[samplename[i]]] <- Ag_high_score[grep("None",Ag_high_score$raw_clonotype_id, invert = TRUE),"raw_clonotype_id"]
}

df_Exp03 <- as.data.frame(unlist(df_list))

### Find out the BEAM-T positive cells 
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/"
exp <- "Exp03"
files_list <- list.files(maindir)
samplename <- files_list[grep(paste("_",exp,"$",sep = ""),files_list)]
req_files <- paste(maindir,samplename,"/outs/per_sample_outs/",samplename,"/antigen_analysis/antigen_specificity_scores.csv",sep = "")

rm(df_list)
df_list <- list()

for (j in 1:length(req_files)) {
  Ag_spec <- read.csv(req_files[j])
  barcodes <- unique(Ag_spec$barcode)
  
  ## Creating an empty dataframe for all Ag information
  rm(df)
  df <- data.frame(matrix(ncol = 7, nrow = length(barcodes)))
  colnames(df) <- c("samplename","barcode","Ag_UMI","Con_UMI","Ag_Con_UMI","max_specificity_score","max_Ag")
  df$samplename <- samplename[j]
  df$barcode <- barcodes
  
  for (i in 1:length(barcodes)) {
    subset_table <- Ag_spec[grep(barcodes[i], Ag_spec$barcode),c("antigen","antigen_umi","control_umi","antigen_specificity_score")]
    df[i,"Ag_UMI"] <- sum(subset_table$antigen_umi)
    df[i,"Con_UMI"] <- sum(subset_table$control_umi)
    df[i,"Ag_Con_UMI"] <- sum(subset_table$antigen_umi,subset_table$control_umi)
    df[i,"max_specificity_score"] <- max(subset_table[,"antigen_specificity_score"])
    df[i,"max_Ag"] <- paste(subset_table[grep(df[i,"max_specificity_score"], subset_table[,"antigen_specificity_score"]),"antigen"], collapse = ", ")
  }
  
  df_list[[samplename[j]]] <- df
}

combined_df_Exp03 <- bind_rows(df_list)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/Table/"
write.table(combined_df_Exp03, paste(savedir,"combined_df_Exp03.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

# head -1 combined_df_Exp03.txt > header
# cat combined_df_Exp01.txt combined_df_Exp02.txt combined_df_Exp03.txt | grep -v "barcode" > combined_Ag_Exp01_02_03.txt
# cat header combined_Ag_Exp01_02_03.txt > combined_Ag_Exp01_02_03_wid_header.txt

### Identifying the bulk and BEAM-T positive cells
## The cut-off that we are using is Ag+Con > 70 and Ag specificity score > 50 for the  
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/"
Ag_bulk_table <- read.table(paste(savedir,"Table/combined_Ag_Exp01_02_03_wid_header.txt",sep = ""), header = TRUE, sep = "\t")
Ag_gt_50 <- Ag_bulk_table[Ag_bulk_table$max_specificity_score > 50,]
Ag_not_gt_50 <- Ag_bulk_table[Ag_bulk_table$max_specificity_score <= 50,]
bulk_table <- Ag_not_gt_50[Ag_not_gt_50$Ag_Con_UMI <= 70,]

write.table(bulk_table, paste(savedir,"Table/Exp01_02_03_bulk_table.txt",sep = ""), sep = "\t", row.names = T, col.names = T, quote = F)
write.table(Ag_gt_50, paste(savedir,"Table/Exp01_02_03_Ag_table.txt",sep = ""),sep = "\t", row.names = T, col.names = T, quote = F)
write.table(Ag_not_gt_50, paste(savedir,"Table/Exp01_02_03_unspecified_table.txt",sep = ""),sep = "\t", row.names = T, col.names = T, quote = F)

### Preprocessing ####
library("Matrix")
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/CR7_scCITESeq_QC.R")

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/"
Exp01 = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01_First_4/post_processing/"
Exp02 = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02_First_4/post_processing/"
Exp03 = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/"

Exp01_files <- list.files(Exp01, pattern = "Exp01")
Exp01_files_2 <- Exp01_files[grep(".csv",Exp01_files, ignore.case = T, invert = T)]
Exp01_2 <- vector()
for (i in 1:length(Exp01_files_2)) {
  Exp01_2[i] <- paste(Exp01,Exp01_files_2[i],"/outs/per_sample_outs/",Exp01_files_2[i],"/",sep = "")
}

Exp02_files <- list.files(Exp02, pattern = "Exp02")
Exp02_files_2 <- Exp02_files[grep(".csv",Exp02_files, ignore.case = T, invert = T)]
Exp02_2 <- vector()
for (i in 1:length(Exp02_files_2)) {
  Exp02_2[i] <- paste(Exp02,Exp02_files_2[i],"/outs/per_sample_outs/",Exp02_files_2[i],"/",sep = "")
}

Exp03_files <- list.files(Exp03, pattern = "Exp03")
Exp03_files_2 <- Exp03_files[grep(".csv",Exp03_files, ignore.case = T, invert = T)]
Exp03_2 <- vector()
for (i in 1:length(Exp03_files_2)) {
  Exp03_2[i] <- paste(Exp03,Exp03_files_2[i],"/outs/per_sample_outs/",Exp03_files_2[i],"/",sep = "")
}

Exp <- c(Exp01_2, Exp02_2, Exp03_2)

id =  basename(Exp)
samplepath = Exp[which(basename(Exp) %in% id)]
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

mat <- matrix(nrow = length(samplename))
df <- as.data.frame(mat)
df$V1 <- table2$Doublets
rownames(df) <- samplename
colnames(df) <- "doublets"
table3 <- merge(table,df,by=0)
table3 <- table3[,-3]

dir.create(paste(savedir,"counts/",sep = ""), showWarnings = FALSE)
write.table(table3, paste(savedir,"counts/Doublets_Removed.txt",sep = ""), sep = "\t", quote = F, row.names = F, col.names = T)

#### Adding the unwanted effect i.e. cell cycle phase and mitochondria effect to each samples ####
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cellcycle_mitoscore.R")
obj1 <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "_RNA_singlet.RDS",full.names = TRUE)

RNA_RDS <- obj1
samplename <- gsub("_RNA_singlet.RDS","",basename(RNA_RDS))

combined <- merge(x=readRDS(RNA_RDS[1]), y = c(readRDS(RNA_RDS[2]),readRDS(RNA_RDS[3]),readRDS(RNA_RDS[4]),readRDS(RNA_RDS[5]),readRDS(RNA_RDS[6]),
                                               readRDS(RNA_RDS[7]),readRDS(RNA_RDS[8]),readRDS(RNA_RDS[9]),readRDS(RNA_RDS[10]),readRDS(RNA_RDS[11]),
                                               readRDS(RNA_RDS[12]),readRDS(RNA_RDS[13]),readRDS(RNA_RDS[14])),
                  add.cell.ids = samplename, project = "combined")

obj_ADT <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "_ADT_singlet.RDS",full.names = TRUE)
ADT_RDS <- c(obj_ADT)
# samplename_ADT <- paste(samplename,"ADT",sep = "_")

## We have to combine the RNA singlet object to perform the normalization for the unwanted the effects.
combined_ADT <- merge(x=readRDS(ADT_RDS[1]), y = c(readRDS(ADT_RDS[2]),readRDS(ADT_RDS[3]),readRDS(ADT_RDS[4]),readRDS(ADT_RDS[5]),readRDS(ADT_RDS[6]),
                                                   readRDS(ADT_RDS[7]),readRDS(ADT_RDS[8]),readRDS(ADT_RDS[9]),readRDS(ADT_RDS[10]),readRDS(ADT_RDS[11]),
                                                   readRDS(ADT_RDS[12]),readRDS(ADT_RDS[13]),readRDS(ADT_RDS[14])),
                      add.cell.ids = samplename, project = "combined_ADT")

combined_ADT <- NormalizeData(combined_ADT, normalization.method = 'CLR', margin = 2)
combined_ADT <- FindVariableFeatures(combined_ADT, selection.method = "vst", nfeatures = 36)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cellcycle_mitoscore.R")
samplename = "combined"
process = "merged"
Assay = "RNA"
combined <- cellcycle_mito(combined, savedir, ngenes = 4000, process = process, Assay = Assay, samplename = samplename)

combined@meta.data$samplename <- gsub("_[A|T|G|C].*","",rownames(combined@meta.data))
combined@meta.data$sample_id = gsub("_.*","",combined@meta.data$orig.ident)
combined@meta.data$Run = gsub(".*_","",combined@meta.data$orig.ident)
combined@meta.data$Age <- gsub("[0-9]","",combined@meta.data$sample_id)
combined@meta.data$gender <- gsub(".*._M_.*.","M",combined@meta.data$orig.ident) %>% gsub(".*._F_.*.","F",.)
combined@meta.data$Age_number <- gsub(".*._EBV_","",combined@meta.data$orig.ident) %>% gsub("_.*.","",.)

saveRDS(combined,paste(savedir,"saveRDS_obj/combined.RDS",sep = ""))
saveRDS(combined_ADT,paste(savedir,"saveRDS_obj/combined_ADT.RDS",sep = ""))

# combined = readRDS(paste(savedir,"saveRDS_obj/combined.RDS",sep = ""))
# combined_ADT = readRDS(paste(savedir,"saveRDS_obj/combined_ADT.RDS",sep = ""))
## In here we are extracting  CD8 positive T cells. 
## We have decided that we are going to remove the double negative and double positive population from CD4 while include it in the CD8
## Ines BEAM-T data
# combined = readRDS(paste(savedir,"saveRDS_obj/combined.RDS",sep = ""))
# combined_ADT = readRDS(paste(savedir,"saveRDS_obj/combined_ADT.RDS",sep = ""))

GEX_CD4_CD8 <- as.data.frame(t(as.data.frame(combined@assays$RNA@data[c("CD4","CD8A","CD8B"),])))
ADT_CD4_CD8 <- as.data.frame(t(as.data.frame(combined_ADT@assays$ADT@data[c("CD4-protein","CD8a-protein"),])))
colnames(ADT_CD4_CD8) <- c("CD4_protein","CD8a_protein")
all(rownames(combined@meta.data) == rownames(ADT_CD4_CD8))
ADT_CD4_CD8$Run <- combined@meta.data$Run
all(rownames(combined@meta.data) == rownames(GEX_CD4_CD8))
GEX_CD4_CD8$Run <- combined@meta.data$Run

library(ggplot2)
p1 <- ggplot(ADT_CD4_CD8, aes(x=CD4_protein, y=CD8a_protein, color = Run)) +
  geom_point(alpha=0.2) + theme_bw() +
  geom_hline(yintercept=0.5, linetype="dashed", color = "red") +
  geom_vline(xintercept=1, linetype="dashed", color = "red") +
  ggtitle("All ADT")

p2 <- ggplot(GEX_CD4_CD8, aes(x=CD4,y=CD8A, color = Run)) +
  geom_point() + theme_bw() +
  geom_hline(yintercept=0.2, linetype="dashed", color = "red") +
  geom_vline(xintercept=0.2, linetype="dashed", color = "red") +
  ggtitle("All GEX CD8A")

p3 <- ggplot(GEX_CD4_CD8, aes(x=CD4,y=CD8B, color = Run)) +
  geom_point() + theme_bw() +
  geom_hline(yintercept=0.2, linetype="dashed", color = "red") +
  geom_vline(xintercept=0.2, linetype="dashed", color = "red") +
  ggtitle("All GEX CD8B")


dir.create(paste(savedir,"scatter_plot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"scatter_plot/ADT_and_GEX_CD4_and_CD8.pdf",sep = ""), width = 15, height = 6)
print(p1+p2+p3)
dev.off()

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/"
# combined = readRDS(paste(savedir,"saveRDS_obj/combined.RDS",sep = ""))
# combined_ADT = readRDS(paste(savedir,"saveRDS_obj/combined_ADT.RDS",sep = ""))
# CD4_cellnames <- read.table(paste(savedir,"Run01_02_CD4_cellnames.txt",sep = ""))[,1]

### RNA #####
library(Seurat)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_subsetting.R")
RNA_obj_path = combined
ADT_obj_path = combined_ADT
ADT_main = combined_ADT
CD8_cluster <- NULL
CD8_cellnames <- NULL
objname = "CD8"
Assay = "RNA"
process = "subset"

CD8_RNA <- RNA_subsetting(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path, ADT_main = ADT_main, subset_cluster = CD8_cluster,
                          cellnames=CD8_cellnames, ngenes = 4000, saveDir = savedir, objname = objname, Assay = Assay, process = process)

CD8_RNA <- RunUMAP(CD8_RNA, dims = 1:20, reduction = "pca")

dir.create(paste(savedir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/CD8_no_integration_samples.pdf",sep = ""),width = 10, height = 5)
print(DimPlot(CD8_RNA, group.by = "orig.ident"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_no_integration_samples_splitted.pdf",sep = ""),width = 15, height = 15)
print(DimPlot(CD8_RNA, split.by = "orig.ident", group.by="orig.ident", ncol = 4) + NoLegend())
dev.off()

dir.create(paste(savedir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/CD8_no_integration_Run.pdf",sep = ""),width = 6, height = 5)
print(DimPlot(CD8_RNA, group.by = "Run"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_no_integration_Run_splitted.pdf",sep = ""),width = 15, height = 5)
print(DimPlot(CD8_RNA, split.by = "Run", group.by="Run", ncol = 3) + NoLegend())
dev.off()

## Since there is a batch effect we have integrate based on the Run
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_sctransform_V2.R")
objname = "CD8"
Assay = "RNA"
process = "sctransform"
CD8_RNA_sctransformed <- sctransform_V2_integration(obj = CD8_RNA, saveDir = savedir, ngenes = 4000,
                                                    regress = c("nCount_RNA"),
                                                    dims = 20,
                                                    Assay = Assay, process = process, objname = objname,
                                                    split_by = "Run",
                                                    reference = NULL,
                                                    sample_tree = NULL)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "CD8"
process = "integration"
Assay = "RNA"
CD8_integrated_RNA <- RNA_integration(CD8_RNA_sctransformed, savedir, dims = 20, RNA_features = c("CD4","CD8A"),
                                      Assay=Assay, process=process, objname=objname, ncol = 4, ndims = 50)

CD8_integrated_RNA_NN <- FindNeighbors(CD8_integrated_RNA, dims = 1:20)

CD8_integrated_RNA_NN@meta.data$samplename <- gsub("_[A|T|G|C].*","",rownames(CD8_integrated_RNA_NN@meta.data))
CD8_integrated_RNA_NN@meta.data$sample_id = gsub("_.*","",CD8_integrated_RNA_NN@meta.data$orig.ident)
CD8_integrated_RNA_NN@meta.data$Run = gsub(".*_","",CD8_integrated_RNA_NN@meta.data$orig.ident)
CD8_integrated_RNA_NN@meta.data$Age <- gsub("[0-9]","",CD8_integrated_RNA_NN@meta.data$sample_id)
CD8_integrated_RNA_NN@meta.data$gender <- gsub(".*._M_.*.","M",CD8_integrated_RNA_NN@meta.data$orig.ident) %>% gsub(".*._F_.*.","F",.)
CD8_integrated_RNA_NN@meta.data$Age_number <- gsub(".*._EBV_","",CD8_integrated_RNA_NN@meta.data$orig.ident) %>% gsub("_.*.","",.)


pdf(paste(savedir,"UMAP/CD8_RNA_Age.pdf",sep = ""))
DimPlot(CD8_integrated_RNA_NN, group.by = "Age")
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_batch.pdf",sep = ""))
DimPlot(CD8_integrated_RNA_NN, group.by = "Run")
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_Age_split.pdf",sep = ""), width = 10, height = 6)
DimPlot(CD8_integrated_RNA_NN, group.by = "Age", split.by = "Age", ncol = 2)
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_batch_split.pdf",sep = ""), width = 15, height = 6)
DimPlot(CD8_integrated_RNA_NN, group.by = "Run", split.by = "Run", ncol = 3)
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_nCount_RNA_SCT.pdf",sep = ""), width = 10, height = 6)
FeaturePlot(CD8_integrated_RNA_NN, c("nCount_SCT", "nCount_RNA"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_nCount_RNA_SCT_Vln.pdf",sep = ""), width = 10, height = 6)
VlnPlot(CD8_integrated_RNA_NN, c("nCount_SCT", "nCount_RNA"), group.by = "orig.ident")
dev.off()

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8,1)
res = 0.4
for (i in 1:length(res)) {
  process = paste("RNA_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  CD8_integrated_RNA_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD8_integrated_RNA_NN, dims = 30, res = res[i], 
                                                       saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_RNA","nFeature_RNA"), 
                                                       objname = "CD8_integrated_RNA_NN_cluster",
                                                       process = process, col_sel = c("Age","Run","orig.ident","gender"))
}

saveRDS(CD8_integrated_RNA_NN_cluster,paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster.RDS",sep = ""))

## Since Magic does not work on the conda seurat_new we have to perform the analysis outside conda.
## ml load r/4.2.2
library(Seurat)
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/"
CD8_integrated_RNA_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster.RDS",sep = ""))
DefaultAssay(CD8_integrated_RNA_NN_cluster) <- "RNA"
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/miniconda3/envs/seurat_new/bin/python3.8")
library(Rmagic)
py_discover_config("magic") # to check
CD8_integrated_RNA_NN_cluster_impute <- magic(CD8_integrated_RNA_NN_cluster, npca=30) ## imputing the RNA data as for RNA PCs are 20
saveRDS(CD8_integrated_RNA_NN_cluster_impute,paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster_impute.RDS",sep = ""))

DefaultAssay(CD8_integrated_RNA_NN_cluster) <- "RNA"
Y3_Y6_cluster <- FindMarkers(CD8_integrated_RNA_NN_cluster, ident.1 = c(7,10,25,26,27), only.pos = FALSE)
Y3_Y6_cluster_between_naive <- FindMarkers(CD8_integrated_RNA_NN_cluster, ident.1 = c(10,25,26,27), ident.2 = c(1,20,22,23), only.pos = FALSE)
Y3_Y6_no_clus7_cluster_between_naive <- FindMarkers(CD8_integrated_RNA_NN_cluster, ident.1 = c(10,25,26,27), ident.2 = c(1,20,22,23), only.pos = FALSE)
Y3_Y6_cluster_no_clus7 <- FindMarkers(CD8_integrated_RNA_NN_cluster, ident.1 = c(10,25,26,27), only.pos = FALSE)

Y3_Y6_cluster_no_clus7_ordered <- Y3_Y6_cluster_no_clus7[order(-Y3_Y6_cluster_no_clus7$avg_log2FC),]

genes <- rownames(Y3_Y6_cluster_no_clus7_ordered)[c(1:10,188:196)]
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p <- featureplot_front(CD8_integrated_RNA_NN_cluster, genes[i],
                         reduction = "umap", x="UMAP_1", y="UMAP_2",size=0.2) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- p
}

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/Y3_Y6_high_expression_no_clus_7.pdf",sep = ""))
print(plot_list)
dev.off()

library(scales)
CD8 <- CD8_integrated_RNA_NN_cluster
dir.create(paste(savedir,"UMAP/clusters/",sep = ""), showWarnings = FALSE)
clusters = 0:27
CD8@meta.data$seurat_clusters <- factor(CD8@meta.data$seurat_clusters, levels = clusters)
rm(plot_list)
plot_list <- list()
for (i in 1:length(clusters)) {
  col <- rep("white",length(clusters))
  col[i] <- "red"
    p <- DimPlot(CD8, group.by = "seurat_clusters", label = TRUE, cols=alpha(col,0.3), reduction = "umap")
    pdf(paste(savedir,"UMAP/clusters/CD8_",clusters[i],"_UMAP.pdf",sep = ""))
    print(p)
    dev.off()
    plot_list[[i]] <- p
}

dir.create(paste(savedir,"UMAP/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"UMAP/clusters/CD8_all_combine.pdf",sep = ""), width = 4.5, height = 4)
plot_list
dev.off()



### CD8 ADT #######
CD8_ADT = readRDS(paste(savedir,"saveRDS_obj/CD8_subset_ADT.RDS",sep = ""))

CD8_ADT@meta.data$samplename <- gsub("_[A|T|G|C].*","",rownames(CD8_ADT@meta.data))
CD8_ADT@meta.data$sample_id = gsub("_.*","",CD8_ADT@meta.data$orig.ident)
CD8_ADT@meta.data$Run = gsub(".*_","",CD8_ADT@meta.data$orig.ident)
CD8_ADT@meta.data$Age <- gsub("[0-9]","",CD8_ADT@meta.data$sample_id)
CD8_ADT@meta.data$gender <- gsub(".*._M_.*.","M",CD8_ADT@meta.data$orig.ident) %>% gsub(".*._F_.*.","F",.)
CD8_ADT@meta.data$Age_number <- gsub(".*._EBV_","",CD8_ADT@meta.data$orig.ident) %>% gsub("_.*.","",.)

### No integration
CD8_ADT <- NormalizeData(CD8_ADT, normalization.method = 'CLR', margin = 2)
CD8_ADT <- FindVariableFeatures(CD8_ADT, selection.method = "vst", nfeatures = 30)
CD8_ADT <- ScaleData(CD8_ADT, verbose = FALSE)
CD8_ADT <- RunPCA(CD8_ADT, npcs = 30, approx=FALSE)
CD8_ADT <- RunUMAP(CD8_ADT, dims = 1:10)

dir.create(paste(savedir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/CD8_ADT_no_integration_samples.pdf",sep = ""),width = 10, height = 5)
print(DimPlot(CD8_ADT, group.by = "orig.ident"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_no_integration_samples_splitted.pdf",sep = ""),width = 15, height = 15)
print(DimPlot(CD8_ADT, split.by = "orig.ident", group.by="orig.ident", ncol = 4) + NoLegend())
dev.off()

dir.create(paste(savedir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/CD8_ADT_no_integration_Run.pdf",sep = ""),width = 6, height = 5)
print(DimPlot(CD8_ADT, group.by = "Run"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_no_integration_Run_splitted.pdf",sep = ""),width = 15, height = 5)
print(DimPlot(CD8_ADT, split.by = "Run", group.by="Run", ncol = 3) + NoLegend())
dev.off()

DefaultAssay(CD8_ADT) <- "ADT"
dir.create(paste(savedir,"heatmap",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"heatmap/CD8_ADT.pdf",sep =""))
DoHeatmap(CD8_ADT, features = rownames(CD8_ADT))
dev.off()

DefaultAssay(CD8_ADT) <- "ADT"
dir.create(paste(savedir,"heatmap",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"heatmap/CD8_ADT_Run.pdf",sep =""))
DoHeatmap(CD8_ADT, features = rownames(CD8_ADT),group.by = "Run")
dev.off()

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")
objname = "CD8"
process = "integrated"
Assay = "ADT"

CD8_ADT_integrated <- ADT_merging(CD8_ADT, savedir, dims = 6, numfeatures=36,
                                  Assay=Assay, process=process, objname=objname,
                                  sample_tree = NULL, split_by = "Run",
                                  reference=NULL)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "CD8"
process = "integration"
Assay = "ADT"
CD8_ADT_integrated <- RNA_integration(obj_path = CD8_ADT_integrated, saveDir = savedir, dims = 10,
                                      RNA_features = c("CD8-protein","CD8a-protein"),Assay = Assay,
                                      process = process, objname = objname, ncol = 4)

### Nearest Neighbour and Clustering
CD8_ADT_integrated_NN <- FindNeighbors(CD8_ADT_integrated, dims = 1:10)

pdf(paste(savedir,"UMAP/CD8_ADT_Age.pdf",sep = ""))
DimPlot(CD8_ADT_integrated_NN, group.by = "Age")
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_batch.pdf",sep = ""))
DimPlot(CD8_ADT_integrated_NN, group.by = "Run")
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_Age_split.pdf",sep = ""), width = 10, height = 6)
DimPlot(CD8_ADT_integrated_NN, group.by = "Age", split.by = "Age", ncol = 2)
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_batch_split.pdf",sep = ""), width = 15, height = 6)
DimPlot(CD8_ADT_integrated_NN, group.by = "Run", split.by = "Run", ncol = 3)
dev.off()

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8)
res = 0.6
for (i in 1:length(res)) {
  process = paste("ADT_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  CD8_ADT_integrated_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD8_ADT_integrated_NN, dims = 10, res = res[i], saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_ADT","nFeature_ADT"), objname = "CD8_ADT_integrated_NN_cluster",
                                                       process = process,col_sel = c("Age","Run","orig.ident","gender"))
}

### Generating the MAGIC ADT and MAGIC RNA
saveRDS(CD8_ADT_integrated_NN_cluster,paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster.RDS",sep = ""))

## Since Magic does not work on the conda seurat_new we have to perform the analysis outside conda.
## ml load r/4.2.2
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/"
CD8_ADT_integrated_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster.RDS",sep = ""))
DefaultAssay(CD8_ADT_integrated_NN_cluster) <- "ADT"
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/miniconda3/envs/seurat_new/bin/python3.8")
library(Rmagic)
py_discover_config("magic") # to check
CD8_ADT_integrated_NN_cluster_impute <- magic(CD8_ADT_integrated_NN_cluster, npca=10) ## imputing the RNA data as for RNA PCs are 20
saveRDS(CD8_ADT_integrated_NN_cluster_impute,paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster_impute.RDS",sep = ""))

# CD8_integrated_RNA_NN_cluster <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run06/Downstream/saveRDS_obj/CD8_integrated_RNA_NN_cluster.RDS")
#### CD8 Modality Integration ######
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITEseq_modality_integration.R")
# Find the last save object for RNA and ADT
# we have RNA already in the path
# CD8_ADT_integrated_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster_impute.RDS",sep = ""))
# CD8_integrated_RNA_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster_impute.RDS",sep = ""))
RNA_obj_path <- CD8_integrated_RNA_NN_cluster_impute
ADT_obj_path <- CD8_ADT_integrated_NN_cluster_impute

objname <- "CD8"
process <- "modality_integrate"
Assay <- "integrated"

CD8_modality_integrate <- modality_integration(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path,
                                               RNA_dims = 20, ADT_dims = 10,
                                               saveDir = savedir, Assay = Assay, process = process, objname = objname)

DefaultAssay(CD8_modality_integrate) <- "integrated"

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/"
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/Modal_cluster_and_UMAP.R")
obj_path <- CD8_modality_integrate
RNA_features = c("CD4","CD8A")
ADT_features = c("CD4-protein","CD8a-protein")

# readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run03/analysis/downstream/saveRDS_obj/CD8_modality_integrate_cluster_1.RDS")
res = c(0.2,0.4,0.6,0.8,1)
res=0.2
for (i in 1:length(res)) {
  objname = "CD8"
  Assay = "integrated"
  process = paste("modality_cluster_and_QC",res[i],sep = "_")
  CD8_modality_integrate_cluster <- Modal_Cluster_and_UMAP_QC(modal_integrate_obj_path = obj_path, saveDir = savedir,
                                                              res = res[i], RNA_features = RNA_features,
                                                              ADT_features = ADT_features,Assay = Assay, protein_assay = "ADT",
                                                              process = process, objname=objname,
                                                              col_sel=c("Age","Run","orig.ident","gender"))
}

DefaultAssay(CD8_modality_integrate_cluster) <- "RNA"
clus6_17 <- FindMarkers(CD8_modality_integrate_cluster, ident.1 =c(6,17), only.pos = TRUE)

saveRDS(CD8_modality_integrate_cluster,paste(savedir,"saveRDS_obj/CD8_modality_integrate_cluster_0.2.RDS",sep = ""))

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/ADT_featureplot.R")
ADT_featureplot(obj = CD8_modality_integrate_cluster, savedir = savedir,objname = "CD8_modality_integrate_cluster",
                process = "featureplots",reduction = "wnn.umap",x = "wnnUMAP_1",y = "wnnUMAP_2")

CD8_samplewise_counts <- table(CD8_modality_integrate_cluster@meta.data$orig.ident) %>% as.data.frame()
write.table(CD8_samplewise_counts, paste(savedir,"Table/CD8_samplewise_counts.txt",sep = ""), sep = "\t",
            row.names = F, col.names = T, quote = F)

Ag_pos_cells <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Table/Antigen_positive_cells.txt", header = TRUE)
Ag_pos_cells$Ag <- gsub(".*.-1_","",Ag_pos_cells$barcode_Ag)

CD8_modality_integrate_cluster@meta.data$Ag_gt_50 <- Ag_pos_cells[match(rownames(CD8_modality_integrate_cluster@meta.data),Ag_pos_cells$sample_id_barcode),"Ag"]

CD8_modality_integrate_cluster@meta.data$Ag_Age <- paste(CD8_modality_integrate_cluster@meta.data$Ag_gt_50, CD8_modality_integrate_cluster@meta.data$Age,sep = "_")

Ag_Age <- unique(CD8_modality_integrate_cluster@meta.data$Ag_Age)
Ag_Age2 <- Ag_Age[grep("NA_O|NA_Y",Ag_Age,invert=TRUE)]

rm(plot_list)
plot_list <- list()
for (i in 1:length(Ag_Age2)) {
  cellnames <- CD8_modality_integrate_cluster@meta.data[grep(Ag_Age2[i], CD8_modality_integrate_cluster@meta.data$Ag_Age),] %>% rownames()
  p <- DimPlot(CD8_modality_integrate_cluster,
               cells.highlight = cellnames, 
               reduction = "wnn.umap", 
               label = FALSE, cols.highlight = "deeppink2",
               sizes.highlight = 0.5,
               cols = "gray92") + 
    ggtitle(paste(Ag_Age2[i])) + 
    NoLegend() + 
    theme(plot.title = element_text(hjust = 0.2))
  # pdf(paste(savedir,"CD4/CD4_TRB_",TRB_clonal_noNA_gt_2$Var1[i],".pdf",sep = ""))
  # print(p)
  # dev.off()
  plot_list[[i]] <- p
}

pdf(paste(savedir,"UMAP/Ag_gt_50_Age.pdf",sep = ""), width = 5.5, height = 5)
plot_list
dev.off()

CD8_modality_integrate_cluster@meta.data$Ag_group <- NA
CD8_modality_integrate_cluster@meta.data[grep("EBV_B",CD8_modality_integrate_cluster@meta.data$Ag_gt_50),"Ag_group"] <- "Lytic"
CD8_modality_integrate_cluster@meta.data[grep("EBV_EB|EBV_LM",CD8_modality_integrate_cluster@meta.data$Ag_gt_50),"Ag_group"] <- "Latent"
CD8_modality_integrate_cluster@meta.data$Ag_group_Age <- paste(CD8_modality_integrate_cluster@meta.data$Ag_group,CD8_modality_integrate_cluster$Age,sep="_")
group_Ag <- grep("NA_",unique(CD8_modality_integrate_cluster@meta.data$Ag_group_Age),invert = TRUE,value = TRUE)

rm(plot_list)
plot_list <- list()
for (i in 1:length(group_Ag)) {
  cellnames <- CD8_modality_integrate_cluster@meta.data[grep(group_Ag[i], CD8_modality_integrate_cluster@meta.data$Ag_group_Age),] %>% rownames()
  p <- DimPlot(CD8_modality_integrate_cluster,
               cells.highlight = cellnames,
               reduction = "wnn.umap",
               label = FALSE, cols.highlight = "deeppink2",
               sizes.highlight = 0.5,
               cols = "gray92") + 
    ggtitle(paste(group_Ag[i])) + 
    NoLegend() + 
    theme(plot.title = element_text(hjust = 0.2))
  # pdf(paste(savedir,"CD4/CD4_TRB_",TRB_clonal_noNA_gt_2$Var1[i],".pdf",sep = ""))
  # print(p)
  # dev.off()
  plot_list[[i]] <- p
}

pdf(paste(savedir,"UMAP/Ag_group_gt_50_Age.pdf",sep = ""), width = 5.5, height = 5)
plot_list
dev.off()

pdf(paste(savedir,"UMAP/Ag_group_gt_50_Age.pdf",sep = ""), width = 11, height = 5)
DimPlot(CD8_modality_integrate_cluster, split.by = "Age", group.by = "Age", ncol = 2, reduction = "wnn.umap")
dev.off()


# CD4_genes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run01/Analysis/downstream/Table/CD4_genes")[,1]
# dir.create(paste(savedir,"dotplot",sep = ""),showWarnings = FALSE)
# pdf(paste(savedir,"dotplot/CD4_genes_CD4_dotplot.pdf",sep = ""))
# DotPlot(CD4_modality_integrate_cluster,"RNA",CD4_genes) + coord_flip()
# dev.off()
# 
# CD8_genes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run01/Analysis/downstream/Table/CD8_genes")[,1]
# dir.create(paste(savedir,"dotplot",sep = ""),showWarnings = FALSE)
# pdf(paste(savedir,"dotplot/CD8_genes_CD4_dotplot.pdf",sep = ""), width = 7, height = 9)
# DotPlot(CD4_modality_integrate_cluster,"RNA",CD8_genes) + coord_flip()
# dev.off()
# 
# ### exhauster genes
# ex_genes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run01/Analysis/downstream/Table/exhausted_genes_2")
# ex_genes_up <- toupper(ex_genes[,1])
# pdf(paste(savedir,"dotplot/exhausted_CD4_dotplot.pdf",sep = ""), height = 10, width = 8)
# DotPlot(CD4_modality_integrate_cluster,"RNA",ex_genes_up) + coord_flip()
# dev.off()
# 
# DefaultAssay(CD4_modality_integrate_cluster) <- "MAGIC_RNA"
# genes <- c("IKZF2","FOXP3","TCF7","LEF1","IRF2","RUNX3","EOMES","TBX21","PRDM1","XBP1",
#            "KLF6","NOTCH1","SLAMF6","FOSL2","FOS","RORC","JUN","JUND","JUNB","STAT3",
#            "Zeb2","ID2","Tox","Myb","TCF7","ID3","BCL6","myc","IL7R","CD84","IL2RB")
# genes <- toupper(genes)
# genes <- rownames(CD4_modality_integrate_cluster)[match(toupper(genes),rownames(CD4_modality_integrate_cluster))]

genes <- c("IL7R", "CCR7", "CD45RA", "CD45RO", "PD1", "TIM3")
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p <- featureplot_front(CD4_modality_integrate_cluster, genes[i],
                         reduction = "wnn.umap", x="wnnUMAP_1", y="wnnUMAP_2",size=0.2) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- p
}

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/CD4_Goronzy.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list)
dev.off()

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/Gene_score.R")
CD4_modality_integrate_cluster <- genescore_abhinav(obj = CD4_modality_integrate_cluster,assay = "RNA",
                                                    savedir = savedir,
                                                    cores = 5,objname = "CD4")

### Bulk RNA ####
library("Matrix")
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/"
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Bulk/regress_TCR/"
dir.create(savedir,showWarnings = FALSE)
bulk_table <- read.table(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Bulk/Table/Exp01_02_03_bulk_table.txt",sep = ""), sep = "\t", header = T)

# combined_ADT <- readRDS(paste(maindir, "saveRDS_obj/combined_ADT.RDS",sep = ""))
# combined <- readRDS(paste(maindir, "saveRDS_obj/combined.RDS",sep = ""))

# TCR_genes <- grep("^TRAV|^TRAJ|^TRBC|^TRBV|^TRBJ|^TRBC",rownames(CD8_modality_integrate_cluster@assays$RNA@counts),value=TRUE)
# write.table(TCR_genes, "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/TCR_genes.txt",
#             sep = "\t", row.names = F, col.names = F, quote = F)

remove_genes <- c("MTRNR2L8",grep("^TRAV|^TRAJ|^TRBC|^TRBV|^TRBJ|^TRBC",
                                  rownames(CD8_modality_integrate_cluster@assays$RNA@counts),
                                  value=TRUE))
keep_genes <- grep(paste(remove_genes,collapse = "|"),
                   rownames(CD8_modality_integrate_cluster@assays$RNA@counts), 
                   invert = TRUE, value = TRUE)
combined_no_TCR <- combined[keep_genes,]

### Just to check
grep("MTRNR2L8|^TRAV|^TRAJ|^TRBC|^TRBV|^TRBJ|^TRBC",
     rownames(combined_no_TCR@assays$RNA@counts),
     value=TRUE)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_subsetting.R")
RNA_obj_path = combined_no_TCR
ADT_obj_path = combined_ADT
ADT_main = combined_ADT
CD8_cluster <- NULL
CD8_cellnames <- bulk_cellnames <- paste(bulk_table$samplename,bulk_table$barcode, sep="_")
objname = "CD8"
Assay = "RNA"
process = "bulk_subset"

CD8_RNA <- RNA_subsetting(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path, ADT_main = ADT_main, subset_cluster = CD8_cluster,
                          cellnames=CD8_cellnames, ngenes = 4000, saveDir = savedir, objname = objname, Assay = Assay, process = process)

CD8_RNA <- RunUMAP(CD8_RNA, dims = 1:20, reduction = "pca")

dir.create(paste(savedir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/CD8_no_integration_samples.pdf",sep = ""),width = 10, height = 5)
print(DimPlot(CD8_RNA, group.by = "orig.ident"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_no_integration_samples_splitted.pdf",sep = ""),width = 15, height = 15)
print(DimPlot(CD8_RNA, split.by = "orig.ident", group.by="orig.ident", ncol = 4) + NoLegend())
dev.off()

dir.create(paste(savedir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/CD8_no_integration_Run.pdf",sep = ""),width = 6, height = 5)
print(DimPlot(CD8_RNA, group.by = "Run"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_no_integration_Run_splitted.pdf",sep = ""),width = 15, height = 5)
print(DimPlot(CD8_RNA, split.by = "Run", group.by="Run", ncol = 3) + NoLegend())
dev.off()

## Since there is a batch effect we have integrate based on the Run
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_sctransform_V2.R")
objname = "CD8"
Assay = "RNA"
process = "sctransform"
CD8_RNA_sctransformed <- sctransform_V2_integration(obj = CD8_RNA, saveDir = savedir, ngenes = 4000,
                                                    regress = c("nCount_RNA"),
                                                    dims = 20,
                                                    Assay = Assay, process = process, objname = objname,
                                                    split_by = "Run",
                                                    reference = NULL,
                                                    sample_tree = NULL)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "CD8"
process = "integration"
Assay = "RNA"
CD8_integrated_RNA <- RNA_integration(CD8_RNA_sctransformed, savedir, dims = 20, RNA_features = c("CD8B","CD8A"),
                                      Assay=Assay, process=process, objname=objname, ncol = 4, ndims = 50)

CD8_integrated_RNA_NN <- FindNeighbors(CD8_integrated_RNA, dims = 1:20)

CD8_integrated_RNA_NN@meta.data$samplename <- gsub("_[A|T|G|C].*","",rownames(CD8_integrated_RNA_NN@meta.data))
CD8_integrated_RNA_NN@meta.data$sample_id = gsub("_.*","",CD8_integrated_RNA_NN@meta.data$orig.ident)
CD8_integrated_RNA_NN@meta.data$Run = gsub(".*_","",CD8_integrated_RNA_NN@meta.data$orig.ident)
CD8_integrated_RNA_NN@meta.data$Age <- gsub("[0-9]","",CD8_integrated_RNA_NN@meta.data$sample_id)
CD8_integrated_RNA_NN@meta.data$gender <- gsub(".*._M_.*.","M",CD8_integrated_RNA_NN@meta.data$orig.ident) %>% gsub(".*._F_.*.","F",.)
CD8_integrated_RNA_NN@meta.data$Age_number <- gsub(".*._EBV_","",CD8_integrated_RNA_NN@meta.data$orig.ident) %>% gsub("_.*.","",.)


pdf(paste(savedir,"UMAP/CD8_RNA_Age.pdf",sep = ""))
DimPlot(CD8_integrated_RNA_NN, group.by = "Age")
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_batch.pdf",sep = ""))
DimPlot(CD8_integrated_RNA_NN, group.by = "Run")
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_Age_split.pdf",sep = ""), width = 10, height = 6)
DimPlot(CD8_integrated_RNA_NN, group.by = "Age", split.by = "Age", ncol = 2)
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_batch_split.pdf",sep = ""), width = 15, height = 6)
DimPlot(CD8_integrated_RNA_NN, group.by = "Run", split.by = "Run", ncol = 3)
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_nCount_RNA_SCT.pdf",sep = ""), width = 10, height = 6)
FeaturePlot(CD8_integrated_RNA_NN, c("nCount_SCT", "nCount_RNA"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_nCount_RNA_SCT_Vln.pdf",sep = ""), width = 10, height = 6)
VlnPlot(CD8_integrated_RNA_NN, c("nCount_SCT", "nCount_RNA"), group.by = "orig.ident")
dev.off()

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8,1)
res = 0.4
for (i in 1:length(res)) {
  process = paste("RNA_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  CD8_integrated_RNA_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD8_integrated_RNA_NN, dims = 30, res = res[i], 
                                                       saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_RNA","nFeature_RNA"), 
                                                       objname = "CD8_integrated_RNA_NN_cluster",
                                                       process = process, col_sel = c("Age","Run","orig.ident","gender"))
}

saveRDS(CD8_integrated_RNA_NN_cluster,paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster_0.4.RDS",sep = ""))

## Since Magic does not work on the conda seurat_new we have to perform the analysis outside conda.
## ml load r/4.2.2
library(Seurat)
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Bulk/regress_TCR/"
CD8_integrated_RNA_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster_0.4.RDS",sep = ""))
DefaultAssay(CD8_integrated_RNA_NN_cluster) <- "RNA"
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/miniconda3/envs/seurat_new/bin/python3.8")
library(Rmagic)
py_discover_config("magic") # to check
CD8_integrated_RNA_NN_cluster_impute <- magic(CD8_integrated_RNA_NN_cluster, npca=30) ## imputing the RNA data as for RNA PCs are 20
saveRDS(CD8_integrated_RNA_NN_cluster_impute,paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster_impute.RDS",sep = ""))

DefaultAssay(CD8_integrated_RNA_NN_cluster) <- "RNA"
Y3_Y6_cluster <- FindMarkers(CD8_integrated_RNA_NN_cluster, ident.1 = c(7,10,25,26,27), only.pos = FALSE)
Y3_Y6_cluster_between_naive <- FindMarkers(CD8_integrated_RNA_NN_cluster, ident.1 = c(10,25,26,27), ident.2 = c(1,20,22,23), only.pos = FALSE)
Y3_Y6_no_clus7_cluster_between_naive <- FindMarkers(CD8_integrated_RNA_NN_cluster, ident.1 = c(10,25,26,27), ident.2 = c(1,20,22,23), only.pos = FALSE)
Y3_Y6_cluster_no_clus7 <- FindMarkers(CD8_integrated_RNA_NN_cluster, ident.1 = c(10,25,26,27), only.pos = FALSE)

Y3_Y6_cluster_no_clus7_ordered <- Y3_Y6_cluster_no_clus7[order(-Y3_Y6_cluster_no_clus7$avg_log2FC),]

genes <- rownames(Y3_Y6_cluster_no_clus7_ordered)[c(1:10,188:196)]
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p <- featureplot_front(CD8_integrated_RNA_NN_cluster, genes[i],
                         reduction = "umap", x="UMAP_1", y="UMAP_2",size=0.2) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- p
}

dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"featureplot/Y3_Y6_high_expression_no_clus_7.pdf",sep = ""))
print(plot_list)
dev.off()

library(scales)
CD8 <- CD8_integrated_RNA_NN_cluster
dir.create(paste(savedir,"UMAP/clusters/",sep = ""), showWarnings = FALSE)
clusters = 0:27
CD8@meta.data$seurat_clusters <- factor(CD8@meta.data$seurat_clusters, levels = clusters)
rm(plot_list)
plot_list <- list()
for (i in 1:length(clusters)) {
  col <- rep("white",length(clusters))
  col[i] <- "red"
    p <- DimPlot(CD8, group.by = "seurat_clusters", label = TRUE, cols=alpha(col,0.3), reduction = "umap")
    pdf(paste(savedir,"UMAP/clusters/CD8_",clusters[i],"_UMAP.pdf",sep = ""))
    print(p)
    dev.off()
    plot_list[[i]] <- p
}

dir.create(paste(savedir,"UMAP/",sep =""),showWarnings = FALSE)
pdf(paste(savedir,"UMAP/clusters/CD8_all_combine.pdf",sep = ""), width = 4.5, height = 4)
plot_list
dev.off()

### CD8 ADT #######
CD8_ADT = readRDS(paste(savedir,"saveRDS_obj/CD8_bulk_subset_ADT.RDS",sep = ""))

CD8_ADT@meta.data$samplename <- gsub("_[A|T|G|C].*","",rownames(CD8_ADT@meta.data))
CD8_ADT@meta.data$sample_id = gsub("_.*","",CD8_ADT@meta.data$orig.ident)
CD8_ADT@meta.data$Run = gsub(".*_","",CD8_ADT@meta.data$orig.ident)
CD8_ADT@meta.data$Age <- gsub("[0-9]","",CD8_ADT@meta.data$sample_id)
CD8_ADT@meta.data$gender <- gsub(".*._M_.*.","M",CD8_ADT@meta.data$orig.ident) %>% gsub(".*._F_.*.","F",.)
CD8_ADT@meta.data$Age_number <- gsub(".*._EBV_","",CD8_ADT@meta.data$orig.ident) %>% gsub("_.*.","",.)

### No integration
CD8_ADT <- NormalizeData(CD8_ADT, normalization.method = 'CLR', margin = 2)
CD8_ADT <- FindVariableFeatures(CD8_ADT, selection.method = "vst", nfeatures = 30)
CD8_ADT <- ScaleData(CD8_ADT, verbose = FALSE)
CD8_ADT <- RunPCA(CD8_ADT, npcs = 30, approx=FALSE)
CD8_ADT <- RunUMAP(CD8_ADT, dims = 1:10)

dir.create(paste(savedir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/CD8_ADT_no_integration_samples.pdf",sep = ""),width = 10, height = 5)
print(DimPlot(CD8_ADT, group.by = "orig.ident"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_no_integration_samples_splitted.pdf",sep = ""),width = 15, height = 15)
print(DimPlot(CD8_ADT, split.by = "orig.ident", group.by="orig.ident", ncol = 4) + NoLegend())
dev.off()

dir.create(paste(savedir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/CD8_ADT_no_integration_Run.pdf",sep = ""),width = 6, height = 5)
print(DimPlot(CD8_ADT, group.by = "Run"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_no_integration_Run_splitted.pdf",sep = ""),width = 15, height = 5)
print(DimPlot(CD8_ADT, split.by = "Run", group.by="Run", ncol = 3) + NoLegend())
dev.off()

DefaultAssay(CD8_ADT) <- "ADT"
dir.create(paste(savedir,"heatmap",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"heatmap/CD8_ADT.pdf",sep =""))
DoHeatmap(CD8_ADT, features = rownames(CD8_ADT))
dev.off()

DefaultAssay(CD8_ADT) <- "ADT"
dir.create(paste(savedir,"heatmap",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"heatmap/CD8_ADT_Run.pdf",sep =""))
DoHeatmap(CD8_ADT, features = rownames(CD8_ADT),group.by = "Run")
dev.off()

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")
objname = "CD8"
process = "integrated"
Assay = "ADT"

CD8_ADT_integrated <- ADT_merging(CD8_ADT, savedir, dims = 6, numfeatures=36,
                                  Assay=Assay, process=process, objname=objname,
                                  sample_tree = NULL, split_by = "Run",
                                  reference=NULL)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "CD8"
process = "integration"
Assay = "ADT"
CD8_ADT_integrated <- RNA_integration(obj_path = CD8_ADT_integrated, saveDir = savedir, dims = 10,
                                      RNA_features = c("CD4-protein","CD8a-protein"),Assay = Assay,
                                      process = process, objname = objname, ncol = 4)

### Nearest Neighbour and Clustering
CD8_ADT_integrated_NN <- FindNeighbors(CD8_ADT_integrated, dims = 1:10)

pdf(paste(savedir,"UMAP/CD8_ADT_Age.pdf",sep = ""))
DimPlot(CD8_ADT_integrated_NN, group.by = "Age")
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_batch.pdf",sep = ""))
DimPlot(CD8_ADT_integrated_NN, group.by = "Run")
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_Age_split.pdf",sep = ""), width = 10, height = 6)
DimPlot(CD8_ADT_integrated_NN, group.by = "Age", split.by = "Age", ncol = 2)
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_batch_split.pdf",sep = ""), width = 15, height = 6)
DimPlot(CD8_ADT_integrated_NN, group.by = "Run", split.by = "Run", ncol = 3)
dev.off()

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8)
res = 0.6
for (i in 1:length(res)) {
  process = paste("ADT_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  CD8_ADT_integrated_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD8_ADT_integrated_NN, dims = 10, res = res[i], saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_ADT","nFeature_ADT"), objname = "CD8_ADT_integrated_NN_cluster",
                                                       process = process,col_sel = c("Age","Run","orig.ident","gender"))
}

### Generating the MAGIC ADT and MAGIC RNA
saveRDS(CD8_ADT_integrated_NN_cluster,paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster.RDS",sep = ""))

## Since Magic does not work on the conda seurat_new we have to perform the analysis outside conda.
## ml load r/4.2.2
library(Seurat)
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Bulk/"
CD8_ADT_integrated_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster.RDS",sep = ""))
DefaultAssay(CD8_ADT_integrated_NN_cluster) <- "ADT"
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/miniconda3/envs/seurat_new/bin/python3.8")
library(Rmagic)
py_discover_config("magic") # to check
CD8_ADT_integrated_NN_cluster_impute <- magic(CD8_ADT_integrated_NN_cluster, npca=10) ## imputing the RNA data as for RNA PCs are 20
saveRDS(CD8_ADT_integrated_NN_cluster_impute,paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster_impute.RDS",sep = ""))

# CD8_integrated_RNA_NN_cluster <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run06/Downstream/saveRDS_obj/CD8_integrated_RNA_NN_cluster.RDS")
#### CD8 Modality Integration ######
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITEseq_modality_integration.R")
# Find the last save object for RNA and ADT
# we have RNA already in the path
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Bulk/regress_TCR/"
# CD8_ADT_integrated_NN_cluster <- readRDS(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Bulk/saveRDS_obj/CD8_ADT_integrated_NN_cluster_impute.RDS",sep = ""))
# CD8_integrated_RNA_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster_impute.RDS",sep = ""))
RNA_obj_path <- CD8_integrated_RNA_NN_cluster
ADT_obj_path <- CD8_ADT_integrated_NN_cluster

objname <- "CD8"
process <- "modality_integrate"
Assay <- "integrated"

CD8_modality_integrate <- modality_integration(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path,
                                               RNA_dims = 20, ADT_dims = 10,
                                               saveDir = savedir, Assay = Assay, process = process, objname = objname)

DefaultAssay(CD8_modality_integrate) <- "integrated"

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/Modal_cluster_and_UMAP.R")
obj_path <- CD8_modality_integrate
RNA_features = c("CD4","CD8A")
ADT_features = c("CD4-protein","CD8a-protein")

# readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run03/analysis/downstream/saveRDS_obj/CD8_modality_integrate_cluster_1.RDS")
res = c(0.2,0.4,0.6,0.8,1)
res=0.6
for (i in 1:length(res)) {
  objname = "CD8"
  Assay = "integrated"
  process = paste("modality_cluster_and_QC",res[i],sep = "_")
  CD8_modality_integrate_cluster <- Modal_Cluster_and_UMAP_QC(modal_integrate_obj_path = obj_path, saveDir = savedir,
                                                              res = res[i], RNA_features = RNA_features,
                                                              ADT_features = ADT_features,Assay = Assay, protein_assay = "ADT",
                                                              process = process, objname=objname,
                                                              col_sel=c("Age","Run","orig.ident","gender"))
}

DefaultAssay(CD8_modality_integrate_cluster) <- "RNA"
clus_7 <- FindMarkers(CD8_modality_integrate_cluster, ident.1 =c(7), only.pos = TRUE)
clus_5 <- FindMarkers(CD8_modality_integrate_cluster, ident.1 =c(5), only.pos = TRUE)
clus_11 <- FindMarkers(CD8_modality_integrate_cluster, ident.1 =c(11), only.pos = TRUE)
clus_14 <- FindMarkers(CD8_modality_integrate_cluster, ident.1 =c(14), only.pos = TRUE)
clus_15 <- FindMarkers(CD8_modality_integrate_cluster, ident.1 =c(15), only.pos = TRUE)

saveRDS(CD8_modality_integrate_cluster,paste(savedir,"saveRDS_obj/CD8_modality_integrate_cluster_1.RDS",sep = ""))

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/ADT_featureplot.R")
ADT_featureplot(obj = CD8_modality_integrate_cluster, savedir = savedir,objname = "CD8_modality_integrate_cluster",
                process = "featureplots",reduction = "wnn.umap",x = "wnnUMAP_1",y = "wnnUMAP_2")

CD8_samplewise_counts <- table(CD8_modality_integrate_cluster@meta.data$orig.ident) %>% as.data.frame()
write.table(CD8_samplewise_counts, paste(savedir,"Table/CD8_samplewise_counts.txt",sep = ""), sep = "\t",
            row.names = F, col.names = T, quote = F)

Ag_pos_cells <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Table/Antigen_positive_cells.txt", header = TRUE)
Ag_pos_cells$Ag <- gsub(".*.-1_","",Ag_pos_cells$barcode_Ag)

CD8_modality_integrate_cluster@meta.data$Ag_gt_50 <- Ag_pos_cells[match(rownames(CD8_modality_integrate_cluster@meta.data),Ag_pos_cells$sample_id_barcode),"Ag"]

CD8_modality_integrate_cluster@meta.data$Ag_Age <- paste(CD8_modality_integrate_cluster@meta.data$Ag_gt_50, CD8_modality_integrate_cluster@meta.data$Age,sep = "_")

Ag_Age <- unique(CD8_modality_integrate_cluster@meta.data$Ag_Age)
Ag_Age2 <- Ag_Age[grep("NA_O|NA_Y",Ag_Age,invert=TRUE)]

rm(plot_list)
plot_list <- list()
for (i in 1:length(Ag_Age2)) {
  cellnames <- CD8_modality_integrate_cluster@meta.data[grep(Ag_Age2[i], CD8_modality_integrate_cluster@meta.data$Ag_Age),] %>% rownames()
  p <- DimPlot(CD8_modality_integrate_cluster,
               cells.highlight = cellnames, 
               reduction = "wnn.umap", 
               label = FALSE, cols.highlight = "deeppink2",
               sizes.highlight = 0.5,
               cols = "gray92") + 
    ggtitle(paste(Ag_Age2[i])) + 
    NoLegend() + 
    theme(plot.title = element_text(hjust = 0.2))
  # pdf(paste(savedir,"CD4/CD4_TRB_",TRB_clonal_noNA_gt_2$Var1[i],".pdf",sep = ""))
  # print(p)
  # dev.off()
  plot_list[[i]] <- p
}

pdf(paste(savedir,"UMAP/Ag_gt_50_Age.pdf",sep = ""), width = 5.5, height = 5)
plot_list
dev.off()

CD8_modality_integrate_cluster@meta.data$Ag_group <- NA
CD8_modality_integrate_cluster@meta.data[grep("EBV_B",CD8_modality_integrate_cluster@meta.data$Ag_gt_50),"Ag_group"] <- "Lytic"
CD8_modality_integrate_cluster@meta.data[grep("EBV_EB|EBV_LM",CD8_modality_integrate_cluster@meta.data$Ag_gt_50),"Ag_group"] <- "Latent"
CD8_modality_integrate_cluster@meta.data$Ag_group_Age <- paste(CD8_modality_integrate_cluster@meta.data$Ag_group,CD8_modality_integrate_cluster$Age,sep="_")
group_Ag <- grep("NA_",unique(CD8_modality_integrate_cluster@meta.data$Ag_group_Age),invert = TRUE,value = TRUE)

rm(plot_list)
plot_list <- list()
for (i in 1:length(group_Ag)) {
  cellnames <- CD8_modality_integrate_cluster@meta.data[grep(group_Ag[i], CD8_modality_integrate_cluster@meta.data$Ag_group_Age),] %>% rownames()
  p <- DimPlot(CD8_modality_integrate_cluster,
               cells.highlight = cellnames,
               reduction = "wnn.umap",
               label = FALSE, cols.highlight = "deeppink2",
               sizes.highlight = 0.5,
               cols = "gray92") + 
    ggtitle(paste(group_Ag[i])) + 
    NoLegend() + 
    theme(plot.title = element_text(hjust = 0.2))
  # pdf(paste(savedir,"CD4/CD4_TRB_",TRB_clonal_noNA_gt_2$Var1[i],".pdf",sep = ""))
  # print(p)
  # dev.off()
  plot_list[[i]] <- p
}

pdf(paste(savedir,"UMAP/Ag_group_gt_50_Age.pdf",sep = ""), width = 5.5, height = 5)
plot_list
dev.off()

pdf(paste(savedir,"UMAP/Ag_group_gt_50_Age.pdf",sep = ""), width = 11, height = 5)
DimPlot(CD8_modality_integrate_cluster, split.by = "Age", group.by = "Age", ncol = 2, reduction = "wnn.umap")
dev.off()


# CD4_genes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run01/Analysis/downstream/Table/CD4_genes")[,1]
# dir.create(paste(savedir,"dotplot",sep = ""),showWarnings = FALSE)
# pdf(paste(savedir,"dotplot/CD4_genes_CD4_dotplot.pdf",sep = ""))
# DotPlot(CD4_modality_integrate_cluster,"RNA",CD4_genes) + coord_flip()
# dev.off()
# 
# CD8_genes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run01/Analysis/downstream/Table/CD8_genes")[,1]
# dir.create(paste(savedir,"dotplot",sep = ""),showWarnings = FALSE)
# pdf(paste(savedir,"dotplot/CD8_genes_CD4_dotplot.pdf",sep = ""), width = 7, height = 9)
# DotPlot(CD4_modality_integrate_cluster,"RNA",CD8_genes) + coord_flip()
# dev.off()
# 
# ### exhauster genes
# ex_genes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run01/Analysis/downstream/Table/exhausted_genes_2")
# ex_genes_up <- toupper(ex_genes[,1])
# pdf(paste(savedir,"dotplot/exhausted_CD4_dotplot.pdf",sep = ""), height = 10, width = 8)
# DotPlot(CD4_modality_integrate_cluster,"RNA",ex_genes_up) + coord_flip()
# dev.off()
# 
# DefaultAssay(CD4_modality_integrate_cluster) <- "MAGIC_RNA"
# genes <- c("IKZF2","FOXP3","TCF7","LEF1","IRF2","RUNX3","EOMES","TBX21","PRDM1","XBP1",
#            "KLF6","NOTCH1","SLAMF6","FOSL2","FOS","RORC","JUN","JUND","JUNB","STAT3",
#            "Zeb2","ID2","Tox","Myb","TCF7","ID3","BCL6","myc","IL7R","CD84","IL2RB")
# genes <- toupper(genes)
# genes <- rownames(CD4_modality_integrate_cluster)[match(toupper(genes),rownames(CD4_modality_integrate_cluster))]

genes <- c("MTRNR2L8","TRBV20-1", "TRAV1-2", "TRBV1", "TRAV7-9", "TRAV22", "TRAV19","TRBV7-9", "TRBV20-1", "TRBV7-9")
DefaultAssay(CD8_modality_integrate_cluster) <- "MAGIC_RNA"
genes <- grep(paste(genes, collapse = "|"),rownames(CD8_modality_integrate_cluster@assays$RNA@counts),value = TRUE)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p <- featureplot_front(CD8_modality_integrate_cluster, genes[i],
                         reduction = "wnn.umap", x="wnnUMAP_1", y="wnnUMAP_2",size=0.2) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- p
}

dir.create(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Bulk/featureplot/",sep =""),showWarnings = FALSE)
pdf(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Bulk/featureplot/CD4_Goronzy.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list)
dev.off()

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/Gene_score.R")
CD4_modality_integrate_cluster <- genescore_abhinav(obj = CD4_modality_integrate_cluster,assay = "RNA",
                                                    savedir = savedir,
                                                    cores = 5,objname = "CD4")

### Bulk Pseudobulk ####
NA_SCM <- c(0,4,14,15,17,20)
CM <- c(5,7,11)
EM1 <- c(1,10,12)
EM2 <- c(8,16,18,19)
EM3 <- c(6,13)
MAIT <- 9
TEMRA <- c(2,3)

memory = c(1,2,3,5,6,7,8,9,10,11,12,13,16,18,19)
celltypes <- c("memory")

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Bulk/regress_TCR/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)

for (i in 1:length(celltypes)) {
  try({
    # cluster = get(celltypes[i])
    cluster = "all"
    group1 = "O"
    group2 = "Y"
    remove_samples = NULL
    remove_samples_name <- paste(remove_samples,collapse="_and_")
    DefaultAssay(CD8_EBV) <- "RNA"
    # savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run01_to_10_analysis/vaccine/S_vs_Z/"
    source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/COVID_pseudobulk_PCA_within_AJ_ADT.R")
    CD8_subset_dds <- COVID_pseudobulk_within_cluster_AJ(obj = CD8_EBV, savedir = savedir, group1 = group1, group2 = group2,
                                                         grouping_by = "Age", cluster = cluster, cell_freq = 20,remove_samples = remove_samples,
                                                         cluster_group = "seurat_clusters",sample_col = "orig.ident", batch_col = "Run",
                                                         gene_min_counts =  5, 
                                                         column_name_split = c("sampleid","virus","age_number","gender","Run","Age"))
    
    cluster2 <- paste(cluster, sep="_", collapse="_")
    savedir2 <- paste(savedir,"pseudobulk/clus_",cluster2,"_removed_",remove_samples_name,"_",group1,"_vs_",group2,"/",sep = "")
    # savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run04/Analysis/pseudobulk/clus_13_removed__DMSO_vs_Sprotein/"
    source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/RNAseq/pipeline_function/RNAseq_limma_EdgeR.R")
    design0 <- model.matrix(~ 0 + Age + Run, data = colData(CD8_subset_dds))
    colnames(design0) <- c(group1,group2,paste("Run",1:(ncol(design0)-2),sep = ""))
    cm <- makeContrasts(O_VS_Y = O-Y,levels = design0)
    desl_clus <- LimmaEdgeR_differential(dds = CD8_subset_dds,
                                         design0 = design0,
                                         cm = cm, 
                                         savedir = savedir2,
                                         logfc = 0.5,
                                         p_value_adj = 0.05)
  })
}


#### ADT
NA_SCM <- c(0,4,14,15,17,20)
CM <- c(5,7,11)
EM1 <- c(1,10,12)
EM2 <- c(8,16,18,19)
EM3 <- c(6,13)
MAIT <- 9
TEMRA <- c(2,3)

memory = c(1,2,3,5,6,7,8,9,10,11,12,13,16,18,19)
celltypes <- c("memory","NA_SCM","CM","EM1","EM2","EM3","MAIT","TEMRA")

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Bulk/regress_TCR/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)

for (i in 1:length(celltypes)) {
  try({
    cluster = "all"
    group1 = "O"
    group2 = "Y"
    remove_samples = NULL
    remove_samples_name <- paste(remove_samples,collapse="_and_")
    DefaultAssay(CD8_EBV) <- "RNA"
    # savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run01_to_10_analysis/vaccine/S_vs_Z/"
    source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/COVID_pseudobulk_PCA_within_AJ_2.R")
    CD8_subset_dds <- COVID_pseudobulk_within_cluster_AJ(obj = CD8_EBV, savedir = savedir, group1 = group1, group2 = group2,
                                                         grouping_by = "Age", cluster = cluster,cell_freq = 20,remove_samples = remove_samples,
                                                         cluster_group = "seurat_clusters",sample_col = "orig.ident", batch_col = "Run",
                                                         gene_min_counts =  5, 
                                                         column_name_split = c("sampleid","virus","age_number","gender","Run","Age"))
    
    cluster2 <- paste(cluster, sep="_", collapse="_")
    savedir2 <- paste(savedir,"pseudobulk/clus_",cluster2,"_removed_",remove_samples_name,"_",group1,"_vs_",group2,"/",sep = "")
    # savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run04/Analysis/pseudobulk/clus_13_removed__DMSO_vs_Sprotein/"
    source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/RNAseq/pipeline_function/RNAseq_limma_EdgeR.R")
    design0 <- model.matrix(~ 0 + Age + Run, data = colData(CD8_subset_dds))
    colnames(design0) <- c(group1,group2,paste("Run",1:(ncol(design0)-2),sep = ""))
    cm <- makeContrasts(O_VS_Y = O-Y,levels = design0)
    desl_clus <- LimmaEdgeR_differential(dds = CD8_subset_dds,
                                         design0 = design0,
                                         cm = cm, 
                                         savedir = savedir2,
                                         logfc = 0.5,
                                         p_value_adj = 0.05)
  })
}



### Antigen ####
library("Matrix")
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.

maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/"
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/regress_TCR/"
dir.create(savedir,showWarnings = FALSE)
Ag_table <- read.table(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/Table/Exp01_02_03_Ag_table.txt",sep = ""), sep = "\t", header = T)

# combined_ADT <- readRDS(paste(maindir, "saveRDS_obj/combined_ADT.RDS",sep = ""))
# combined <- readRDS(paste(maindir, "saveRDS_obj/combined.RDS",sep = ""))

### Since the TCR genes and the MTRNR2L8 are making various clusters that is not required so we will remove these genes from the object
remove_genes <- c("MTRNR2L8",grep("^TRAV|^TRAJ|^TRBC|^TRBV|^TRBJ|^TRBC",
                                  rownames(CD8_modality_integrate_cluster@assays$RNA@counts),
                                  value=TRUE))
keep_genes <- grep(paste(remove_genes,collapse = "|"),
                   rownames(CD8_modality_integrate_cluster@assays$RNA@counts), 
                   invert = TRUE, value = TRUE)
combined_no_TCR <- combined[keep_genes,]

### Just to check
grep("MTRNR2L8|^TRAV|^TRAJ|^TRBC|^TRBV|^TRBJ|^TRBC",
                  rownames(combined_no_TCR@assays$RNA@counts),
                  value=TRUE)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_subsetting.R")
RNA_obj_path = combined_no_TCR
ADT_obj_path = combined_ADT
ADT_main = combined_ADT
CD8_cluster <- NULL
Ag_cellnames <- paste(Ag_table$samplename,Ag_table$barcode,sep="_")
objname = "CD8"
Assay = "RNA"
process = "Ag_subset"

CD8_RNA <- RNA_subsetting(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path, ADT_main = ADT_main, subset_cluster = CD8_cluster,
                          cellnames=Ag_cellnames, ngenes = 4000, saveDir = savedir, objname = objname, Assay = Assay, process = process)

CD8_RNA <- RunUMAP(CD8_RNA, dims = 1:20, reduction = "pca")

dir.create(paste(savedir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/CD8_no_integration_samples.pdf",sep = ""),width = 10, height = 5)
print(DimPlot(CD8_RNA, group.by = "orig.ident"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_no_integration_samples_splitted.pdf",sep = ""),width = 15, height = 15)
print(DimPlot(CD8_RNA, split.by = "orig.ident", group.by="orig.ident", ncol = 4) + NoLegend())
dev.off()

dir.create(paste(savedir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/CD8_no_integration_Run.pdf",sep = ""),width = 6, height = 5)
print(DimPlot(CD8_RNA, group.by = "Run"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_no_integration_Run_splitted.pdf",sep = ""),width = 15, height = 5)
print(DimPlot(CD8_RNA, split.by = "Run", group.by="Run", ncol = 3) + NoLegend())
dev.off()

## Since there is a batch effect we have integrate based on the Run
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_sctransform_V2.R")
objname = "CD8"
Assay = "RNA"
process = "sctransform"
CD8_RNA_sctransformed <- sctransform_V2_integration(obj = CD8_RNA, saveDir = savedir, ngenes = 4000,
                                                    regress = c("nCount_RNA"),
                                                    dims = 20,
                                                    Assay = Assay, 
                                                    process = process, 
                                                    objname = objname,
                                                    split_by = "Run",
                                                    reference = NULL,
                                                    sample_tree = NULL)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "CD8"
process = "integration"
Assay = "RNA"
CD8_integrated_RNA <- RNA_integration(CD8_RNA_sctransformed, savedir, dims = 20, RNA_features = c("CD4","CD8A"),
                                      Assay=Assay, process=process, objname=objname, ncol = 4, ndims = 50)

CD8_integrated_RNA_NN <- FindNeighbors(CD8_integrated_RNA, dims = 1:20)

CD8_integrated_RNA_NN@meta.data$samplename <- gsub("_[A|T|G|C].*","",rownames(CD8_integrated_RNA_NN@meta.data))
CD8_integrated_RNA_NN@meta.data$sample_id = gsub("_.*","",CD8_integrated_RNA_NN@meta.data$orig.ident)
CD8_integrated_RNA_NN@meta.data$Run = gsub(".*_","",CD8_integrated_RNA_NN@meta.data$orig.ident)
CD8_integrated_RNA_NN@meta.data$Age <- gsub("[0-9]","",CD8_integrated_RNA_NN@meta.data$sample_id)
CD8_integrated_RNA_NN@meta.data$gender <- gsub(".*._M_.*.","M",CD8_integrated_RNA_NN@meta.data$orig.ident) %>% gsub(".*._F_.*.","F",.)
CD8_integrated_RNA_NN@meta.data$Age_number <- gsub(".*._EBV_","",CD8_integrated_RNA_NN@meta.data$orig.ident) %>% gsub("_.*.","",.)


pdf(paste(savedir,"UMAP/CD8_RNA_Age.pdf",sep = ""))
DimPlot(CD8_integrated_RNA_NN, group.by = "Age")
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_batch.pdf",sep = ""))
DimPlot(CD8_integrated_RNA_NN, group.by = "Run")
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_Age_split.pdf",sep = ""), width = 10, height = 6)
DimPlot(CD8_integrated_RNA_NN, group.by = "Age", split.by = "Age", ncol = 2)
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_batch_split.pdf",sep = ""), width = 15, height = 6)
DimPlot(CD8_integrated_RNA_NN, group.by = "Run", split.by = "Run", ncol = 3)
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_nCount_RNA_SCT.pdf",sep = ""), width = 10, height = 6)
FeaturePlot(CD8_integrated_RNA_NN, c("nCount_SCT", "nCount_RNA"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_nCount_RNA_SCT_Vln.pdf",sep = ""), width = 10, height = 6)
VlnPlot(CD8_integrated_RNA_NN, c("nCount_SCT", "nCount_RNA"), group.by = "orig.ident")
dev.off()

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8,1)
res = 0.4
for (i in 1:length(res)) {
  process = paste("RNA_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  CD8_integrated_RNA_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD8_integrated_RNA_NN, dims = 30, res = res[i], 
                                                       saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_RNA","nFeature_RNA"), 
                                                       objname = "CD8_integrated_RNA_NN_cluster",
                                                       process = process, col_sel = c("Age","Run","orig.ident","gender"))
}

saveRDS(CD8_integrated_RNA_NN_cluster,paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster.RDS",sep = ""))

## Since Magic does not work on the conda seurat_new we have to perform the analysis outside conda.
## ml load r/4.2.2
library(Seurat)
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/regress_TCR/"
CD8_integrated_RNA_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster.RDS",sep = ""))
DefaultAssay(CD8_integrated_RNA_NN_cluster) <- "RNA"
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/miniconda3/envs/seurat_new/bin/python3.8")
library(Rmagic)
py_discover_config("magic") # to check
CD8_integrated_RNA_NN_cluster_impute <- magic(CD8_integrated_RNA_NN_cluster, npca=30) ## imputing the RNA data as for RNA PCs are 20
saveRDS(CD8_integrated_RNA_NN_cluster_impute,paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster_impute.RDS",sep = ""))

# DefaultAssay(CD8_integrated_RNA_NN_cluster) <- "RNA"
# Y3_Y6_cluster <- FindMarkers(CD8_integrated_RNA_NN_cluster, ident.1 = c(7,10,25,26,27), only.pos = FALSE)
# Y3_Y6_cluster_between_naive <- FindMarkers(CD8_integrated_RNA_NN_cluster, ident.1 = c(10,25,26,27), ident.2 = c(1,20,22,23), only.pos = FALSE)
# Y3_Y6_no_clus7_cluster_between_naive <- FindMarkers(CD8_integrated_RNA_NN_cluster, ident.1 = c(10,25,26,27), ident.2 = c(1,20,22,23), only.pos = FALSE)
# Y3_Y6_cluster_no_clus7 <- FindMarkers(CD8_integrated_RNA_NN_cluster, ident.1 = c(10,25,26,27), only.pos = FALSE)
# Y3_Y6_cluster_no_clus7_ordered <- Y3_Y6_cluster_no_clus7[order(-Y3_Y6_cluster_no_clus7$avg_log2FC),]
# 
# genes <- rownames(Y3_Y6_cluster_no_clus7_ordered)[c(1:10,188:196)]
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
# rm(plot_list)
# plot_list <- list()
# for (i in 1:length(genes)) {
#   p <- featureplot_front(CD8_integrated_RNA_NN_cluster, genes[i],
#                          reduction = "umap", x="UMAP_1", y="UMAP_2",size=0.2) +
#     scale_color_gradientn(colours = ArchRPalettes$solarExtra)
#   plot_list[[i]] <- p
# }
# 
# dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
# pdf(paste(savedir,"featureplot/Y3_Y6_high_expression_no_clus_7.pdf",sep = ""))
# print(plot_list)
# dev.off()
# 
# library(scales)
# CD8 <- CD8_integrated_RNA_NN_cluster
# dir.create(paste(savedir,"UMAP/clusters/",sep = ""), showWarnings = FALSE)
# clusters = 0:27
# CD8@meta.data$seurat_clusters <- factor(CD8@meta.data$seurat_clusters, levels = clusters)
# rm(plot_list)
# plot_list <- list()
# for (i in 1:length(clusters)) {
#   col <- rep("white",length(clusters))
#   col[i] <- "red"
#     p <- DimPlot(CD8, group.by = "seurat_clusters", label = TRUE, cols=alpha(col,0.3), reduction = "umap")
#     pdf(paste(savedir,"UMAP/clusters/CD8_",clusters[i],"_UMAP.pdf",sep = ""))
#     print(p)
#     dev.off()
#     plot_list[[i]] <- p
# }
# 
# dir.create(paste(savedir,"UMAP/",sep =""),showWarnings = FALSE)
# pdf(paste(savedir,"UMAP/clusters/CD8_all_combine.pdf",sep = ""), width = 4.5, height = 4)
# plot_list
# dev.off()

### CD8 ADT #######
CD8_ADT = readRDS(paste(savedir,"saveRDS_obj/CD8_Ag_subset_ADT.RDS",sep = ""))

CD8_ADT@meta.data$samplename <- gsub("_[A|T|G|C].*","",rownames(CD8_ADT@meta.data))
CD8_ADT@meta.data$sample_id = gsub("_.*","",CD8_ADT@meta.data$orig.ident)
CD8_ADT@meta.data$Run = gsub(".*_","",CD8_ADT@meta.data$orig.ident)
CD8_ADT@meta.data$Age <- gsub("[0-9]","",CD8_ADT@meta.data$sample_id)
CD8_ADT@meta.data$gender <- gsub(".*._M_.*.","M",CD8_ADT@meta.data$orig.ident) %>% gsub(".*._F_.*.","F",.)
CD8_ADT@meta.data$Age_number <- gsub(".*._EBV_","",CD8_ADT@meta.data$orig.ident) %>% gsub("_.*.","",.)

### No integration
CD8_ADT <- NormalizeData(CD8_ADT, normalization.method = 'CLR', margin = 2)
CD8_ADT <- FindVariableFeatures(CD8_ADT, selection.method = "vst", nfeatures = 30)
CD8_ADT <- ScaleData(CD8_ADT, verbose = FALSE)
CD8_ADT <- RunPCA(CD8_ADT, npcs = 30, approx=FALSE)
CD8_ADT <- RunUMAP(CD8_ADT, dims = 1:10)

dir.create(paste(savedir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/CD8_ADT_no_integration_samples.pdf",sep = ""),width = 10, height = 5)
print(DimPlot(CD8_ADT, group.by = "orig.ident"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_no_integration_samples_splitted.pdf",sep = ""),width = 15, height = 15)
print(DimPlot(CD8_ADT, split.by = "orig.ident", group.by="orig.ident", ncol = 4) + NoLegend())
dev.off()

dir.create(paste(savedir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/CD8_ADT_no_integration_Run.pdf",sep = ""),width = 6, height = 5)
print(DimPlot(CD8_ADT, group.by = "Run"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_no_integration_Run_splitted.pdf",sep = ""),width = 15, height = 5)
print(DimPlot(CD8_ADT, split.by = "Run", group.by="Run", ncol = 3) + NoLegend())
dev.off()

DefaultAssay(CD8_ADT) <- "ADT"
dir.create(paste(savedir,"heatmap",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"heatmap/CD8_ADT.pdf",sep =""))
DoHeatmap(CD8_ADT, features = rownames(CD8_ADT))
dev.off()

DefaultAssay(CD8_ADT) <- "ADT"
dir.create(paste(savedir,"heatmap",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"heatmap/CD8_ADT_Run.pdf",sep =""))
DoHeatmap(CD8_ADT, features = rownames(CD8_ADT),group.by = "Run")
dev.off()

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")
objname = "CD8"
process = "integrated"
Assay = "ADT"

CD8_ADT_integrated <- ADT_merging(CD8_ADT, savedir, dims = 6, numfeatures=36,
                                  Assay=Assay, process=process, objname=objname,
                                  sample_tree = NULL, split_by = "Run",
                                  reference=NULL)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "CD8"
process = "integration"
Assay = "ADT"
CD8_ADT_integrated <- RNA_integration(obj_path = CD8_ADT_integrated, saveDir = savedir, dims = 10,
                                      RNA_features = c("CD8-protein","CD8a-protein"),Assay = Assay,
                                      process = process, objname = objname, ncol = 4)

### Nearest Neighbour and Clustering
CD8_ADT_integrated_NN <- FindNeighbors(CD8_ADT_integrated, dims = 1:10)

pdf(paste(savedir,"UMAP/CD8_ADT_Age.pdf",sep = ""))
DimPlot(CD8_ADT_integrated_NN, group.by = "Age")
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_batch.pdf",sep = ""))
DimPlot(CD8_ADT_integrated_NN, group.by = "Run")
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_Age_split.pdf",sep = ""), width = 10, height = 6)
DimPlot(CD8_ADT_integrated_NN, group.by = "Age", split.by = "Age", ncol = 2)
dev.off()

pdf(paste(savedir,"UMAP/CD8_ADT_batch_split.pdf",sep = ""), width = 15, height = 6)
DimPlot(CD8_ADT_integrated_NN, group.by = "Run", split.by = "Run", ncol = 3)
dev.off()

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8)
res = 0.6
for (i in 1:length(res)) {
  process = paste("ADT_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  CD8_ADT_integrated_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD8_ADT_integrated_NN, dims = 10, res = res[i], saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_ADT","nFeature_ADT"), objname = "CD8_ADT_integrated_NN_cluster",
                                                       process = process,col_sel = c("Age","Run","orig.ident","gender"))
}

### Generating the MAGIC ADT and MAGIC RNA
saveRDS(CD8_ADT_integrated_NN_cluster,paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster.RDS",sep = ""))

## Since Magic does not work on the conda seurat_new we have to perform the analysis outside conda.
## ml load r/4.2.2
library(Seurat)
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/"
CD8_ADT_integrated_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster.RDS",sep = ""))
DefaultAssay(CD8_ADT_integrated_NN_cluster) <- "ADT"
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/miniconda3/envs/seurat_new/bin/python3.8")
library(Rmagic)
py_discover_config("magic") # to check
CD8_ADT_integrated_NN_cluster_impute <- magic(CD8_ADT_integrated_NN_cluster, npca=10) ## imputing the RNA data as for RNA PCs are 20
saveRDS(CD8_ADT_integrated_NN_cluster_impute,paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster_impute.RDS",sep = ""))

# CD8_integrated_RNA_NN_cluster <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run06/Downstream/saveRDS_obj/CD8_integrated_RNA_NN_cluster.RDS")
#### CD8 Modality Integration ######
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITEseq_modality_integration.R")
# Find the last save object for RNA and ADT
# we have RNA already in the path
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/regress_TCR/"
## Since ADT we donot have to regress the TCR we can take it from the previous analysis
# CD8_ADT_integrated_NN_cluster <- readRDS(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/saveRDS_obj/CD8_ADT_integrated_NN_cluster_impute.RDS",sep = ""))
# CD8_integrated_RNA_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster_impute.RDS",sep = ""))
RNA_obj_path <- CD8_integrated_RNA_NN_cluster
ADT_obj_path <- CD8_ADT_integrated_NN_cluster

objname <- "CD8"
process <- "modality_integrate"
Assay <- "integrated"
CD8_modality_integrate <- modality_integration(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path,
                                               RNA_dims = 20, ADT_dims = 10,
                                               saveDir = savedir, Assay = Assay, process = process, objname = objname)

DefaultAssay(CD8_modality_integrate) <- "integrated"
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/Modal_cluster_and_UMAP.R")
obj_path <- CD8_modality_integrate
RNA_features = c("CD4","CD8A")
ADT_features = c("CD4-protein","CD8a-protein")

# readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run03/analysis/downstream/saveRDS_obj/CD8_modality_integrate_cluster_1.RDS")
res = c(1.2,1.4,1.6,1.8,2)
res=1
for (i in 1:length(res)) {
  objname = "CD8"
  Assay = "integrated"
  process = paste("modality_cluster_and_QC",res[i],sep = "_")
  CD8_modality_integrate_cluster <- Modal_Cluster_and_UMAP_QC(modal_integrate_obj_path = obj_path, saveDir = savedir,
                                                              res = res[i], RNA_features = RNA_features,
                                                              ADT_features = ADT_features,Assay = Assay, protein_assay = "ADT",
                                                              process = process, objname=objname,
                                                              col_sel=c("Age","Run","orig.ident","gender"))
}

genes <- c("MTRNR2L8","TRBV20-1", "TRAV1-2", "TRBV1", "TRAV7-9", "TRAV22", "TRAV19","TRBV7-9", "TRBV20-1", "TRBV7-9")
DefaultAssay(CD8_modality_integrate_cluster) <- "MAGIC_RNA"
genes <- grep(paste(genes, collapse = "|"),rownames(CD8_modality_integrate_cluster@assays$RNA@counts),value = TRUE)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
rm(plot_list)
plot_list <- list()
for (i in 1:length(genes)) {
  p <- featureplot_front(CD8_modality_integrate_cluster, genes[i],
                         reduction = "wnn.umap", x="wnnUMAP_1", y="wnnUMAP_2",size=0.2) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra)
  plot_list[[i]] <- p
}

dir.create(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/featureplot/",sep =""),showWarnings = FALSE)
pdf(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/featureplot/CD4_Goronzy.pdf",sep = ""), width = 4.5, height = 4)
print(plot_list)
dev.off()


DefaultAssay(CD8_modality_integrate_cluster) <- "RNA"
clus6_17 <- FindMarkers(CD8_modality_integrate_cluster, ident.1 =c(6,17), only.pos = TRUE)

clus_9 <- FindMarkers(CD8_modality_integrate_cluster, ident.1 =c(9), only.pos = TRUE)
clus_11 <- FindMarkers(CD8_modality_integrate_cluster, ident.1 =c(11), only.pos = TRUE)
clus_12 <- FindMarkers(CD8_modality_integrate_cluster, ident.1 =c(12), only.pos = TRUE)
clus_13 <- FindMarkers(CD8_modality_integrate_cluster, ident.1 =c(13), only.pos = TRUE)
clus_14 <- FindMarkers(CD8_modality_integrate_cluster, ident.1 =c(14), only.pos = TRUE)
clus_6 <- FindMarkers(CD8_modality_integrate_cluster, ident.1 =c(6), only.pos = TRUE)

CD8_EBV@meta.data$Age <- gsub("[0-9]_.*.","",CD8_EBV@meta.data$orig.ident)
CD8_EBV@meta.data$orig.ident <- paste(CD8_EBV@meta.data$samplename,CD8_EBV@meta.data$Age,sep = "_")

cluster <- c(0,4,14,15,17,20,5,7,11,1,10,12,8,16,18,19,6,13,9,2,3)
celltypes <- c(rep("naive_SCM",6),rep("CM",3),rep("EM1",3),rep("EM2",4),"EM3","EM3","MAIT","TEMRA","TEMRA")

library(stringr)
patterns = cluster
replacements = celltypes
# names(replacement) <- patterns
CD8_EBV@meta.data$celltypes <- CD8_EBV@meta.data$seurat_clusters
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  CD8_EBV@meta.data$celltypes <- str_replace_all(CD8_EBV@meta.data$celltypes, pattern, replacements[i])
}

pdf("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Abhinav/sc_aging/vlnplot/CD8_EBV/CD8_EBV_celltypes.pdf")
DimPlot(CD8_EBV, group.by = "celltypes", reduction = "wnn.umap", label = TRUE, label.size = 5)
dev.off()

pdf("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Abhinav/sc_aging/vlnplot/CD8_EBV/CD8_EBV_dimplot.pdf")
DimPlot(CD8_EBV, reduction = "wnn.umap", label = TRUE, label.size = 5)
dev.off()

saveRDS(CD8_EBV,paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Bulk/regress_TCR/saveRDS_obj/CD8_modality_integrate_cluster.RDS",sep = ""))

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/ADT_featureplot.R")
ADT_featureplot(obj = CD8_modality_integrate_cluster, savedir = savedir,objname = "CD8_modality_integrate_cluster",
                process = "featureplots",reduction = "wnn.umap",x = "wnnUMAP_1",y = "wnnUMAP_2")

CD8_samplewise_counts <- table(CD8_modality_integrate_cluster@meta.data$orig.ident) %>% as.data.frame()
write.table(CD8_samplewise_counts, paste(savedir,"Table/CD8_samplewise_counts.txt",sep = ""), sep = "\t",
            row.names = F, col.names = T, quote = F)

Ag_pos_cells <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/Table/Exp01_02_03_Ag_table.txt", 
                           header = TRUE, sep = "\t")

Ag_pos_cells$sample_id_barcode <- paste(Ag_pos_cells$samplename, Ag_pos_cells$barcode, sep = "_")
CD8_modality_integrate_cluster@meta.data$max_Ag_gt_50 <- Ag_pos_cells[match(rownames(CD8_modality_integrate_cluster@meta.data),Ag_pos_cells$sample_id_barcode),"max_Ag"]
CD8_modality_integrate_cluster@meta.data$max_Ag_gt_50 <- gsub(", ","_",CD8_modality_integrate_cluster@meta.data$max_Ag_gt_50)

CD8_modality_integrate_cluster@meta.data$Ag_UMI <- Ag_pos_cells[match(rownames(CD8_modality_integrate_cluster@meta.data),Ag_pos_cells$sample_id_barcode),"Ag_UMI"]
CD8_modality_integrate_cluster@meta.data$Con_UMI <- Ag_pos_cells[match(rownames(CD8_modality_integrate_cluster@meta.data),Ag_pos_cells$sample_id_barcode),"Con_UMI"]
CD8_modality_integrate_cluster@meta.data$max_specificity_score <- Ag_pos_cells[match(rownames(CD8_modality_integrate_cluster@meta.data),Ag_pos_cells$sample_id_barcode),"max_specificity_score"]
CD8_modality_integrate_cluster@meta.data$sample_Ag <- paste(CD8_modality_integrate_cluster@meta.data$samplename, CD8_modality_integrate_cluster@meta.data$max_Ag_gt_50,sep = "_")

CD8_modality_integrate_cluster@meta.data$max_specificity_score <- as.numeric(CD8_modality_integrate_cluster@meta.data$max_specificity_score)
CD8_modality_integrate_cluster@meta.data$Con_UMI <- as.numeric(CD8_modality_integrate_cluster@meta.data$Con_UMI)
CD8_modality_integrate_cluster@meta.data$Ag_UMI <- as.numeric(CD8_modality_integrate_cluster@meta.data$Ag_UMI)

Ag_sample_table <- table(CD8_modality_integrate_cluster@meta.data$seurat_clusters,CD8_modality_integrate_cluster@meta.data$sample_Ag)
write.table(Ag_sample_table, paste(savedir,"Table/sample_antigen_cluster_table.txt",sep = ""),
            quote = F, row.names = T, col.names = T, sep = "\t")

CD8_modality_integrate_cluster@meta.data$Ag_Age <- paste(CD8_modality_integrate_cluster@meta.data$max_Ag_gt_50, CD8_modality_integrate_cluster@meta.data$Age,sep = "_")
Ag_Age2 <- unique(CD8_modality_integrate_cluster@meta.data$Ag_Age)

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/"
pdf(paste(savedir,"UMAP/Ag_cluster.pdf",sep = ""))
DimPlot(CD8_modality_integrate_cluster, group.by = "Ag_Age", reduction = "wnn.umap") + NoLegend()
dev.off()

library(ArchR)
pdf(paste(savedir,"UMAP/max_specificity_score.pdf",sep = ""))
FeaturePlot(CD8_modality_integrate_cluster, features = "max_specificity_score", reduction = "wnn.umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra)
dev.off()

pdf(paste(savedir,"UMAP/Ag_UMI.pdf",sep = ""))
FeaturePlot(CD8_modality_integrate_cluster, features = "Ag_UMI", reduction = "wnn.umap") + scale_color_gradientn(colors = ArchRPalettes$solarExtra, limits=c(0,500))
dev.off()

pdf(paste(savedir,"UMAP/Con_UMI.pdf",sep = ""))
FeaturePlot(CD8_modality_integrate_cluster, features = "Con_UMI", reduction = "wnn.umap") + 
  scale_color_gradientn(colors = ArchRPalettes$solarExtra, limits=c(0,50))
dev.off()

dir.create(paste(savedir,"Vlnplot",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"Vlnplot/max_specificity_score.pdf",sep = ""), width = 12, height = 8)
VlnPlot(CD8_modality_integrate_cluster, features = "max_specificity_score")
dev.off()

pdf(paste(savedir,"Vlnplot/Ag_UMI.pdf",sep = ""))
VlnPlot(CD8_modality_integrate_cluster, features = "Ag_UMI") + ylim(c(50,500))
dev.off()

pdf(paste(savedir,"Vlnplot/Con_UMI.pdf",sep = ""))
VlnPlot(CD8_modality_integrate_cluster, features = "Con_UMI") + ylim(c(5,200))
dev.off()

rm(plot_list)
plot_list <- list()
for (i in 1:length(Ag_Age2)) {
  cellnames <- CD8_modality_integrate_cluster@meta.data[grep(Ag_Age2[i], CD8_modality_integrate_cluster@meta.data$Ag_Age),] %>% rownames()
  p <- DimPlot(CD8_modality_integrate_cluster,
               cells.highlight = cellnames, 
               reduction = "wnn.umap", 
               label = FALSE, cols.highlight = "deeppink2",
               sizes.highlight = 0.5,
               cols = "gray92") + 
    ggtitle(paste(Ag_Age2[i])) + 
    NoLegend() + 
    theme(plot.title = element_text(hjust = 0.2))
  # pdf(paste(savedir,"CD4/CD4_TRB_",TRB_clonal_noNA_gt_2$Var1[i],".pdf",sep = ""))
  # print(p)
  # dev.off()
  plot_list[[i]] <- p
}

pdf(paste(savedir,"UMAP/Ag_gt_50_Age.pdf",sep = ""), width = 5.5, height = 5)
plot_list
dev.off()

CD8_modality_integrate_cluster@meta.data$Ag_group <- NA
CD8_modality_integrate_cluster@meta.data[grep("EBV_B",CD8_modality_integrate_cluster@meta.data$max_Ag_gt_50),"Ag_group"] <- "Lytic"
CD8_modality_integrate_cluster@meta.data[grep("EBV_EB|EBV_LM",CD8_modality_integrate_cluster@meta.data$max_Ag_gt_50),"Ag_group"] <- "Latent"
CD8_modality_integrate_cluster@meta.data[grep("VZV",CD8_modality_integrate_cluster@meta.data$max_Ag_gt_50),"Ag_group"] <- "VZV"
CD8_modality_integrate_cluster@meta.data$Ag_group_Age <- paste(CD8_modality_integrate_cluster@meta.data$Ag_group,CD8_modality_integrate_cluster$Age,sep="_")
# group_Ag <- grep("NA_",unique(CD8_modality_integrate_cluster@meta.data$Ag_group_Age),invert = TRUE,value = TRUE)
group_Ag <- unique(CD8_modality_integrate_cluster@meta.data$Ag_group_Age)

rm(plot_list)
plot_list <- list()
for (i in 1:length(group_Ag)) {
  cellnames <- CD8_modality_integrate_cluster@meta.data[grep(group_Ag[i], CD8_modality_integrate_cluster@meta.data$Ag_group_Age),] %>% rownames()
  p <- DimPlot(CD8_modality_integrate_cluster,
               cells.highlight = cellnames,
               reduction = "wnn.umap",
               label = FALSE, cols.highlight = "deeppink2",
               sizes.highlight = 0.5,
               cols = "gray92") + 
    ggtitle(paste(group_Ag[i])) + 
    NoLegend() + 
    theme(plot.title = element_text(hjust = 0.2))
  # pdf(paste(savedir,"CD4/CD4_TRB_",TRB_clonal_noNA_gt_2$Var1[i],".pdf",sep = ""))
  # print(p)
  # dev.off()
  plot_list[[i]] <- p
}

pdf(paste(savedir,"UMAP/Ag_group_gt_50_Age.pdf",sep = ""), width = 5.5, height = 5)
plot_list
dev.off()

# require(gridExtra)
pdf(paste(savedir,"UMAP/Ag_group_gt_50_Age.pdf",sep = ""), width = 12, height = 18)
grid.arrange(plot_list[[1]], plot_list[[5]], 
             plot_list[[2]],plot_list[[4]],
             plot_list[[3]],plot_list[[6]], 
             nrow = 3, ncol = 2)
dev.off()

pdf(paste(savedir,"UMAP/Ag_group_gt_50_Age.pdf",sep = ""), width = 11, height = 5)
DimPlot(CD8_modality_integrate_cluster, split.by = "Age", group.by = "Age", ncol = 2, reduction = "wnn.umap")
dev.off()


# CD4_genes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run01/Analysis/downstream/Table/CD4_genes")[,1]
# dir.create(paste(savedir,"dotplot",sep = ""),showWarnings = FALSE)
# pdf(paste(savedir,"dotplot/CD4_genes_CD4_dotplot.pdf",sep = ""))
# DotPlot(CD4_modality_integrate_cluster,"RNA",CD4_genes) + coord_flip()
# dev.off()
# 
# CD8_genes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run01/Analysis/downstream/Table/CD8_genes")[,1]
# dir.create(paste(savedir,"dotplot",sep = ""),showWarnings = FALSE)
# pdf(paste(savedir,"dotplot/CD8_genes_CD4_dotplot.pdf",sep = ""), width = 7, height = 9)
# DotPlot(CD4_modality_integrate_cluster,"RNA",CD8_genes) + coord_flip()
# dev.off()
# 
# ### exhauster genes
# ex_genes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run01/Analysis/downstream/Table/exhausted_genes_2")
# ex_genes_up <- toupper(ex_genes[,1])
# pdf(paste(savedir,"dotplot/exhausted_CD4_dotplot.pdf",sep = ""), height = 10, width = 8)
# DotPlot(CD4_modality_integrate_cluster,"RNA",ex_genes_up) + coord_flip()
# dev.off()
# 
# DefaultAssay(CD4_modality_integrate_cluster) <- "MAGIC_RNA"
# genes <- c("IKZF2","FOXP3","TCF7","LEF1","IRF2","RUNX3","EOMES","TBX21","PRDM1","XBP1",
#            "KLF6","NOTCH1","SLAMF6","FOSL2","FOS","RORC","JUN","JUND","JUNB","STAT3",
#            "Zeb2","ID2","Tox","Myb","TCF7","ID3","BCL6","myc","IL7R","CD84","IL2RB")
# genes <- toupper(genes)
# genes <- rownames(CD4_modality_integrate_cluster)[match(toupper(genes),rownames(CD4_modality_integrate_cluster))]

# genes <- c("IL7R", "CCR7", "CD45RA", "CD45RO", "PD1", "TIM3")
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
# rm(plot_list)
# plot_list <- list()
# for (i in 1:length(genes)) {
#   p <- featureplot_front(CD4_modality_integrate_cluster, genes[i],
#                          reduction = "wnn.umap", x="wnnUMAP_1", y="wnnUMAP_2",size=0.2) +
#     scale_color_gradientn(colours = ArchRPalettes$solarExtra)
#   plot_list[[i]] <- p
# }
# 
# dir.create(paste(savedir,"featureplot/",sep =""),showWarnings = FALSE)
# pdf(paste(savedir,"featureplot/CD4_Goronzy.pdf",sep = ""), width = 4.5, height = 4)
# print(plot_list)
# dev.off()
# 
# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/Gene_score.R")
# CD4_modality_integrate_cluster <- genescore_abhinav(obj = CD4_modality_integrate_cluster,assay = "RNA",
#                                                     savedir = savedir,
#                                                     cores = 5,objname = "CD4")
# 
# 

### Performing PCA for the Pseudobulk
library(Seurat)
library(dplyr)
setwd("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/regress_TCR/")
# CD4 <- readRDS("saveRDS_obj/CD4_modality_integrate_cluster_0.8.RDS")
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/regress_TCR/"
cell_freq = 20
remove_samples = NULL
grouping_by = "max_Ag_gt_50"
column_name_split <- c("samplename","Ag","age_number","gender","Run", "Ag2")

print("all")
dir.create(savedir,showWarnings = FALSE)
savedir2 <- paste(savedir,"pseudobulk/",sep = "")
dir.create(savedir2, showWarnings = FALSE)

DefaultAssay(obj) <- "RNA"
cell_samples <- table(obj@meta.data[,paste(sample_col,grouping_by,sep = "_")]) %>% as.data.frame()
sample_low <- cell_samples[cell_samples$Freq < cell_freq,1]
sample_low <- gsub("_","_",sample_low)
sample_low_2 <- c(sample_low,remove_samples)

message(paste("Removing this sample",sample_low_2,"\n",sep = " "))
Treg_cts <- AggregateExpression(obj, group.by = paste(sample_col,grouping_by,sep = "_"),
                                assays = 'RNA', slot = "counts", return.seurat = FALSE)

Treg_cts_2 <- as.data.frame(Treg_cts$RNA)
if(length(sample_low_2) != 0){
  index <- grep(paste(sample_low_2,collapse="|"),colnames(Treg_cts_2), invert=TRUE)
  Treg_cts_3 <- Treg_cts_2[,index]
} else {
  Treg_cts_3 <- Treg_cts_2
}

Treg_metadata <- data.frame(samples = colnames(Treg_cts_3))
# Function to split and convert to dataframe
Treg_metadata$sample_id = gsub("_.*","",Treg_metadata$samples)
Treg_metadata$Run = gsub(".*.F_|.*.M_","",Treg_metadata$samples) %>% gsub("_.*.","",.)
Treg_metadata$Age <- gsub("[0-9]","",Treg_metadata$sample_id)
Treg_metadata$gender <- gsub(".*._M_.*.","M",Treg_metadata$samples) %>% gsub(".*._F_.*.","F",.)
Treg_metadata$Ag <- gsub(".*._","",Treg_metadata$samples)

stopifnot(all(Treg_metadata$samples == colnames(Treg_cts_3))) # if TRUE move forward

Treg_dds <- DESeqDataSetFromMatrix(countData = Treg_cts_3,
                                   colData = Treg_metadata,
                                   design = ~Run + Ag)

keep = filterByExpr(Treg_dds, group=colData(Treg_dds)[,"samples"], min.count=10)
table(keep)
Treg_dds_2 <- Treg_dds[keep,]

rld<-vst(Treg_dds_2)
PCA <- DESeq2::plotPCA(object = rld,
                       intgroup = colnames(Treg_metadata),
                       returnData=TRUE,ntop = 5000)
percentVar <- round(100 * attr(PCA, "percentVar"))

savedir4 <- paste(savedir2,"PCA/",sep = "")
dir.create(savedir4, showWarnings = FALSE)

PCA[grep("B",PCA$Ag),"Ag_group"] <- "Lytic"
PCA[grep("EB|LM",PCA$Ag),"Ag_group"] <- "Latent"
PCA[grep("I",PCA$Ag),"Ag_group"] <- "VZV"

group_centers <- aggregate(cbind(PC1, PC2) ~ Ag, data = as.data.frame(PCA), FUN = mean)

PC1 <-ggplot(PCA, aes_string("PC1", "PC2", label = "Ag" , shape = "Age", color="Ag")) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=3) +
  # annotate("text", x = group_centers$PC1, y = group_centers$PC2, label = group_centers$Ag, vjust = -1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir4,"PCA_Ag_cell_gt_20_samples_2.pdf",sep = ""),width = 9, height = 8)
print(PC1)
dev.off()

PC_Ag_group <- PC1 + facet_grid(~Ag_group)
pdf(paste(savedir4,"PCA_Ag_cell_gt_20_samples_sep_Ag_group.pdf",sep = ""),width = 20, height = 8)
print(PC_Ag_group)
dev.off()

PCA$Age_Ag <- paste(PCA$Age, PCA$Ag, sep = "_")
PCA$Age_Ag <- factor(PCA$Age_Ag, levels = rev(c("O_BMLF1","Y_BMLF1","O_BRLF1","Y_BRLF1","O_BMRF1",
                                            "Y_BMRF1","O_EBNA3C","Y_EBNA3C","O_EBNA1","Y_EBNA1",
                                            "O_LMP2","O_LMP1","Y_LMP1","Y_LMP2","O_IE62","Y_IE62")))

p <- ggplot(PCA, aes(PC1, y=Age_Ag, fill = Ag_group)) + geom_boxplot() + geom_point() + theme_bw()
pdf(paste(savedir4,"PCA_Ag_cell_gt_20_samples_sep_Ag_group_box.pdf",sep = ""),width = 8, height = 8)
p
dev.off()

p <- ggplot(PCA, aes(PC1, y=Age_Ag, fill = Ag_group), color = "black") +
  # geom_bar(position=position_dodge(), stat="identity") +
  geom_bar(position = "dodge", stat = "summary", fun.x = "mean") +
  # geom_errorbar(aes(xmin=PC1-se, xmax=PC1+se)) +
  theme_bw()


std_PC1 <- vector()
mean_PC1 <- vector()

Age_Ag <- levels(PCA$Age_Ag)
for (i in 1:length(Age_Ag)) {
  mean_PC1[i] <- mean(PCA[grep(Age_Ag[i], PCA$Age_Ag),"PC1"])
  std_PC1[i] <- sd(PCA[grep(Age_Ag[i], PCA$Age_Ag),"PC1"])
}

std_PC1[is.na(std_PC1)] <- 0

df <- data.frame(matrix(nrow = length(Age_Ag), ncol = 3))
colnames(df) <- c("Age_Ag","mean_PC1","std_PC1")

df$Age_Ag <- Age_Ag
df$mean_PC1 <- mean_PC1
df$std_PC1 <- std_PC1

df[grep("_B",df$Age_Ag),"Ag_group"] <- "Lytic"
df[grep("_EB|LM",df$Age_Ag),"Ag_group"] <- "Latent"
df[grep("_I",df$Age_Ag),"Ag_group"] <- "VZV"

df$Age_Ag <- factor(df$Age_Ag, levels = rev(c("O_BMLF1","Y_BMLF1","O_BRLF1","Y_BRLF1","O_BMRF1",
                                              "Y_BMRF1","O_EBNA3C","Y_EBNA3C","O_EBNA1","Y_EBNA1",
                                              "O_LMP2","O_LMP1","Y_LMP1","Y_LMP2","O_IE62","Y_IE62")))

p2 <- ggplot(df, aes(mean_PC1, y=Age_Ag),color = "black") +
  # geom_bar(position=position_dodge(), stat="identity") +
  geom_bar(aes( fill = Ag_group), position = "dodge", stat = "identity") + geom_point() +
  geom_errorbar(aes(x=mean_PC1, y=Age_Ag, xmin=mean_PC1-std_PC1, xmax=mean_PC1+std_PC1)) +
  theme_bw()

pdf(paste(savedir4,"PC1_Ag_cell_gt_20_samples_sep_Ag_group_bar_box.pdf",sep = ""),width = 16, height = 8)
p1+p2
dev.off()


pdf(paste(savedir4,"PCA_Ag_cell_gt_20_samples_sep_Ag_group_bar.pdf",sep = ""),width = 8, height = 8)
p
dev.off()


p1 <- ggplot(PCA, aes(PC2, y=Age_Ag, fill = Ag_group)) + geom_boxplot() + geom_point() + theme_bw()

# p2 <- ggplot(PCA, aes(PC2, y=Age_Ag, fill = Ag_group), color = "black") +
#   # geom_bar(position=position_dodge(), stat="identity") +
#   geom_bar(position = "dodge", stat = "summary", fun.x = "mean") +
#   # geom_errorbar(aes(xmin=PC2-se, xmax=PC2+se)) +
#   theme_bw()

std_PC2 <- vector()
mean_PC2 <- vector()

Age_Ag <- levels(PCA$Age_Ag)
for (i in 1:length(Age_Ag)) {
  mean_PC2[i] <- mean(PCA[grep(Age_Ag[i], PCA$Age_Ag),"PC2"])
  std_PC2[i] <- sd(PCA[grep(Age_Ag[i], PCA$Age_Ag),"PC2"])
}

std_PC2[is.na(std_PC2)] <- 0

df <- data.frame(matrix(nrow = length(Age_Ag), ncol = 3))
colnames(df) <- c("Age_Ag","mean_PC2","std_PC2")

df$Age_Ag <- Age_Ag
df$mean_PC2 <- mean_PC2
df$std_PC2 <- std_PC2

df[grep("_B",df$Age_Ag),"Ag_group"] <- "Lytic"
df[grep("_EB|LM",df$Age_Ag),"Ag_group"] <- "Latent"
df[grep("_I",df$Age_Ag),"Ag_group"] <- "VZV"

df$Age_Ag <- factor(df$Age_Ag, levels = rev(c("O_BMLF1","Y_BMLF1","O_BRLF1","Y_BRLF1","O_BMRF1",
                                 "Y_BMRF1","O_EBNA3C","Y_EBNA3C","O_EBNA1","Y_EBNA1",
                                 "O_LMP2","O_LMP1","Y_LMP1","Y_LMP2","O_IE62","Y_IE62")))

p2 <- ggplot(df, aes(mean_PC2, y=Age_Ag),color = "black") +
  # geom_bar(position=position_dodge(), stat="identity") +
  geom_bar(aes( fill = Ag_group), position = "dodge", stat = "identity") + geom_point() +
  geom_errorbar(aes(x=mean_PC2, y=Age_Ag, xmin=mean_PC2-std_PC2, xmax=mean_PC2+std_PC2)) +
  theme_bw()

pdf(paste(savedir4,"PC2_Ag_cell_gt_20_samples_sep_Ag_group_bar_box.pdf",sep = ""),width = 16, height = 8)
p1+p2
dev.off()

write.table(PCA, paste(paste(savedir,"PCA.txt",sep = "")), quote = F, row.names = T, col.names = T, sep = "\t")


### If we want to have only one group
group_centers <- aggregate(cbind(PC1, PC2) ~ Ag, data = as.data.frame(PCA), FUN = mean)


PC1 <-ggplot(PCA, aes_string("PC1", "PC2", shape = "Age", color="Ag_group")) +
  geom_point(size=5, alpha=0.7) +
  annotate("text", x = group_centers$PC1, y = group_centers$PC2, label = group_centers$Ag, vjust = -1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir4,"PCA_Ag_cell_gt_20_samples_2.pdf",sep = ""),width = 9, height = 8)
print(PC1)
dev.off()

write.table(PCA , paste(savedir4,"PCA_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".txt",sep = ""),
            quote = F, row.names = T, col.names = T, sep = "\t")

assay(rld) <- limma::removeBatchEffect(assay(rld),
                                       batch=rld@colData[,batch_col])

PCA_batch <- DESeq2::plotPCA(object = rld,
                             intgroup = colnames(Treg_metadata),
                             returnData=TRUE,ntop = 5000)

percentVar <- round(100 * attr(PCA_batch, "percentVar"))

group_centers <- aggregate(cbind(PC1, PC2) ~ Ag, data = as.data.frame(PCA_batch), FUN = mean)

PCA_batch[grep("B",PCA_batch$Ag),"Ag_group"] <- "Lytic"
PCA_batch[grep("EB|LM",PCA_batch$Ag),"Ag_group"] <- "Latent"
PCA_batch[grep("I",PCA_batch$Ag),"Ag_group"] <- "VZV"

PC1 <-ggplot(PCA_batch, aes_string("PC1", "PC2", shape = "Age", color="Ag_group")) +
  geom_point(size=5, alpha=0.7) +
  annotate("text", x = group_centers$PC1, y = group_centers$PC2, label = group_centers$Ag, vjust = -1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir4,"PCA_Ag_cell_gt_20_samples_2_remove_batch.pdf.pdf",sep = ""),width = 9, height = 8)
print(PC1)
dev.off()

write.table(PCA_batch, paste(savedir4,"PCA_remove_batch_effect.txt",sep = ""),
            quote = F, row.names = T, col.names = T, sep = "\t")


#### Genotyping Demultiplex #####
# singularity exec -B /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav:/mnt /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/souporcell_latest.sif \
# souporcell_pipeline.py -i /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01A/outs/possorted_genome_bam.bam \
# -b /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01A/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
# -f /mnt/Resources/refdata-gex-GRCh38-2020-A/fasta/genome.fa -t 8 -o /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01A/post_processing/soupercell_VCF -k 4
# 
# singularity exec -B /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav:/mnt /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/souporcell_latest.sif \
# souporcell_pipeline.py -i /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01B/outs/possorted_genome_bam.bam \
# -b /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01B/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
# -f /mnt/Resources/refdata-gex-GRCh38-2020-A/fasta/genome.fa -t 8 -o /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01B/post_processing/soupercell_VCF -k 4

### Identifying the number of cells for each individuals identified
# [m256617@mforgehn1 genotype_demux]$ grep -wFf Y1_HTO.txt Exp01B_singlet.tsv | awk '{print $3}' |  sort | uniq -c
# 6 0
# 3477 1
# 42 2
# 9 3
# [m256617@mforgehn1 genotype_demux]$ grep -wFf Y1_HTO.txt Exp01A_singlet.tsv | awk '{print $3}' |  sort | uniq -c
# 4 0
# 9 1
# 33 2
# 3653 3
# [m256617@mforgehn1 genotype_demux]$ grep -wFf Y2_HTO.txt Exp01A_singlet.tsv | awk '{print $3}' |  sort | uniq -c
# 4 1
# 3015 2
# 6 3
# [m256617@mforgehn1 genotype_demux]$ grep -wFf Y2_HTO.txt Exp01B_singlet.tsv | awk '{print $3}' |  sort | uniq -c
# 11 1
# 2938 2
# 8 3
# [m256617@mforgehn1 genotype_demux]$ grep -wFf O1_HTO.txt Exp01B_singlet.tsv | awk '{print $3}' |  sort | uniq -c
# 3502 0
# 8 1
# 25 2
# 8 3
# [m256617@mforgehn1 genotype_demux]$ grep -wFf O1_HTO.txt Exp01A_singlet.tsv | awk '{print $3}' |  sort | uniq -c
# 3725 0
# 10 1
# 32 2
# 9 3
# [m256617@mforgehn1 genotype_demux]$ grep -wFf O2_HTO.txt Exp01A_singlet.tsv | awk '{print $3}' |  sort | uniq -c
# 2 0
# 2751 1
# 24 2
# 10 3
# [m256617@mforgehn1 genotype_demux]$ grep -wFf O2_HTO.txt Exp01B_singlet.tsv | awk '{print $3}' |  sort | uniq -c
# 4 0
# 7 1
# 18 2
# 2532 3

### Thinking to do confusion matrix and the Venn Diagram
### Lets automate this process
Exp1_HTO <- c("Y1","Y2","O1","O2")
Exp1A_gen <- c(3,2,0,1)
Exp1B_gen <- c(1,2,0,3)

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01_First_4/post_processing/genotype_demux/"
Exp01A <- read_tsv(paste(savedir,"Exp01A_singlet.tsv",sep = ""), col_names = FALSE)

library(stringr)
patterns = Exp1A_gen
replacements = Exp1_HTO
# names(replacement) <- patterns
Exp01A$sample <- Exp01A$X3
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  Exp01A$sample <- str_replace_all(Exp01A$sample, pattern, replacements[i])
}

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01_First_4/post_processing/genotype_demux/"
Exp01B <- read_tsv(paste(savedir,"Exp01B_singlet.tsv",sep = ""), col_names = FALSE)

library(stringr)
patterns = Exp1B_gen
replacements = Exp1_HTO
# names(replacement) <- patterns
Exp01B$sample <- Exp01B$X3
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  Exp01B$sample <- str_replace_all(Exp01B$sample, pattern, replacements[i])
}

Exp01A_req <- as.data.frame(Exp01A[,c("X1","sample")])
Exp01B_req <- as.data.frame(Exp01B[,c("X1","sample")])

Exp01_req <- rbind(Exp01A_req, Exp01B_req)
colnames(Exp01_req) <- c("barcode","sample")

HTO_filename <- list.files(savedir, pattern = "_HTO.txt")
samplename <- gsub("_HTO.txt","",HTO_filename)
for (i in 1:length(HTO_filename)) {
  HTO_files <- read.table(paste(savedir,HTO_filename[i],sep = ""), header = FALSE)
  HTO_files$sample <- samplename[i]
  assign(paste(samplename[i],"_file",sep = ""), HTO_files)
}

#### Put the number of files to merge
sample_file <- ls(pattern = "[0-9]_file")
merge_samples <- rbind(get(sample_file[1]),get(sample_file[2]),get(sample_file[3]),get(sample_file[4]))
colnames(merge_samples) <- c("barcode","sample")

# Assuming merge_samples is the predicted dataframe and Exp01_req is the actual dataframe

# Merge the two dataframes based on the 'barcode' column
merged_df <- merge(merge_samples, Exp01_req, by = "barcode", suffixes = c("_predicted", "_actual"))

# Create a confusion matrix
conf_matrix <- table(merged_df$sample_actual, merged_df$sample_predicted)

# Print the confusion matrix
conf_matrix_df <- as.data.frame.matrix(conf_matrix)
conf_matrix_df$sample <- rownames(conf_matrix_df) 
conf_matrix_df_melted <- melt(conf_matrix_df, by = "sample")
colnames(conf_matrix_df_melted) <- c("HTO","Geno","cells")

## Finding the sensitivity and specificity
# Merge the two dataframes based on the 'barcode' column
merged_df <- merge(merge_samples, Exp01_req, by = "barcode", suffixes = c("_predicted", "_actual"))

# Ensure the 'sample_actual' and 'sample_predicted' columns are factors with the same levels
merged_df$sample_actual <- factor(merged_df$sample_actual, levels = unique(c(merged_df$sample_actual, merged_df$sample_predicted)))
merged_df$sample_predicted <- factor(merged_df$sample_predicted, levels = unique(c(merged_df$sample_actual, merged_df$sample_predicted)))

# Load caret package if not already loaded
# install.packages("caret")
library(caret)

# Create a confusion matrix with caret
merged_df$sample_predicted <- factor(merged_df$sample_predicted) 
caret_conf_matrix <- caret::confusionMatrix(merged_df$sample_predicted, merged_df$sample_actual)
conf_matrix_by_class <- as.data.frame(caret_conf_matrix$byClass)

sensitivity <- mean(conf_matrix_by_class$Sensitivity)
specificity <- mean(conf_matrix_by_class$Specificity)
accuracy <- mean(conf_matrix_by_class$`Balanced Accuracy`)

### Making a plot
library(ggplot2)
p <- ggplot(data = conf_matrix_df_melted, aes(x = HTO, y = Geno, fill = cells, label = cells)) +
  geom_tile(color = "white") +
  geom_text(vjust = 1) +
  scale_fill_gradientn(colors = ArchRPalettes$solarExtra) +
  theme_minimal() +
  labs(x = "HTO", y = "Geno") +
  ggtitle("Exp01 HTO and Genotype Matrix") +
  theme(legend.position = "right") +  # Move the legend to the right
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +  # Position legend at the top
  annotate("text", x = 0.5, y = 4.55, label = paste("Sensitivity:", round(sensitivity, 2), 
                                                    "  Specificity:", round(specificity, 2),
                                                    "  Accuracy:", round(accuracy, 2)),
           color = "black", size = 3, hjust = 0)

pdf(paste(savedir,"Exp01_HTO_geno_confusion_matx.pdf",sep = ""), width = 4.5, height = 4.5)
p
dev.off()


### Making the Venn diagram for each celltypes
library(dplyr)

diff2 <- list(Geno=Exp01_req$barcode, 
              HTO=merge_samples$barcode)

ggvenn(diff2, text_size = 8)

pdf(paste(savedir,"Exp01_venn_all.pdf",sep = ""))
print(ggvenn(diff2, text_size = 8) + ggtitle("Exp01 all") + theme(plot.title = element_text(hjust = 0.5)))
dev.off()


for (i in 1:length(Exp1_HTO)) {
  geno_sample <- Exp01_req[grep(Exp1_HTO[i],Exp01_req$sample),"barcode"]
  HTO_sample <- merge_samples[grep(Exp1_HTO[i],merge_samples$sample),"barcode"]
  
  diff2 <- list(Geno=geno_sample, 
                HTO=HTO_sample)
  
  pdf(paste(savedir,"Exp01_venn_",Exp1_HTO[i],".pdf",sep = ""))
  print(ggvenn(diff2, text_size = 8) + ggtitle(paste("Exp01",Exp1_HTO[i])) + theme(plot.title = element_text(hjust = 0.5))) 
  dev.off()
}

### Since we need to run the cellranger on the barcode extracted
write.table(merge_samples, paste(savedir,"HTO_multiplexed.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(Exp01_req, paste(savedir,"Genotype_multiplexed.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")


### Cell same between Exp01A and Exp01B
library(dplyr)
library(ggvenn)

diff2 <- list(Exp01A=Exp01A_req$X1,
              Exp01B=Exp01B_req$X1)

# ggvenn(diff2, text_size = 8)

pdf(paste(savedir,"Exp01A_and_B_venn_all.pdf",sep = ""))
print(ggvenn(diff2, text_size = 8) + ggtitle("Exp01 A and B") + theme(plot.title = element_text(hjust = 0.5)))
dev.off()


### Exp02 
# singularity exec -B /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav:/mnt /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/souporcell_latest.sif \
# souporcell_pipeline.py -i /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02A/outs/possorted_genome_bam.bam \
# -b /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02A/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
# -f /mnt/Resources/refdata-gex-GRCh38-2020-A/fasta/genome.fa -t 8 -o /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02A/post_processing/soupercell_VCF -k 5
# 
# singularity exec -B /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav:/mnt /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/souporcell_latest.sif \
# souporcell_pipeline.py -i /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02B/outs/possorted_genome_bam.bam \
# -b /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02B/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
# -f /mnt/Resources/refdata-gex-GRCh38-2020-A/fasta/genome.fa -t 8 -o /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02B/post_processing/soupercell_VCF -k 5

### Thinking to do confusion matrix and the Venn Diagram
### Lets automate this process
library(readr)
library(caret)
library(dplyr)
library(ggplot2)
Exp2_HTO <- c("Y3","Y4","Y5","O3","O4")
Exp2A_gen <- c(3,4,0,1,2)
Exp2B_gen <- c(3,2,0,4,1)

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02_First_4/post_processing/genotype_demux/"
Exp02A <- read_tsv(paste(savedir,"Exp02A_singlet.tsv",sep = ""), col_names = FALSE)

library(stringr)
patterns = Exp2A_gen
replacements = Exp2_HTO
# names(replacement) <- patterns
Exp02A$sample <- Exp02A$X3
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  Exp02A$sample <- str_replace_all(Exp02A$sample, pattern, replacements[i])
}

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02_First_4/post_processing/genotype_demux/"
Exp02B <- read_tsv(paste(savedir,"Exp02B_singlet.tsv",sep = ""), col_names = FALSE)

library(stringr)
patterns = Exp2B_gen
replacements = Exp2_HTO
# names(replacement) <- patterns
Exp02B$sample <- Exp02B$X3
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  Exp02B$sample <- str_replace_all(Exp02B$sample, pattern, replacements[i])
}

Exp02A_req <- as.data.frame(Exp02A[,c("X1","sample")])
Exp02B_req <- as.data.frame(Exp02B[,c("X1","sample")])

Exp02_req <- rbind(Exp02A_req, Exp02B_req)
colnames(Exp02_req) <- c("barcode","sample")

HTO_filename <- list.files(savedir, pattern = "_HTO.txt")
samplename <- gsub("_HTO.txt","",HTO_filename)
for (i in 1:length(HTO_filename)) {
  HTO_files <- read.table(paste(savedir,HTO_filename[i],sep = ""), header = FALSE)
  HTO_files$sample <- samplename[i]
  assign(paste(samplename[i],"_file",sep = ""), HTO_files)
}

#### Put the number of files to merge
sample_file <- ls(pattern = "[0-9]_file")
merge_samples <- rbind(get(sample_file[1]),get(sample_file[2]),get(sample_file[3]),get(sample_file[4]),get(sample_file[5]))
colnames(merge_samples) <- c("barcode","sample")

# Assuming merge_samples is the predicted dataframe and Exp02_req is the actual dataframe

# Merge the two dataframes based on the 'barcode' column
merged_df <- merge(merge_samples, Exp02_req, by = "barcode", suffixes = c("_predicted", "_actual"))

# Create a confusion matrix
conf_matrix <- table(merged_df$sample_actual, merged_df$sample_predicted)

# Print the confusion matrix
library(reshape2)
conf_matrix_df <- as.data.frame.matrix(conf_matrix)
conf_matrix_df$sample <- rownames(conf_matrix_df) 
conf_matrix_df_melted <- melt(conf_matrix_df, by = "sample")
colnames(conf_matrix_df_melted) <- c("HTO","Geno","cells")

## Finding the sensitivity and specificity
# Merge the two dataframes based on the 'barcode' column
merged_df <- merge(merge_samples, Exp02_req, by = "barcode", suffixes = c("_predicted", "_actual"))

# # Ensure the 'sample_actual' and 'sample_predicted' columns are factors with the same levels
# merged_df$sample_actual <- factor(merged_df$sample_actual, levels = unique(c(merged_df$sample_actual, merged_df$sample_predicted)))
# merged_df$sample_predicted <- factor(merged_df$sample_predicted, levels = unique(c(merged_df$sample_actual, merged_df$sample_predicted)))

# Load caret package if not already loaded
# install.packages("caret")
library(caret)

# Create a confusion matrix with caret
merged_df$sample_predicted <- as.factor(merged_df$sample_predicted)
merged_df$sample_actual <- as.factor(merged_df$sample_actual)
caret_conf_matrix <- confusionMatrix(merged_df$sample_predicted, merged_df$sample_actual)
conf_matrix_by_class <- as.data.frame(caret_conf_matrix$byClass)

sensitivity <- mean(conf_matrix_by_class$Sensitivity)
specificity <- mean(conf_matrix_by_class$Specificity)
accuracy <- mean(conf_matrix_by_class$`Balanced Accuracy`)

### Making a plot
library(ggplot2)
library(ArchR)
p <- ggplot(data = conf_matrix_df_melted, aes(x = HTO, y = Geno, fill = cells, label = cells)) +
  geom_tile(color = "white") +
  geom_text(vjust = 1) +
  scale_fill_gradientn(colors = ArchRPalettes$solarExtra) +
  theme_minimal() +
  labs(x = "HTO", y = "Geno") +
  ggtitle("Exp02 HTO and Genotype Matrix") +
  theme(legend.position = "right") +  # Move the legend to the right
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +  # Position legend at the top
  annotate("text", x = 0.5, y = 5.55, label = paste("Sensitivity:", round(sensitivity, 2), 
                                                    "  Specificity:", round(specificity, 2),
                                                    "  Accuracy:", round(accuracy, 2)),
           color = "black", size = 3, hjust = 0)

pdf(paste(savedir,"Exp02_HTO_geno_confusion_matx.pdf",sep = ""), width = 4.5, height = 4.5)
p
dev.off()


### Making the Venn diagram for each celltypes
library(dplyr)
library(ggvenn)

diff2 <- list(Geno=Exp02_req$barcode, 
              HTO=merge_samples$barcode)

ggvenn(diff2, text_size = 8)

pdf(paste(savedir,"Exp02_venn_all.pdf",sep = ""))
print(ggvenn(diff2, text_size = 8) + ggtitle("Exp02 all") + theme(plot.title = element_text(hjust = 0.5)))
dev.off()


for (i in 1:length(Exp2_HTO)) {
  geno_sample <- Exp02_req[grep(Exp2_HTO[i],Exp02_req$sample),"barcode"]
  HTO_sample <- merge_samples[grep(Exp2_HTO[i],merge_samples$sample),"barcode"]
  
  diff2 <- list(Geno=geno_sample, 
                HTO=HTO_sample)
  
  pdf(paste(savedir,"Exp02_venn_",Exp2_HTO[i],".pdf",sep = ""))
  print(ggvenn(diff2, text_size = 8) + ggtitle(paste("Exp02",Exp2_HTO[i])) + theme(plot.title = element_text(hjust = 0.5))) 
  dev.off()
}


### Making the Venn diagram for each celltypes
library(dplyr)
library(ggvenn)

diff2 <- list(Geno=Exp02_req$barcode,
              HTO=merge_samples$barcode)

ggvenn(diff2, text_size = 8)

pdf(paste(savedir,"Exp02_venn_all.pdf",sep = ""))
print(ggvenn(diff2, text_size = 8) + ggtitle("Exp02 all") + theme(plot.title = element_text(hjust = 0.5)))
dev.off()

for (i in 1:length(Exp3_HTO)) {
  geno_sample <- Exp02_req[grep(Exp3_HTO[i],Exp02_req$sample),"barcode"]
  HTO_sample <- merge_samples[grep(Exp3_HTO[i],merge_samples$sample),"barcode"]
  
  diff2 <- list(Geno=geno_sample,
                HTO=HTO_sample)
  
  pdf(paste(savedir,"Exp02_venn_",Exp3_HTO[i],".pdf",sep = ""))
  print(ggvenn(diff2, text_size = 8) + ggtitle(paste("Exp02",Exp3_HTO[i])) + theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
}

### Since we need to run the cellranger on the barcode extracted
write.table(merge_samples, paste(savedir,"HTO_multiplexed.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(Exp02_req, paste(savedir,"Genotype_multiplexed.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

### Cell same between Exp02A and Exp02B
library(dplyr)
library(ggvenn)

diff2 <- list(Exp02A=Exp02A_req$X1,
              Exp02B=Exp02B_req$X1)

# ggvenn(diff2, text_size = 8)

pdf(paste(savedir,"Exp02A_and_B_venn_all.pdf",sep = ""))
print(ggvenn(diff2, text_size = 8) + ggtitle("Exp02 A and B") + theme(plot.title = element_text(hjust = 0.5)))
dev.off()

# singularity exec -B /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav:/mnt /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/souporcell_latest.sif \
# souporcell_pipeline.py -i /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03A/outs/possorted_genome_bam.bam \
# -b /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03A/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
# -f /mnt/Resources/refdata-gex-GRCh38-2020-A/fasta/genome.fa -t 8 -o /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03A/post_processing/soupercell_VCF -k 5
# 
# singularity exec -B /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav:/mnt /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/tools/souporcell_latest.sif \
# souporcell_pipeline.py -i /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03B/outs/possorted_genome_bam.bam \
# -b /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03B/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
# -f /mnt/Resources/refdata-gex-GRCh38-2020-A/fasta/genome.fa -t 8 -o /mnt/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03B/post_processing/soupercell_VCF -k 5

# cd /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/
# ls *_Exp03/outs/per_sample_outs/*/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz | sed 's:/:\t:g' | awk '{print "zcat "$1"/outs/per_sample_outs/"$1"/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz > /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/genotype_demux/"$1"_HTO.txt"}'
# grep "singlet" /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03B/post_processing/soupercell_VCF/clusters.tsv > /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/genotype_demux/Exp03A_singlet.tsv
# ls *.txt | awk '{print "grep -wFf "$1" Exp03A_singlet.tsv | cut -f3 | sort | uniq -c > "$1"_Exp03A_summary.txt"}' > Exp03A_summary.sh

### Thinking to do confusion matrix and the Venn Diagram
### Lets automate this process
library(readr)
library(caret)
library(dplyr)
library(ggplot2)
Exp3_HTO <- c("O5","O6","O7","Y6","Y7")
Exp3A_gen <- c(4,3,2,0,1)
Exp3B_gen <- c(4,3,1,2,0)

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/genotype_demux/"
Exp03A <- read_tsv(paste(savedir,"Exp03A_singlet.tsv",sep = ""), col_names = FALSE)

library(stringr)
patterns = Exp3A_gen
replacements = Exp3_HTO
# names(replacement) <- patterns
Exp03A$sample <- Exp03A$X3
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  Exp03A$sample <- str_replace_all(Exp03A$sample, pattern, replacements[i])
}

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/genotype_demux/"
Exp03B <- read_tsv(paste(savedir,"Exp03B_singlet.tsv",sep = ""), col_names = FALSE)

library(stringr)
patterns = Exp3B_gen
replacements = Exp3_HTO
# names(replacement) <- patterns
Exp03B$sample <- Exp03B$X3
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  Exp03B$sample <- str_replace_all(Exp03B$sample, pattern, replacements[i])
}

Exp03A_req <- as.data.frame(Exp03A[,c("X1","sample")])
Exp03B_req <- as.data.frame(Exp03B[,c("X1","sample")])

Exp03_req <- rbind(Exp03A_req, Exp03B_req)
colnames(Exp03_req) <- c("barcode","sample")

HTO_filename <- list.files(savedir, pattern = "_HTO.txt")
HTO_filename <- grep("summary.txt",HTO_filename,invert = TRUE,value=TRUE)
samplename <- gsub("_.*.","",HTO_filename)
for (i in 1:length(HTO_filename)) {
  HTO_files <- read.table(paste(savedir,HTO_filename[i],sep = ""), header = FALSE)
  HTO_files$sample <- samplename[i]
  assign(paste(samplename[i],"_file",sep = ""), HTO_files)
}

#### Put the number of files to merge
sample_file <- ls(pattern = "[0-9]_file")
merge_samples <- rbind(get(sample_file[1]),get(sample_file[2]),get(sample_file[3]),get(sample_file[4]),get(sample_file[5]))
colnames(merge_samples) <- c("barcode","sample")

# Assuming merge_samples is the predicted dataframe and Exp03_req is the actual dataframe

# Merge the two dataframes based on the 'barcode' column
merged_df <- merge(merge_samples, Exp03_req, by = "barcode", suffixes = c("_predicted", "_actual"))

# Create a confusion matrix
conf_matrix <- table(merged_df$sample_actual, merged_df$sample_predicted)

# Print the confusion matrix
library(reshape2)
conf_matrix_df <- as.data.frame.matrix(conf_matrix)
conf_matrix_df$sample <- rownames(conf_matrix_df) 
conf_matrix_df_melted <- melt(conf_matrix_df, by = "sample")
colnames(conf_matrix_df_melted) <- c("HTO","Geno","cells")

## Finding the sensitivity and specificity
# Merge the two dataframes based on the 'barcode' column
merged_df <- merge(merge_samples, Exp03_req, by = "barcode", suffixes = c("_predicted", "_actual"))

# # Ensure the 'sample_actual' and 'sample_predicted' columns are factors with the same levels
# merged_df$sample_actual <- factor(merged_df$sample_actual, levels = unique(c(merged_df$sample_actual, merged_df$sample_predicted)))
# merged_df$sample_predicted <- factor(merged_df$sample_predicted, levels = unique(c(merged_df$sample_actual, merged_df$sample_predicted)))

# Load caret package if not already loaded
# install.packages("caret")
library(caret)
# Create a confusion matrix with caret
merged_df$sample_predicted <- as.factor(merged_df$sample_predicted)
merged_df$sample_actual <- as.factor(merged_df$sample_actual)
caret_conf_matrix <- caret::confusionMatrix(merged_df$sample_predicted, merged_df$sample_actual)
conf_matrix_by_class <- as.data.frame(caret_conf_matrix$byClass)

sensitivity <- mean(conf_matrix_by_class$Sensitivity)
specificity <- mean(conf_matrix_by_class$Specificity)
accuracy <- mean(conf_matrix_by_class$`Balanced Accuracy`)

### Making a plot
library(ggplot2)
library(ArchR)
p <- ggplot(data = conf_matrix_df_melted, aes(x = HTO, y = Geno, fill = cells, label = cells)) +
  geom_tile(color = "white") +
  geom_text(vjust = 1) +
  scale_fill_gradientn(colors = ArchRPalettes$solarExtra) +
  theme_minimal() +
  labs(x = "HTO", y = "Geno") +
  ggtitle("Exp03 HTO and Genotype Matrix") +
  theme(legend.position = "right") +  # Move the legend to the right
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +  # Position legend at the top
  annotate("text", x = 0.5, y = 5.55, label = paste("Sensitivity:", round(sensitivity, 2), 
                                                    "  Specificity:", round(specificity, 2),
                                                    "  Accuracy:", round(accuracy, 2)),
           color = "black", size = 3, hjust = 0)

pdf(paste(savedir,"Exp03_HTO_geno_confusion_matx.pdf",sep = ""), width = 4.5, height = 4.5)
p
dev.off()


### Making the Venn diagram for each celltypes
library(dplyr)
library(ggvenn)

diff2 <- list(Geno=Exp03_req$barcode,
              HTO=merge_samples$barcode)

ggvenn(diff2, text_size = 8)

pdf(paste(savedir,"Exp03_venn_all.pdf",sep = ""))
print(ggvenn(diff2, text_size = 8) + ggtitle("Exp03 all") + theme(plot.title = element_text(hjust = 0.5)))
dev.off()

for (i in 1:length(Exp3_HTO)) {
  geno_sample <- Exp03_req[grep(Exp3_HTO[i],Exp03_req$sample),"barcode"]
  HTO_sample <- merge_samples[grep(Exp3_HTO[i],merge_samples$sample),"barcode"]
  
  diff2 <- list(Geno=geno_sample,
                HTO=HTO_sample)
  
  pdf(paste(savedir,"Exp03_venn_",Exp3_HTO[i],".pdf",sep = ""))
  print(ggvenn(diff2, text_size = 8) + ggtitle(paste("Exp03",Exp3_HTO[i])) + theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
}

### Since we need to run the cellranger on the barcode extracted
write.table(merge_samples, paste(savedir,"HTO_multiplexed.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(Exp03_req, paste(savedir,"Genotype_multiplexed.txt",sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

### Cell same between Exp03A and Exp03B
library(dplyr)
library(ggvenn)

diff2 <- list(Exp03A=Exp03A_req$X1,
              Exp03B=Exp03B_req$X1)

# ggvenn(diff2, text_size = 8)

pdf(paste(savedir,"Exp03A_and_B_venn_all.pdf",sep = ""))
print(ggvenn(diff2, text_size = 8) + ggtitle("Exp03 A and B") + theme(plot.title = element_text(hjust = 0.5)))
dev.off()

### Exp03 Ag and Bulk #####
#### Identifying the cut off for the Ag UMI, Bulk UMI and Ag specificity score
### Find out the BEAM-T positive cells 
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/genotype_demux/Exp03_geno/outs/per_sample_outs/Exp03_geno/antigen_analysis/"
Ag_spec <- read.csv(paste(maindir, "antigen_specificity_scores.csv", sep = ""))
barcodes <- unique(Ag_spec$barcode)
  
## Creating an empty dataframe for all Ag information
rm(df)
df <- data.frame(matrix(ncol = 7, nrow = length(barcodes)))
colnames(df) <- c("samplename","barcode","Ag_UMI","Con_UMI","Ag_Con_UMI","max_specificity_score","max_Ag")
df$barcode <- barcodes

#### Which barcode match the samples
Exp03_demux <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/genotype_demux/Exp03_cellbarcode_samplename.txt", 
                          sep = ",", header = TRUE)

df$samplename <- Exp03_demux[match(df$barcode, Exp03_demux$Barcode),"Assignment"]

for (i in 1:length(barcodes)) {
  subset_table <- Ag_spec[grep(barcodes[i], Ag_spec$barcode),c("antigen","antigen_umi","control_umi","antigen_specificity_score")]
  df[i,"Ag_UMI"] <- sum(subset_table$antigen_umi)
  df[i,"Con_UMI"] <- sum(subset_table$control_umi)
  df[i,"Ag_Con_UMI"] <- sum(subset_table$antigen_umi,subset_table$control_umi)
  df[i,"max_specificity_score"] <- max(subset_table[,"antigen_specificity_score"])
  df[i,"max_Ag"] <- paste(subset_table[grep(df[i,"max_specificity_score"], subset_table[,"antigen_specificity_score"]),"antigen"], collapse = ", ")
  }

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/"
write.table(df, paste(savedir,"Table/Exp02_Ag_UMI.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

#### Exp02
### Ag and Bulk #####
#### Identifying the cut off for the Ag UMI, Bulk UMI and Ag specificity score
### Find out the BEAM-T positive cells 
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02_First_4/post_processing/genotype_demux/Exp02_geno/outs/per_sample_outs/Exp02_geno/antigen_analysis/"

Ag_spec <- read.csv(paste(maindir, "antigen_specificity_scores.csv", sep = ""))
barcodes <- unique(Ag_spec$barcode)

## Creating an empty dataframe for all Ag information
rm(df)
df <- data.frame(matrix(ncol = 7, nrow = length(barcodes)))
colnames(df) <- c("samplename","barcode","Ag_UMI","Con_UMI","Ag_Con_UMI","max_specificity_score","max_Ag")
df$barcode <- barcodes

#### Which barcode match the samples
Exp02_demux <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02_First_4/post_processing/genotype_demux/Exp02_cellbarcode_samplename.txt", 
                          sep = ",", header = TRUE)

df$samplename <- Exp02_demux[match(df$barcode, Exp02_demux$Barcode),"Assignment"]

for (i in 1:length(barcodes)) {
  subset_table <- Ag_spec[grep(barcodes[i], Ag_spec$barcode),c("antigen","antigen_umi","control_umi","antigen_specificity_score")]
  df[i,"Ag_UMI"] <- sum(subset_table$antigen_umi)
  df[i,"Con_UMI"] <- sum(subset_table$control_umi)
  df[i,"Ag_Con_UMI"] <- sum(subset_table$antigen_umi,subset_table$control_umi)
  df[i,"max_specificity_score"] <- max(subset_table[,"antigen_specificity_score"])
  df[i,"max_Ag"] <- paste(subset_table[grep(df[i,"max_specificity_score"], subset_table[,"antigen_specificity_score"]),"antigen"], collapse = ", ")
}

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/Table/"
write.table(combined_df_Exp02, paste(savedir,"combined_df_Exp02.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

#### Exp01
### Ag and Bulk #####
#### Identifying the cut off for the Ag UMI, Bulk UMI and Ag specificity score
### Find out the BEAM-T positive cells 
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01_First_4/post_processing/genotype_demux/Exp01_geno/outs/per_sample_outs/Exp01_geno/antigen_analysis/"
Ag_spec <- read.csv(paste(maindir, "antigen_specificity_scores.csv", sep = ""))
barcodes <- unique(Ag_spec$barcode)

## Creating an empty dataframe for all Ag information
rm(df)
df <- data.frame(matrix(ncol = 7, nrow = length(barcodes)))
colnames(df) <- c("samplename","barcode","Ag_UMI","Con_UMI","Ag_Con_UMI","max_specificity_score","max_Ag")
df$barcode <- barcodes

#### Which barcode match the samples
Exp01_demux <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01_First_4/post_processing/genotype_demux/Exp01_cellbarcode_samplename.txt",
                          sep = ",", header = TRUE)

df$samplename <- Exp01_demux[match(df$barcode, Exp01_demux$Barcode),"Assignment"]

for (i in 1:length(barcodes)) {
  subset_table <- Ag_spec[grep(barcodes[i], Ag_spec$barcode),c("antigen","antigen_umi","control_umi","antigen_specificity_score")]
  df[i,"Ag_UMI"] <- sum(subset_table$antigen_umi)
  df[i,"Con_UMI"] <- sum(subset_table$control_umi)
  df[i,"Ag_Con_UMI"] <- sum(subset_table$antigen_umi,subset_table$control_umi)
  df[i,"max_specificity_score"] <- max(subset_table[,"antigen_specificity_score"])
  df[i,"max_Ag"] <- paste(subset_table[grep(df[i,"max_specificity_score"], subset_table[,"antigen_specificity_score"]),"antigen"], collapse = ", ")
}

#### We wanted to keep the range of less than 10 if the Antigen specificity score is greater than 50

rm(df)
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01_First_4/post_processing/genotype_demux/Exp01_geno/outs/per_sample_outs/Exp01_geno/antigen_analysis/"

Ag_spec <- read.csv(paste(maindir, "antigen_specificity_scores.csv", sep = ""))
barcodes <- unique(Ag_spec$barcode)

## Creating an empty dataframe for all Ag information
rm(df)
df <- data.frame(matrix(ncol = 7, nrow = length(barcodes)))
colnames(df) <- c("samplename","barcode","Ag_UMI","Con_UMI","Ag_Con_UMI","max_specificity_score","max_Ag")
df$barcode <- barcodes

#### Which barcode match the samples
Exp01_demux <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01_First_4/post_processing/genotype_demux/Exp01_cellbarcode_samplename.txt", 
                          sep = ",", header = TRUE)

df$samplename <- Exp01_demux[match(df$barcode, Exp01_demux$Barcode),"Assignment"]

for (i in 1:length(barcodes)) {
  subset_table <- Ag_spec[grep(barcodes[i], Ag_spec$barcode),c("antigen","antigen_umi","control_umi","antigen_specificity_score")]
  df[i,"Ag_UMI"] <- sum(subset_table$antigen_umi)
  df[i,"Con_UMI"] <- sum(subset_table$control_umi)
  df[i,"Ag_Con_UMI"] <- sum(subset_table$antigen_umi,subset_table$control_umi)
  max_specificity_score <- max(subset_table[,"antigen_specificity_score"])
  if(max_specificity_score >= 50){
    # print("score greater than 50")
    ### Finding the score within -10 range
    specificity_score_range <- max_specificity_score - 10
    subset_table_gt_score_range <- subset_table[subset_table[,"antigen_specificity_score"] > specificity_score_range,]
    subset_table_gt_score_range_ordered <- subset_table_gt_score_range[order(-subset_table_gt_score_range$antigen_specificity_score),]
    df[i,"max_specificity_score"] <- paste(subset_table_gt_score_range_ordered$antigen_specificity_score,collapse=", ")
    df[i,"max_Ag"] <- paste(subset_table_gt_score_range_ordered[,"antigen"],collapse=", ")
  }
  else{
    # print("score less than 50")
    df[i,"max_specificity_score"] <- max(subset_table[,"antigen_specificity_score"])
    df[i,"max_Ag"] <- paste(subset_table[grep(df[i,"max_specificity_score"], subset_table[,"antigen_specificity_score"]),"antigen"], collapse = ", ")
  } 
}

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Table/"
write.table(df, paste(savedir,"Antigen_score_Exp01_with_range_10.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

### Exp02
rm(df)
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02_First_4/post_processing/genotype_demux/Exp02_geno/outs/per_sample_outs/Exp02_geno/antigen_analysis/"

Ag_spec <- read.csv(paste(maindir, "antigen_specificity_scores.csv", sep = ""))
barcodes <- unique(Ag_spec$barcode)

## Creating an empty dataframe for all Ag information
rm(df)
df <- data.frame(matrix(ncol = 7, nrow = length(barcodes)))
colnames(df) <- c("samplename","barcode","Ag_UMI","Con_UMI","Ag_Con_UMI","max_specificity_score","max_Ag")
df$barcode <- barcodes

#### Which barcode match the samples
Exp02_demux <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02_First_4/post_processing/genotype_demux/Exp02_cellbarcode_samplename.txt", 
                          sep = ",", header = TRUE)

df$samplename <- Exp02_demux[match(df$barcode, Exp02_demux$Barcode),"Assignment"]

for (i in 1:length(barcodes)) {
  subset_table <- Ag_spec[grep(barcodes[i], Ag_spec$barcode),c("antigen","antigen_umi","control_umi","antigen_specificity_score")]
  df[i,"Ag_UMI"] <- sum(subset_table$antigen_umi)
  df[i,"Con_UMI"] <- sum(subset_table$control_umi)
  df[i,"Ag_Con_UMI"] <- sum(subset_table$antigen_umi,subset_table$control_umi)
  max_specificity_score <- max(subset_table[,"antigen_specificity_score"])
  if(max_specificity_score >= 50){
    # print("score greater than 50")
    ### Finding the score within -10 range
    specificity_score_range <- max_specificity_score - 10
    subset_table_gt_score_range <- subset_table[subset_table[,"antigen_specificity_score"] > specificity_score_range,]
    subset_table_gt_score_range_ordered <- subset_table_gt_score_range[order(-subset_table_gt_score_range$antigen_specificity_score),]
    df[i,"max_specificity_score"] <- paste(subset_table_gt_score_range_ordered$antigen_specificity_score,collapse=", ")
    df[i,"max_Ag"] <- paste(subset_table_gt_score_range_ordered[,"antigen"],collapse=", ")
  }
  else{
    # print("score less than 50")
    df[i,"max_specificity_score"] <- max(subset_table[,"antigen_specificity_score"])
    df[i,"max_Ag"] <- paste(subset_table[grep(df[i,"max_specificity_score"], subset_table[,"antigen_specificity_score"]),"antigen"], collapse = ", ")
  } 
}

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Table/"
write.table(df, paste(savedir,"Antigen_score_Exp02_with_range_10.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

### Exp03
rm(df)
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/genotype_demux/Exp03_geno/outs/per_sample_outs/Exp03_geno/antigen_analysis/"

Ag_spec <- read.csv(paste(maindir, "antigen_specificity_scores.csv", sep = ""))
barcodes <- unique(Ag_spec$barcode)

## Creating an empty dataframe for all Ag information
rm(df)
df <- data.frame(matrix(ncol = 7, nrow = length(barcodes)))
colnames(df) <- c("samplename","barcode","Ag_UMI","Con_UMI","Ag_Con_UMI","max_specificity_score","max_Ag")
df$barcode <- barcodes

#### Which barcode match the samples
Exp03_demux <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/genotype_demux/Exp03_cellbarcode_samplename.txt", 
                          sep = ",", header = TRUE)

df$samplename <- Exp03_demux[match(df$barcode, Exp03_demux$Barcode),"Assignment"]

for (i in 1:length(barcodes)) {
  subset_table <- Ag_spec[grep(barcodes[i], Ag_spec$barcode),c("antigen","antigen_umi","control_umi","antigen_specificity_score")]
  df[i,"Ag_UMI"] <- sum(subset_table$antigen_umi)
  df[i,"Con_UMI"] <- sum(subset_table$control_umi)
  df[i,"Ag_Con_UMI"] <- sum(subset_table$antigen_umi,subset_table$control_umi)
  max_specificity_score <- max(subset_table[,"antigen_specificity_score"])
  if(max_specificity_score >= 50){
    # print("score greater than 50")
    ### Finding the score within -10 range
    specificity_score_range <- max_specificity_score - 10
    subset_table_gt_score_range <- subset_table[subset_table[,"antigen_specificity_score"] > specificity_score_range,]
    subset_table_gt_score_range_ordered <- subset_table_gt_score_range[order(-subset_table_gt_score_range$antigen_specificity_score),]
    df[i,"max_specificity_score"] <- paste(subset_table_gt_score_range_ordered$antigen_specificity_score,collapse=", ")
    df[i,"max_Ag"] <- paste(subset_table_gt_score_range_ordered[,"antigen"],collapse=", ")
  }
  else{
    # print("score less than 50")
    df[i,"max_specificity_score"] <- max(subset_table[,"antigen_specificity_score"])
    df[i,"max_Ag"] <- paste(subset_table[grep(df[i,"max_specificity_score"], subset_table[,"antigen_specificity_score"]),"antigen"], collapse = ", ")
  }
}

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Table/"
write.table(df, paste(savedir,"Antigen_score_Exp03_with_range_10.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

### Making a matrix for Exp01 02 and 03
### Extract out the cells for each run which we need to add the columns
Exp01 <- read.csv("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01_First_4/post_processing/genotype_demux/Exp01_geno/outs/per_sample_outs/Exp01_geno/antigen_analysis/Exp01_antigen_specificity_scores.csv", header = FALSE)
Exp02 <- read.csv("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02_First_4/post_processing/genotype_demux/Exp02_geno/outs/per_sample_outs/Exp02_geno/antigen_analysis/Exp02_antigen_specificity_scores.csv", header = FALSE)
Exp03 <- read.csv("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/genotype_demux/Exp03_geno/outs/per_sample_outs/Exp03_geno/antigen_analysis/Exp03_antigen_specificity_scores.csv", header = FALSE)

Exp01$V1 <- paste("Exp01_geno_",Exp01$V1,sep="")
Exp02$V1 <- paste("Exp02_geno_",Exp02$V1,sep="")
Exp03$V1 <- paste("Exp03_geno_",Exp03$V1,sep="")

all_exp <- rbind(Exp01, Exp02, Exp03)
colnames(all_exp) <- c("barcode","antigen","antigen_umi","control","control_umi","antigen_specificity_score","mhc_allele","raw_clonotype_id","exact_subclonotype_id")

### cellids for the CD8 object so we can merge it afterwards
cellid <- rownames(CD8@meta.data)
i=1
rm(df)
df <- data.frame(matrix(ncol = 11, nrow = length(cellid)))
colnames(df) <- c("EBV_BALF4","EBV_BMLF1","EBV_BMRF1","EBV_BRLF1","EBV_EBNA1","EBV_EBNA3C","EBV_LMP1","EBV_LMP2","VZV_IE62","VZV_IE63","cellid")
df$cellid <- cellid

for (i in 1:length(cellid)) {
  # print(i)
  cellid_exp <- all_exp[grep(cellid[i], all_exp$barcode),]
  for (j in 1:nrow(cellid_exp)) {
    Ag_subset <- cellid_exp[j,]
    df[i,j] <- paste(Ag_subset[c(5,3,6)],collapse="-")
  }
}

df$max_Ag <- CD8@meta.data$max_Ag
df$Ag_UMI <- CD8@meta.data$Ag_UMI
df$Con_UMI <- CD8@meta.data$Con_UMI
df$celltypes <- CD8@meta.data$celltypes
df$individuals <- CD8@meta.data$orig.ident
df$individuals <- CD8@meta.data$CDR3
df$max_specificity_score <- CD8@meta.data$max_specificity_score

write.table(df, "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Table/all_con_Ag_UMI_specificity_score.txt",
            quote=F, row.names = F, col.names = T, sep = "\t")

### We are combining all the tables which include the max specificity score with the range of +/- 10
### Ag with all the detail
Ag_UMI <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Table/all_con_Ag_UMI_specificity_score.txt",
                     header = TRUE, sep = "\t")

### with max specificity svcore + / - 10
max_specificity_score <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Table/Antigen_score_all_exp_with_range_10.txt", 
           header = TRUE, sep = "\t")

max_specificity_score$Run <- gsub(".*._","",max_specificity_score$samplename)
max_specificity_score$Run_barcode <- paste(max_specificity_score$Run,"geno" ,max_specificity_score$barcode, sep= "_")
Ag_UMI$specificity_score_range_10 <- max_specificity_score[match( Ag_UMI$cellid, max_specificity_score$Run_barcode),"max_specificity_score"]
Ag_UMI$Ag_range_10 <- max_specificity_score[match( Ag_UMI$cellid, max_specificity_score$Run_barcode),"max_Ag"]

stopifnot(all(Ag_UMI$cellid ==  rownames((CD8@meta.data))))
Ag_UMI$cdr_nt <- CD8@meta.data$cdr_nt
Ag_UMI$cdr_aa <- CD8@meta.data$cdr_aa

write.table(Ag_UMI, paste(savedir,"Table/all_con_Ag_UMI_specificity_score_with_clone_nt_aa_separate.txt",sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)

CD8@meta.data$EBV_BALF4 <- Ag_UMI$EBV_BALF4
CD8@meta.data$EBV_BMLF1 <- Ag_UMI$EBV_BMLF1
CD8@meta.data$EBV_BMRF1 <- Ag_UMI$EBV_BMRF1
CD8@meta.data$EBV_BRLF1 <- Ag_UMI$EBV_BRLF1
CD8@meta.data$EBV_EBNA1 <- Ag_UMI$EBV_EBNA1
CD8@meta.data$EBV_EBNA3C <- Ag_UMI$EBV_EBNA3C
CD8@meta.data$EBV_LMP2 <- Ag_UMI$EBV_LMP2
CD8@meta.data$VZV_IE62 <- Ag_UMI$VZV_IE62
CD8@meta.data$VZV_IE63 <- Ag_UMI$VZV_IE63
CD8@meta.data$Ag_range_10 <- Ag_UMI$Ag_range_10
CD8@meta.data$specificity_score_range_10 <- Ag_UMI$specificity_score_range_10

ID_with_cdrnt <- table(paste(CD8$cdr_nt,CD8$orig.ident,sep="_")) %>% as.data.frame()
write.table(ID_with_cdrnt, paste(savedir, "Table/ID_with_cdrnt_table.txt",sep = ""), sep = "\t", row.names = F, col.names = T, quote = F)

saveRDS(CD8, paste(savedir, "saveRDS_obj/CD8_modality_integrate_cluster_0.6.RDS",sep = ""))

### Making a table for TCR including all the TRA, TRB colummn to it.
## It turns out that the TRA has maxium 2 as well as TRB also has 2
Ag_spec <- read.csv("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/Ag_specific_filtered_contig_annotation_wid_header.csv")

## Simplifying the table
Ag_simple <- as.data.frame(matrix(nrow = length(rownames(Ag_spec))))
colnames(Ag_simple) <- "barcode"
Ag_simple$barcode <- Ag_spec$barcode
Ag_simple$chain <- Ag_spec$chain
Ag_simple$length <- Ag_spec$length
Ag_simple$vdjc <- paste(Ag_spec$v_gene, Ag_spec$d_gene, Ag_spec$j_gene, Ag_spec$c_gene, sep="_")
Ag_simple$cdr1_cdr2_cdr3 <- paste(Ag_spec$cdr1, Ag_spec$cdr2, Ag_spec$cdr3,sep="_")
Ag_simple$cdr1_cdr2_cdr3_nt <- paste(Ag_spec$cdr1_nt, Ag_spec$cdr2_nt, Ag_spec$cdr3_nt,sep="_")
Ag_simple$fwr1_fwr2_fwr3 <- paste(Ag_spec$fwr1, Ag_spec$fwr2, Ag_spec$fwr3,sep="_")
Ag_simple$fwr1_fwr2_fwr3_nt <- paste(Ag_spec$fwr1_nt, Ag_spec$fwr2_nt, Ag_spec$fwr3_nt,sep="_")
Ag_simple$reads_umis <- paste( Ag_spec$reads, Ag_spec$umis,sep="_")

## For the table columns and rows convering all the TCRs into one for both the TCRs
## Since maximum TCRA and TCRB could be 2 so based on that we will make it
same_col <- c("vdjc","cdr_aa","cdr_nt","fwr_aa","fwr_nt","reads_umis")
columns <- c("barcodes",paste("TRA1",same_col,sep = "_"), paste("TRA2",same_col,sep = "_"), paste("TRB1",same_col,sep = "_"), paste("TRB2",same_col,sep = "_"))
rows <- unique(Ag_simple$barcode)
Ag_converge <- as.data.frame(matrix(nrow = length(rows), ncol = length(columns)))
colnames(Ag_converge) <- columns
Ag_converge$barcodes <- rows

for (i in 1:length(rows)) {
  barcode_extract <- Ag_simple[grep(paste("^",Ag_converge$barcodes[i],"$",sep=""),Ag_simple$barcode),]
  TRA_TRB_table <- as.data.frame(table(barcode_extract$chain))
  my_vector <- TRA_TRB_table$Var1
  
  if (any(grepl("TRA", my_vector))) {
    ### Since it has both the TRA and TRB it could have 1/2 TRA or 1/2 TRB so we will use the for loop in here
    ### TRA
    if (TRA_TRB_table$Freq[1] == 1) {
      barcode_extract_TRA <- barcode_extract[grep("TRA",barcode_extract$chain),]
      Ag_converge[i,"TRA1_vdjc"] <- barcode_extract_TRA[1,"vdjc"]
      Ag_converge[i,"TRA1_cdr_aa"] <- barcode_extract_TRA[1,"cdr1_cdr2_cdr3"]
      Ag_converge[i,"TRA1_cdr_nt"] <- barcode_extract_TRA[1,"cdr1_cdr2_cdr3_nt"]
      Ag_converge[i,"TRA1_fwr_aa"] <- barcode_extract_TRA[1,"fwr1_fwr2_fwr3"]
      Ag_converge[i,"TRA1_fwr_nt"] <- barcode_extract_TRA[1,"fwr1_fwr2_fwr3_nt"]
      Ag_converge[i,"TRA1_reads_umis"] <- barcode_extract_TRA[1,"reads_umis"]
      
    } else if (TRA_TRB_table$Freq[1] == 2){
      barcode_extract_TRA <- barcode_extract[grep("TRA",barcode_extract$chain),]
      Ag_converge[i,"TRA1_vdjc"] <- barcode_extract_TRA[1,"vdjc"]
      Ag_converge[i,"TRA2_vdjc"] <- barcode_extract_TRA[2,"vdjc"]
      
      ### cdr amino acids
      Ag_converge[i,"TRA1_cdr_aa"] <- barcode_extract_TRA[1,"cdr1_cdr2_cdr3"]
      Ag_converge[i,"TRA2_cdr_aa"] <- barcode_extract_TRA[2,"cdr1_cdr2_cdr3"]
      
      ### cdr nucleotide
      Ag_converge[i,"TRA1_cdr_nt"] <- barcode_extract_TRA[1,"cdr1_cdr2_cdr3_nt"]
      Ag_converge[i,"TRA2_cdr_nt"] <- barcode_extract_TRA[2,"cdr1_cdr2_cdr3_nt"]
      
      ### fwr aa
      Ag_converge[i,"TRA1_fwr_aa"] <- barcode_extract_TRA[1,"fwr1_fwr2_fwr3"]
      Ag_converge[i,"TRA2_fwr_aa"] <- barcode_extract_TRA[2,"fwr1_fwr2_fwr3"]
      
      ### fwr nucleotide
      Ag_converge[i,"TRA1_fwr_nt"] <- barcode_extract_TRA[1,"fwr1_fwr2_fwr3_nt"]
      Ag_converge[i,"TRA2_fwr_nt"] <- barcode_extract_TRA[2,"fwr1_fwr2_fwr3_nt"]
      
      ### reads_umis
      Ag_converge[i,"TRA1_reads_umis"] <- barcode_extract_TRA[1,"reads_umis"]
      Ag_converge[i,"TRA2_reads_umis"] <- barcode_extract_TRA[2,"reads_umis"]
    }
  }
  if (any(grepl("TRB", my_vector))) {
    if (TRA_TRB_table$Freq[1] == 1) {
      barcode_extract_TRB <- barcode_extract[grep("TRB",barcode_extract$chain),]
      Ag_converge[i,"TRB1_vdjc"] <- barcode_extract_TRB[1,"vdjc"]
      Ag_converge[i,"TRB1_cdr_aa"] <- barcode_extract_TRB[1,"cdr1_cdr2_cdr3"]
      Ag_converge[i,"TRB1_cdr_nt"] <- barcode_extract_TRB[1,"cdr1_cdr2_cdr3_nt"]
      Ag_converge[i,"TRB1_fwr_aa"] <- barcode_extract_TRB[1,"fwr1_fwr2_fwr3"]
      Ag_converge[i,"TRB1_fwr_nt"] <- barcode_extract_TRB[1,"fwr1_fwr2_fwr3_nt"]
      Ag_converge[i,"TRB1_reads_umis"] <- barcode_extract_TRB[1,"reads_umis"]
      
    } else if (TRA_TRB_table$Freq[1] == 2){
      barcode_extract_TRB <- barcode_extract[grep("TRB",barcode_extract$chain),]
      Ag_converge[i,"TRB1_vdjc"] <- barcode_extract_TRB[1,"vdjc"]
      Ag_converge[i,"TRB2_vdjc"] <- barcode_extract_TRB[2,"vdjc"]
      
      ### cdr amino acids
      Ag_converge[i,"TRB1_cdr_aa"] <- barcode_extract_TRB[1,"cdr1_cdr2_cdr3"]
      Ag_converge[i,"TRB2_cdr_aa"] <- barcode_extract_TRB[2,"cdr1_cdr2_cdr3"]
      
      ### cdr nucleotide
      Ag_converge[i,"TRB1_cdr_nt"] <- barcode_extract_TRB[1,"cdr1_cdr2_cdr3_nt"]
      Ag_converge[i,"TRB2_cdr_nt"] <- barcode_extract_TRB[2,"cdr1_cdr2_cdr3_nt"]
      
      ### fwr aa
      Ag_converge[i,"TRB1_fwr_aa"] <- barcode_extract_TRB[1,"fwr1_fwr2_fwr3"]
      Ag_converge[i,"TRB2_fwr_aa"] <- barcode_extract_TRB[2,"fwr1_fwr2_fwr3"]
      
      ### fwr nucleotide
      Ag_converge[i,"TRB1_fwr_nt"] <- barcode_extract_TRB[1,"fwr1_fwr2_fwr3_nt"]
      Ag_converge[i,"TRB2_fwr_nt"] <- barcode_extract_TRB[2,"fwr1_fwr2_fwr3_nt"]
      
      ### reads_umis
      Ag_converge[i,"TRB1_reads_umis"] <- barcode_extract_TRB[1,"reads_umis"]
      Ag_converge[i,"TRB2_reads_umis"] <- barcode_extract_TRB[2,"reads_umis"]
    } 
  }
}

remaining_cells <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/remaining_cells",header=FALSE)
same_col <- c("vdjc","cdr_aa","cdr_nt","fwr_aa","fwr_nt","reads_umis")
columns <- c("barcodes",paste("TRA1",same_col,sep = "_"), paste("TRA2",same_col,sep = "_"), paste("TRB1",same_col,sep = "_"), paste("TRB2",same_col,sep = "_"))
rows <- unique(remaining_cells$V1)
remaining_cells_df <- as.data.frame(matrix(nrow = length(rows), ncol = length(columns)))
colnames(remaining_cells_df) <- columns
remaining_cells_df$barcodes <- rows

Ag_converge_all_cells <- rbind(Ag_converge, remaining_cells_df)

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/"
write.table(Ag_converge, paste(savedir,"VDJ/TCR_converged.txt",sep = ""), sep = "\t", row.names = F, col.names = T, quote = F)

write.table(Ag_converge_all_cells, paste(savedir,"VDJ/TCR_converged_all.txt",sep = ""), sep = "\t", row.names = F, col.names = T, quote = F)

### Including everything CD8 seurat object
Ag_converge$TRA1_TRA2_TRB1_TRB2_cdrnt <- paste(Ag_converge$TRA1_cdr_nt,"-", Ag_converge$TRA2_cdr_nt,"--", Ag_converge$TRB1_cdr_nt,"-" ,Ag_converge$TRB2_cdr_nt,sep="")
Ag_converge$TRA1_TRA2_TRB1_TRB2_cdraa <- paste(Ag_converge$TRA1_cdr_aa,"-", Ag_converge$TRA2_cdr_aa,"--", Ag_converge$TRB1_cdr_aa,"-" ,Ag_converge$TRB2_cdr_aa,sep="")
Ag_converge$TRA1_TRA2_TRB1_TRB2_vdjc <- paste(Ag_converge$TRA1_vdjc,"-", Ag_converge$TRA2_vdjc,"--", Ag_converge$TRB1_vdjc,"-" ,Ag_converge$TRB2_vdjc,sep="")
Ag_converge$TRA1_TRA2_TRB1_TRB2_fwrnt <- paste(Ag_converge$TRA1_fwr_nt,"-", Ag_converge$TRA2_fwr_nt,"--", Ag_converge$TRB1_fwr_nt,"-" ,Ag_converge$TRB2_fwr_nt,sep="")
Ag_converge$TRA1_TRA2_TRB1_TRB2_fwr_aa <- paste(Ag_converge$TRA1_fwr_aa,"-", Ag_converge$TRA2_fwr_aa,"--", Ag_converge$TRB1_fwr_aa,"-" ,Ag_converge$TRB2_fwr_aa,sep="")
Ag_converge$TRA1_TRA2_TRB1_TRB2_reads_umis <- paste(Ag_converge$TRA1_reads_umis,"-", Ag_converge$TRA2_reads_umis,"--", Ag_converge$TRB1_reads_umis,"-" ,Ag_converge$TRB2_reads_umis,sep="")

CD8@meta.data$TRA1_TRA2_TRB1_TRB2_cdrnt <- Ag_converge[match(rownames(CD8@meta.data),Ag_converge$barcodes),"TRA1_TRA2_TRB1_TRB2_cdrnt"]
CD8@meta.data$TRA1_TRA2_TRB1_TRB2_cdraa <- Ag_converge[match(rownames(CD8@meta.data),Ag_converge$barcodes),"TRA1_TRA2_TRB1_TRB2_cdraa"]
CD8@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc <- Ag_converge[match(rownames(CD8@meta.data),Ag_converge$barcodes),"TRA1_TRA2_TRB1_TRB2_vdjc"]
CD8@meta.data$TRA1_TRA2_TRB1_TRB2_fwrnt <- Ag_converge[match(rownames(CD8@meta.data),Ag_converge$barcodes),"TRA1_TRA2_TRB1_TRB2_fwrnt"]
CD8@meta.data$TRA1_TRA2_TRB1_TRB2_fwr_aa <- Ag_converge[match(rownames(CD8@meta.data),Ag_converge$barcodes),"TRA1_TRA2_TRB1_TRB2_fwr_aa"]
CD8@meta.data$TRA1_TRA2_TRB1_TRB2_reads_umis <- Ag_converge[match(rownames(CD8@meta.data),Ag_converge$barcodes),"TRA1_TRA2_TRB1_TRB2_reads_umis"]


clonal_expansion <- table(paste(CD8@meta.data$TRA1_TRA2_TRB1_TRB2_cdrnt, CD8@meta.data$orig.ident,sep="_")) %>% as.data.frame()
colnames(clonal_expansion) <- c("TRA1_cdr1_cdr2_chr3-TRA2_cdr1_cdr2_chr3--TRB1_cdr1_cdr2_chr3-TRB2_cdr1_cdr2_chr3_individual", "clones")
write.table(clonal_expansion, paste(savedir,"VDJ/clonal_expanion.txt",sep = ""), sep = "\t", row.names = F, col.names = T, quote = F)

TCR_ind_nt_Ag <- table(paste(CD8@meta.data$TRA1_TRA2_TRB1_TRB2_cdrnt,CD8@meta.data$orig.ident,sep="_"),CD8@meta.data$Ag_range_10)
write.table(TCR_ind_nt_Ag, paste(savedir,"VDJ/TCR_ind_nt_Ag.txt",sep = ""), sep = "\t", row.names = T, col.names = T, quote = F)


#### Identifying the standard deviation for the different Ag to decide the cutoff for the inclusion criteria
### After the discussion we decided to take the TRA1 TRA2 and TRB1 for the cdr nucleotide
CD8_TCR_ind_ID <- dplyr::select(CD8@meta.data, TRA1_TRA2_TRB1_TRB2_cdrnt,orig.ident)
CD8_TCR_ind_ID$TRA <- gsub("--.*.","",CD8_TCR_ind_ID$TRA1_TRA2_TRB1_TRB2_cdrnt)
CD8_TCR_ind_ID$TRA[is.na(CD8_TCR_ind_ID$TRA)] <- "Notfound"
CD8_TCR_ind_ID$TRB1 <- gsub(".*.--","",CD8_TCR_ind_ID$TRA1_TRA2_TRB1_TRB2_cdrnt) %>% gsub("-.*.","",.)
CD8_TCR_ind_ID$TRB1[is.na(CD8_TCR_ind_ID$TRB1)] <- "Notfound"
CD8_TCR_ind_ID$TRA_TRB1_IndID <- paste(CD8_TCR_ind_ID$TRA, CD8_TCR_ind_ID$TRB1, CD8_TCR_ind_ID$orig.ident, sep="-")
CD8_TCR_ind_ID_table <- table(CD8_TCR_ind_ID$TRA_TRB1_IndID) %>% as.data.frame()
CD8_TCR_ind_ID_table_ordered <- CD8_TCR_ind_ID_table[order(-CD8_TCR_ind_ID_table$Freq),]
write.table(CD8_TCR_ind_ID_table_ordered, paste(savedir,"VDJ/TCRA1_TRA2_TRB1_indID_clonal_expansion.txt",sep = ""), sep = "\t", row.names = F, col.names = T, quote = F)
CD8_TCR_ind_ID_table_ordered_gt_50 <- CD8_TCR_ind_ID_table_ordered[CD8_TCR_ind_ID_table_ordered$Freq > 50,]
CD8_TCR_ind_ID_table_ordered_gt_50_req <- CD8_TCR_ind_ID_table_ordered_gt_50[grep("Notfound",CD8_TCR_ind_ID_table_ordered_gt_50$Var1,invert=TRUE),]
CD8_TCR_ind_ID_table_ordered_gt_50_req$Var1 <- as.character(CD8_TCR_ind_ID_table_ordered_gt_50_req$Var1)

stopifnot(all(rownames(CD8_TCR_ind_ID) == rownames(CD8@meta.data)))
CD8@meta.data$TRA_TRB1_IndID <- CD8_TCR_ind_ID$TRA_TRB1_IndID

### Extracting out all the Ags and the clones
Ag_score <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Table/all_con_Ag_UMI_specificity_score_with_clone_nt_aa_separate.txt", header = TRUE, sep = "\t")
stopifnot(all(rownames(CD8@meta.data) == Ag_score$cellid))
# forgot to add LMP1
CD8@meta.data$EBV_LMP1 <- Ag_score$EBV_LMP1
Ag_TRA_TRB1 <- select(CD8@meta.data, EBV_BALF4, EBV_BMLF1, EBV_BMRF1, EBV_BRLF1, EBV_EBNA1, EBV_EBNA3C, EBV_LMP1, EBV_LMP2, VZV_IE62, VZV_IE63, TRA_TRB1_IndID)

### Extracting out the Clones greater than 50
Ag_TRA_TRB1_req <- Ag_TRA_TRB1[grep(paste(paste(CD8_TCR_ind_ID_table_ordered_gt_50_req$Var1, collapse="|")), Ag_TRA_TRB1$TRA_TRB1_IndID),]

### Taking the specificity score
for (i in 1:10) {
  Ag_TRA_TRB1_req[,i] <- gsub(".*.-","",Ag_TRA_TRB1_req[,i])
}

### Making teh gaussian plot for each of the clones
clones <- CD8_TCR_ind_ID_table_ordered_gt_50_req$Var1
Ag <- colnames(Ag_TRA_TRB1_req)[1:10]
 
# guassian_plot <- function(x,y, mean_val, std_dev){
#   plot(x,y, type = "l", lwd = 2)
#   points(x, y, col = "blue", pch = 19)
#   abline(v = mean_val, col = "red", lty = 2, lwd = 2)  # Mean
#   abline(v = mean_val + std_dev, col = "blue", lty = 2, lwd = 2)  # Mean + 1 std deviation
#   abline(v = mean_val - std_dev, col = "blue", lty = 2, lwd = 2)  # Mean - 1 std deviation
#   abline(v = mean_val + 2 * std_dev, col = "green", lty = 2, lwd = 2)  # Mean + 2 std deviations
#   abline(v = mean_val - 2 * std_dev, col = "green", lty = 2, lwd = 2)  # Mean - 2 std deviations
# }

rm(df_out)
df_out = as.data.frame(matrix(nrow = length(clones), ncol =length(Ag)))
colnames(df_out) <- Ag
rownames(df_out) <- clones

require(gridExtra)
for (i in 1:length(clones)) {
  rm(plot_list)
  plot_list <- list()
  Ag_TRA_TRB1_req_subset <- Ag_TRA_TRB1_req[grep(clones[i],Ag_TRA_TRB1_req$TRA_TRB1_IndID),]
  ## We have to make it for each antigen gaussian, the speecificity score needs to be greater than 50
  for (j in 1:length(Ag)) {
    specificity_score <- as.numeric(Ag_TRA_TRB1_req_subset[,Ag[j]])
    ss_gt_50 <- specificity_score[specificity_score > 50]
    if (length(ss_gt_50) > 0) {
      # Calculate mean and standard deviation
      mean_val <- mean(ss_gt_50)
      std_dev <- sd(ss_gt_50)
      
      x=ss_gt_50[order(ss_gt_50)]
      # Create a sequence of y values for the Gaussian curve
      y <- dnorm(x, mean = mean_val, sd = std_dev)
      
      df <- data.frame(SS = x, density=y)
      plot_list[[j]] <- ggplot(df, aes(SS, density)) + 
        geom_point(color = "mediumpurple3", size = 2) + 
        geom_smooth(se = FALSE, color = "black") +
        geom_vline(xintercept = mean_val, linetype="solid", color="red", size = 1) +
        geom_vline(xintercept = mean_val + std_dev, linetype="solid", color="blue", size = 1) +
        geom_vline(xintercept = mean_val - std_dev, linetype="solid", color="blue", size = 1) +
        geom_vline(xintercept = mean_val + 2*std_dev, linetype="solid", color="magenta", size = 1) +
        geom_vline(xintercept = mean_val - 2*std_dev, linetype="solid", color="magenta", size = 1) + 
        scale_x_continuous(limits = c(50,100)) +
        ggtitle(paste(Ag[j],"clones=",nrow(df),"u=",round(mean_val,1),"std=",round(std_dev,1))) + 
        theme(plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.grid.minor = element_line(colour = "white"),
              panel.grid.major = element_line(colour = "white"))
      df_out[i,j] <- paste(nrow(df),mean_val,std_dev,sep = "-")
    }
    else{
      df2 <- data.frame()
      plot_list[[j]] <- ggplot(df2) + geom_point() + xlim(50, 100) + ylim(0, 1) +
        ggtitle(Ag[j]) + 
        theme(plot.title = element_text(hjust = 0.5))
    }
  }
  dir.create(paste(savedir,"Antigen/Table/gaussian_plot/",sep = ""), showWarnings = FALSE)
  pdf(paste(savedir, "Antigen/Table/gaussian_plot/clone_",i,".pdf",sep = ""), width = 16, height = 14)
  grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]],
               plot_list[[6]], plot_list[[7]], plot_list[[8]], plot_list[[9]], plot_list[[10]],
               nrow=3, ncol=4)
  dev.off()
}

write.table(df_out, paste(savedir, "Antigen/Table/gaussian_plot/clones_mean_std_dev_each_Ag.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

Ag_scores <- dplyr::select(CD8_metadata, EBV_BALF4, EBV_BMLF1, EBV_BMRF1, EBV_BRLF1, EBV_EBNA1, EBV_EBNA3C, EBV_LMP1, EBV_LMP2, VZV_IE62,VZV_IE63, Ag_range_10)

write.table(Ag_scores, paste())








