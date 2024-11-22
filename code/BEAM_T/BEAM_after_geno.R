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
## We have three experiments starting 
Exp02 <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Table/Exp02_Ag_UMI_gt_50.txt", header = TRUE, sep = "\t")
Exp01 <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Table/Exp01_Ag_UMI_gt_50.txt", header = TRUE, sep = "\t")
Exp03 <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Table/Exp03_Ag_UMI_gt_50.txt", header = TRUE, sep = "\t")

Exp01_02_03_Ag <- rbind(Exp01, Exp02, Exp03)
write.table(Exp01_02_03_Ag, "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Table/Exp01_02_03_Ag_cutoff.txt",
            sep = "\t", row.names = F, col.names = T, quote = F)

# Ag_number <- table(Exp01_02_03_Ag$samplename, Exp01_02_03_Ag$max_Ag)
# Ag_specific <- Ag_number[,grep(",",colnames(Ag_number), invert=TRUE)] %>% as.matrix()
# Ag_spec_df <- as.data.frame(Ag_specific)
# 
# split_to_dataframe <- function(x) {
#   split_elements <- strsplit(x, "_")[[1]]
#   data.frame(t(split_elements))
# }
# Ag_spec_df$Var1 <- as.character(Ag_spec_df$Var1)
# 
# # Convert split elements to dataframe
# Ag_spec <- do.call(rbind, lapply(Ag_spec_df$Var1, split_to_dataframe))
# colnames(Ag_spec) <- c("SampleID","Virus","Age","Age_number","Gender","Batch")
# Ag_spec$orig.ident <- Ag_spec_df$Var1
# Ag_spec$Ag <- Ag_spec_df$Var2
# Ag_spec$cell_number <- Ag_spec_df$Freq
# 
# Ag_spec %>% group_by(Ag,Age) %>% summarize(total_cell_number = sum(cell_number)) %>% as.data.frame() -> Ag_Age_number
# Ag_Age_number_casted <- dcast(Ag_Age_number,Ag~Age)
# savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/"
# write.table(Ag_Age_number_casted, paste(savedir,"Table/Ag_Age.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
# 
# ### Stats and Cell number 
# library(ggpubr)
# library(rstatix)
# 
# stat.test <- Ag_spec %>%
#   group_by(Ag) %>%
#   t_test(cell_number ~ Age) %>%
#   adjust_pvalue(method = "holm") %>%
#   add_significance("p.adj")
# stat.test
# 
# bxp <- ggboxplot(
#   Ag_spec, x = "Ag", y = "cell_number", color = "Age", palette = c("#00AFBB", "#E7B800")
# ) 
# 
# bxp2 <- bxp +   geom_dotplot(
#   aes(fill = Age, color = Age), trim = FALSE,
#   binaxis='y', stackdir='center', dotsize = 0.5,
#   position = position_dodge(0.8)
# )+
#   scale_fill_manual(values = c("#00AFBB", "#E7B800"))+
#   scale_color_manual(values = c("#00AFBB", "#E7B800"))
# 
# 
# # Add p-values onto the box plots
# stat.test <- stat.test %>%
#   add_xy_position(x = "Ag", dodge = 0.8)
# bxp3 <- bxp2 + stat_pvalue_manual(
#   stat.test,  label = "p", tip.length = 0
# ) +   theme(plot.title = element_text(hjust = 0.5),
#             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
#             panel.background = element_rect(fill = 'white', colour = 'white'),
#             panel.grid.minor = element_line(colour = "white"),
#             panel.grid.major = element_line(colour = "white"))
# 
# # Add 10% spaces between the p-value labels and the plot border
# 
# bxp4 <- bxp3 + stat_pvalue_manual(
#   stat.test,  label = "p", tip.length = 0
# ) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
# 
# 
# pdf(paste(savedir,"Antigen/Table/Ag_O_vs_Y_bxplot.pdf",sep = ""), width = 8, height = 6)
# bxp4
# dev.off()

### Preprocessing All ####
library("Matrix")
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/CR7_scCITESeq_QC.R")

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/all/"
Exp01 = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01_First_4/post_processing/genotype_demux/Exp01_geno/outs/per_sample_outs/Exp01_geno/"
Exp02 = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02_First_4/post_processing/genotype_demux/Exp02_geno/outs/per_sample_outs/Exp02_geno/"
Exp03 = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/genotype_demux/Exp03_geno/outs/per_sample_outs/Exp03_geno/"

samplepath <- c(Exp01,Exp02,Exp03)
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
combined <- merge(x=readRDS(RNA_RDS[1]), y = c(readRDS(RNA_RDS[2]),readRDS(RNA_RDS[3])), add.cell.ids = samplename, project = "combined")

Exp01_demux <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01_First_4/post_processing/genotype_demux/Exp01_cellbarcode_samplename.txt",
                          sep = ",", header = TRUE)
Exp02_demux <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp02_First_4/post_processing/genotype_demux/Exp02_cellbarcode_samplename.txt", 
                          sep = ",", header = TRUE)
Exp03_demux <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp03_First_4/post_processing/genotype_demux/Exp03_cellbarcode_samplename.txt", 
                          sep = ",", header = TRUE)

Exp01_demux$Barcode <- paste("Exp01_geno_",Exp01_demux$Barcode,sep="")
Exp02_demux$Barcode <- paste("Exp02_geno_",Exp02_demux$Barcode,sep="")
Exp03_demux$Barcode <- paste("Exp03_geno_",Exp03_demux$Barcode,sep="")

Exp1_2_3 <- rbind(Exp01_demux, Exp02_demux, Exp03_demux)
cellname_remove <- combined@meta.data[which(is.na(combined@meta.data$samplename)),] %>% rownames()
`%notin%` <- Negate(`%in%`)
cell_req <- rownames(combined@meta.data)[which(rownames(combined@meta.data) %notin% cellname_remove)]
combined_demux <- subset(combined, cells= cell_req)

obj_ADT <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "_ADT_singlet.RDS",full.names = TRUE)
ADT_RDS <- c(obj_ADT)

# samplename_ADT <- paste(samplename,"ADT",sep = "_")
## We have to combine the RNA singlet object to perform the normalization for the unwanted the effects.
combined_ADT <- merge(x=readRDS(ADT_RDS[1]), y = c(readRDS(ADT_RDS[2]),readRDS(ADT_RDS[3])), add.cell.ids = samplename, project = "combined_ADT")
combined_ADT_demux <- subset(combined_ADT, cells= cell_req)

#### Further separating out the bulk and Ag
### Ag subset
Exp01_Ag <- read.table(paste(maindir,"Antigen/Table/Exp01_Ag_UMI_gt_50.txt",sep=""), header=TRUE, sep="\t")
Exp02_Ag <- read.table(paste(maindir,"Antigen/Table/Exp02_Ag_UMI_gt_50.txt",sep=""), header=TRUE, sep="\t")
Exp03_Ag <- read.table(paste(maindir,"Antigen/Table/Exp03_Ag_UMI_gt_50.txt",sep=""), header=TRUE, sep="\t")
Ag_cells <- rbind(Exp01_Ag, Exp02_Ag, Exp03_Ag)
Ag_cells_noNA <- Ag_cells[!is.na(Ag_cells$samplename),]
Ag_cells_noNA$Exp <- paste(gsub(".*._","",Ag_cells_noNA$samplename),"_geno",sep="")
Ag_cells_noNA$Run_barcode <- paste(Ag_cells_noNA$Exp, Ag_cells_noNA$barcode,sep="_")
Ag_specific_cells <- Ag_cells_noNA$Run_barcode
write.table(Ag_specific_cells, paste(maindir,"Antigen/Table/Ag_specific_cells.txt",sep = ""), sep = "\t", row.names = F, col.names = F, quote = F)

Ag_obj <- subset(combined_demux, cells = Ag_specific_cells)
Ag_ADT_obj <- subset(combined_ADT_demux, cells = Ag_specific_cells)

dir.create(paste(maindir,"Antigen/saveRDS_obj",sep = ""), showWarnings = FALSE)
saveRDS(Ag_obj, paste(maindir,"Antigen/saveRDS_obj/Ag_obj.RDS",sep = ""))
saveRDS(Ag_ADT_obj, paste(maindir,"Antigen/saveRDS_obj/Ag_ADT_obj.RDS",sep = ""))

### Bulk Obj
bulk_cells <- read.table(paste(maindir,"Antigen/Table/Exp01_bulk_cells.txt",sep=""), header=FALSE, sep="\t")
colnames(bulk_cells) <- colnames(Ag_cells)
bulk_cells_noNA <- bulk_cells[!is.na(bulk_cells$samplename),]
bulk_cells_noNA$Exp <- paste(gsub(".*._","",bulk_cells_noNA$samplename),"_geno",sep="")
bulk_cells_noNA$Run_barcode <- paste(bulk_cells_noNA$Exp, bulk_cells_noNA$barcode,sep="_")
bulk_specific_cells <- bulk_cells_noNA$Run_barcode
dir.create(paste(maindir,"Bulk/Table",sep = ""), showWarnings = F)
write.table(bulk_specific_cells, paste(maindir,"Bulk/Table/bulk_specific_cells.txt",sep = ""), sep = "\t", row.names = F, col.names = F, quote = F)

bulk_obj <- subset(combined_demux, cells = bulk_specific_cells)
bulk_ADT_obj <- subset(combined_ADT_demux, cells = bulk_specific_cells)

dir.create(paste(maindir,"Bulk/saveRDS_obj",sep = ""), showWarnings = FALSE)
saveRDS(Ag_obj, paste(maindir,"Bulk/saveRDS_obj/bulk_obj.RDS",sep = ""))
saveRDS(Ag_ADT_obj, paste(maindir,"Bulk/saveRDS_obj/bulk_ADT_obj.RDS",sep = ""))

### Bulk ####
library("Matrix")
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.

maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/"
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/"
dir.create(savedir,showWarnings = FALSE)

# bulk_obj <- readRDS(paste(savedir, "saveRDS_obj/bulk_obj.RDS",sep = ""))
# bulk_ADT_obj <- readRDS(paste(savedir, "saveRDS_obj/bulk_ADT_obj.RDS",sep = ""))

### Since the TCR genes and the MTRNR2L8 are making various clusters that is not required so we will remove these genes from the object
remove_genes <- c("MTRNR2L8",grep("^TRAV|^TRAJ|^TRBC|^TRBV|^TRBJ|^TRBC",
                                  rownames(bulk_obj@assays$RNA@counts),
                                  value=TRUE))
keep_genes <- grep(paste(remove_genes,collapse = "|"),
                   rownames(bulk_obj@assays$RNA@counts), 
                   invert = TRUE, value = TRUE)
bulk_obj_no_TCR <- bulk_obj[keep_genes,]

### Just to check
grep("MTRNR2L8|^TRAV|^TRAJ|^TRBC|^TRBV|^TRBJ|^TRBC",
     rownames(bulk_obj_no_TCR@assays$RNA@counts),
     value=TRUE)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_subsetting.R")
RNA_obj_path = bulk_obj_no_TCR
ADT_obj_path = bulk_ADT_obj
ADT_main = bulk_ADT_obj
CD8_cluster <- NULL
bulk_cellnames <- rownames(bulk_ADT_obj@meta.data)
objname = "CD8"
Assay = "RNA"
process = "bulk_subset"

CD8_RNA <- RNA_subsetting(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path, ADT_main = ADT_main, subset_cluster = CD8_cluster,
                          cellnames=bulk_cellnames, ngenes = 4000, saveDir = savedir, objname = objname, Assay = Assay, process = process)

CD8_RNA <- RunUMAP(CD8_RNA, dims = 1:25, reduction = "pca")

CD8_RNA@meta.data$Run <- CD8_RNA@meta.data$orig.ident
CD8_RNA@meta.data$orig.ident <- CD8_RNA@meta.data$samplename

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
                                                    dims = 30,
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
CD8_integrated_RNA <- RNA_integration(CD8_RNA_sctransformed, savedir, dims = 30, RNA_features = c("CD4","CD8A"),
                                      Assay=Assay, process=process, objname=objname, ncol = 4, ndims = 50)

dir.create(paste(savedir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/CD8_integration_Run.pdf",sep = ""),width = 6, height = 5)
print(DimPlot(CD8_integrated_RNA, group.by = "Run"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_integration_Run_splitted.pdf",sep = ""),width = 15, height = 5)
print(DimPlot(CD8_integrated_RNA, split.by = "Run", group.by="Run", ncol = 3) + NoLegend())
dev.off()

CD8_integrated_RNA_NN <- FindNeighbors(CD8_integrated_RNA, dims = 1:30)

CD8_integrated_RNA_NN@meta.data$sample_id = gsub("_.*","",CD8_integrated_RNA_NN@meta.data$orig.ident)
CD8_integrated_RNA_NN@meta.data$Run = gsub(".*_","",CD8_integrated_RNA_NN@meta.data$orig.ident)
CD8_integrated_RNA_NN@meta.data$Age <- gsub("[0-9]","",CD8_integrated_RNA_NN@meta.data$sample_id)
CD8_integrated_RNA_NN@meta.data$gender <- gsub(".*._M_.*.","M",CD8_integrated_RNA_NN@meta.data$orig.ident) %>% gsub(".*._F_.*.","F",.)
CD8_integrated_RNA_NN@meta.data$Age_number <- gsub(".*._O_|.*._Y_|_M_.*.|_F_.*.","",CD8_integrated_RNA_NN@meta.data$orig.ident)

pdf(paste(savedir,"UMAP/CD8_RNA_Age.pdf",sep = ""))
DimPlot(CD8_integrated_RNA_NN, group.by = "Age")
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_batch.pdf",sep = ""))
DimPlot(CD8_integrated_RNA_NN, group.by = "Run")
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_Age_split.pdf",sep = ""), width = 10, height = 6)
DimPlot(CD8_integrated_RNA_NN, group.by = "Age", split.by = "Age", ncol = 2)
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_samples.pdf",sep = ""))
DimPlot(CD8_integrated_RNA_NN, group.by = "orig.ident")
dev.off()

pdf(paste(savedir,"UMAP/CD8_RNA_samples_split.pdf",sep = ""), width = 15, height = 12)
DimPlot(CD8_integrated_RNA_NN, group.by = "orig.ident", split.by = "orig.ident", ncol = 4)
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
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/"
CD8_integrated_RNA_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster.RDS",sep = ""))
DefaultAssay(CD8_integrated_RNA_NN_cluster) <- "RNA"
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/miniconda3/envs/seurat_new/bin/python3.8")
library(Rmagic)
py_discover_config("magic") # to check
CD8_integrated_RNA_NN_cluster_impute <- magic(CD8_integrated_RNA_NN_cluster, npca=30) ## imputing the RNA data as for RNA PCs are 20
saveRDS(CD8_integrated_RNA_NN_cluster_impute,paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster_impute.RDS",sep = ""))

### CD8 ADT #######
CD8_ADT = readRDS(paste(savedir,"saveRDS_obj/CD8_bulk_subset_ADT.RDS",sep = ""))
all(rownames(CD8_integrated_RNA_NN_cluster_impute@meta.data) == rownames(CD8_ADT@meta.data))

CD8_ADT@meta.data$Run <- CD8_ADT@meta.data$orig.ident
CD8_ADT@meta.data$orig.ident <- CD8_integrated_RNA_NN_cluster@meta.data$samplename

CD8_ADT@meta.data$sample_id = gsub("_.*","",CD8_ADT@meta.data$orig.ident)
CD8_ADT@meta.data$Run = gsub(".*_","",CD8_ADT@meta.data$orig.ident)
CD8_ADT@meta.data$Age <- gsub("[0-9]","",CD8_ADT@meta.data$sample_id)
CD8_ADT@meta.data$gender <- gsub(".*._M_.*.","M",CD8_ADT@meta.data$orig.ident) %>% gsub(".*._F_.*.","F",.)
CD8_ADT@meta.data$Age_number <- gsub(".*._O_|.*._Y_|_M_.*.|_F_.*.","",CD8_ADT@meta.data$orig.ident)

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

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")
objname = "CD8"
process = "integrated"
Assay = "ADT"

CD8_ADT_integrated <- ADT_merging(CD8_ADT, savedir, dims = 10, numfeatures=36,
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
# saveRDS(CD8_ADT_integrated_NN_cluster,paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster.RDS",sep = ""))

## Since Magic does not work on the conda seurat_new we have to perform the analysis outside conda.
## ml load r/4.2.2
library(Seurat)
# savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/Antigen/"
# CD8_ADT_integrated_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster.RDS",sep = ""))
DefaultAssay(CD8_ADT_integrated_NN_cluster) <- "ADT"
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/miniconda3/envs/seurat_new/bin/python3.8")
library(Rmagic)
py_discover_config("magic") # to check
CD8_ADT_integrated_NN_cluster_impute <- magic(CD8_ADT_integrated_NN_cluster, npca=10) ## imputing the RNA data as for RNA PCs are 20
saveRDS(CD8_ADT_integrated_NN_cluster_impute,paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster_impute.RDS",sep = ""))

### CD8 Modality Integration #####
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITEseq_modality_integration.R")
# Find the last save object for RNA and ADT
# we have RNA already in the path
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/"
## Since ADT we donot have to regress the TCR we can take it from the previous analysis
# CD8_ADT_integrated_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster_impute.RDS",sep = ""))
# CD8_integrated_RNA_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster_impute.RDS",sep = ""))
RNA_obj_path <- CD8_integrated_RNA_NN_cluster_impute
ADT_obj_path <- CD8_ADT_integrated_NN_cluster_impute

objname <- "CD8"
process <- "modality_integrate"
Assay <- "integrated"
CD8_modality_integrate <- modality_integration(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path,
                                               RNA_dims = 30, ADT_dims = 10,
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

# adding the celltypes to the object
clus_celltypes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/Table/clus_celltypes", header = TRUE)

### integrating the cluster 23 into cluster 1 and cluster 24 into cluster 0
CD8_modality_integrate_cluster@meta.data$seurat_clusters2 <- CD8_modality_integrate_cluster@meta.data$seurat_clusters
gsub("^23$","1",CD8_modality_integrate_cluster@meta.data$seurat_clusters2) %>% gsub("^24","0",.) -> CD8_modality_integrate_cluster@meta.data$seurat_clusters2
Idents(CD8_modality_integrate_cluster) <- CD8_modality_integrate_cluster@meta.data$seurat_clusters2

patterns = clus_celltypes$clus
replacements = clus_celltypes$celltypes
# names(replacement) <- patterns
CD8_modality_integrate_cluster@meta.data$celltypes <- CD8_modality_integrate_cluster@meta.data$seurat_clusters2
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  CD8_modality_integrate_cluster@meta.data$celltypes <- str_replace_all(CD8_modality_integrate_cluster@meta.data$celltypes, pattern, replacements[i])
}

pdf(paste(savedir,"UMAP/celltypes_wnnumap.pdf",sep = ""))
DimPlot(CD8_modality_integrate_cluster, reduction = "wnn.umap", group.by = "celltypes", label = TRUE, label.size = 3)
dev.off()

saveRDS(CD8_modality_integrate_cluster, paste(savedir, "saveRDS_obj/CD8_modality_integrate_cluster_0.6.RDS",sep = ""))

### Antigen ####
library("Matrix")
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.

maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/"
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/"
dir.create(savedir,showWarnings = FALSE)

# Ag_obj <- readRDS(paste(savedir, "saveRDS_obj/Ag_obj.RDS",sep = ""))
# Ag_ADT_obj <- readRDS(paste(savedir, "saveRDS_obj/Ag_ADT_obj.RDS",sep = ""))

### Since the TCR genes and the MTRNR2L8 are making various clusters that is not required so we will remove these genes from the object
remove_genes <- c("MTRNR2L8",grep("^TRAV|^TRAJ|^TRBC|^TRBV|^TRBJ|^TRBC",
                                  rownames(Ag_obj@assays$RNA@counts),
                                  value=TRUE))
keep_genes <- grep(paste(remove_genes,collapse = "|"),
                   rownames(Ag_obj@assays$RNA@counts), 
                   invert = TRUE, value = TRUE)
Ag_obj_no_TCR <- Ag_obj[keep_genes,]

### Just to check
grep("MTRNR2L8|^TRAV|^TRAJ|^TRBC|^TRBV|^TRBJ|^TRBC",
     rownames(Ag_obj_no_TCR@assays$RNA@counts),
     value=TRUE)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_subsetting.R")
RNA_obj_path = Ag_obj_no_TCR
ADT_obj_path = Ag_ADT_obj
ADT_main = Ag_ADT_obj
CD8_cluster <- NULL
Ag_cellnames <- rownames(Ag_ADT_obj@meta.data)
objname = "CD8"
Assay = "RNA"
process = "Ag_subset"

CD8_RNA <- RNA_subsetting(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path, ADT_main = ADT_main, subset_cluster = CD8_cluster,
                          cellnames=Ag_cellnames, ngenes = 4000, saveDir = savedir, objname = objname, Assay = Assay, process = process)

CD8_RNA <- RunUMAP(CD8_RNA, dims = 1:25, reduction = "pca")

CD8_RNA@meta.data$Run <- CD8_RNA@meta.data$orig.ident
CD8_RNA@meta.data$orig.ident <- CD8_RNA@meta.data$samplename

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
CD8_integrated_RNA <- RNA_integration(CD8_RNA_sctransformed, savedir, dims = 30, RNA_features = c("CD4","CD8A"),
                                      Assay=Assay, process=process, objname=objname, ncol = 4, ndims = 50)

dir.create(paste(savedir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/CD8_integration_Run.pdf",sep = ""),width = 6, height = 5)
print(DimPlot(CD8_integrated_RNA, group.by = "Run"))
dev.off()

pdf(paste(savedir,"UMAP/CD8_integration_Run_splitted.pdf",sep = ""),width = 15, height = 5)
print(DimPlot(CD8_integrated_RNA, split.by = "Run", group.by="Run", ncol = 3) + NoLegend())
dev.off()

CD8_integrated_RNA_NN <- FindNeighbors(CD8_integrated_RNA, dims = 1:20)


CD8_integrated_RNA_NN@meta.data$sample_id = gsub("_.*","",CD8_integrated_RNA_NN@meta.data$orig.ident)
CD8_integrated_RNA_NN@meta.data$Run = gsub(".*_","",CD8_integrated_RNA_NN@meta.data$orig.ident)
CD8_integrated_RNA_NN@meta.data$Age <- gsub("[0-9]","",CD8_integrated_RNA_NN@meta.data$sample_id)
CD8_integrated_RNA_NN@meta.data$gender <- gsub(".*._M_.*.","M",CD8_integrated_RNA_NN@meta.data$orig.ident) %>% gsub(".*._F_.*.","F",.)
CD8_integrated_RNA_NN@meta.data$Age_number <- gsub(".*._O_|.*._Y_|_M_.*.|_F_.*.","",CD8_integrated_RNA_NN@meta.data$orig.ident)

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
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/"
CD8_integrated_RNA_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster.RDS",sep = ""))
DefaultAssay(CD8_integrated_RNA_NN_cluster) <- "RNA"
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/miniconda3/envs/seurat_new/bin/python3.8")
library(Rmagic)
py_discover_config("magic") # to check
CD8_integrated_RNA_NN_cluster_impute <- magic(CD8_integrated_RNA_NN_cluster, npca=30) ## imputing the RNA data as for RNA PCs are 20
saveRDS(CD8_integrated_RNA_NN_cluster_impute,paste(savedir,"saveRDS_obj/CD8_integrated_RNA_NN_cluster_impute.RDS",sep = ""))

### CD8 ADT #######
CD8_ADT = readRDS(paste(savedir,"saveRDS_obj/CD8_Ag_subset_ADT.RDS",sep = ""))
all(rownames(CD8_integrated_RNA_NN_cluster@meta.data) == rownames(CD8_ADT@meta.data))

CD8_ADT@meta.data$Run <- CD8_ADT@meta.data$orig.ident
CD8_ADT@meta.data$orig.ident <- CD8_integrated_RNA_NN_cluster@meta.data$samplename

CD8_ADT@meta.data$sample_id = gsub("_.*","",CD8_ADT@meta.data$orig.ident)
CD8_ADT@meta.data$Run = gsub(".*_","",CD8_ADT@meta.data$orig.ident)
CD8_ADT@meta.data$Age <- gsub("[0-9]","",CD8_ADT@meta.data$sample_id)
CD8_ADT@meta.data$gender <- gsub(".*._M_.*.","M",CD8_ADT@meta.data$orig.ident) %>% gsub(".*._F_.*.","F",.)
CD8_ADT@meta.data$Age_number <- gsub(".*._O_|.*._Y_|_M_.*.|_F_.*.","",CD8_ADT@meta.data$orig.ident)

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

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")
objname = "CD8"
process = "integrated"
Assay = "ADT"

CD8_ADT_integrated <- ADT_merging(CD8_ADT, savedir, dims = 10, numfeatures=36,
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
# savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/"
# CD8_ADT_integrated_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster.RDS",sep = ""))
DefaultAssay(CD8_ADT_integrated_NN_cluster) <- "ADT"
library(reticulate)
use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/miniconda3/envs/seurat_new/bin/python3.8")
library(Rmagic)
py_discover_config("magic") # to check
CD8_ADT_integrated_NN_cluster_impute <- magic(CD8_ADT_integrated_NN_cluster, npca=10) ## imputing the RNA data as for RNA PCs are 20
saveRDS(CD8_ADT_integrated_NN_cluster_impute,paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster_impute.RDS",sep = ""))

# CD8_integrated_RNA_NN_cluster <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run06/Downstream/saveRDS_obj/CD8_integrated_RNA_NN_cluster.RDS")
### CD8 Modality Integration #####
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITEseq_modality_integration.R")
# Find the last save object for RNA and ADT
# we have RNA already in the path
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/"
## Since ADT we donot have to regress the TCR we can take it from the previous analysis
# CD8_ADT_integrated_NN_cluster <- readRDS(paste(savedir,"saveRDS_obj/CD8_ADT_integrated_NN_cluster_impute.RDS",sep = ""))
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

### Testing Leiden algorithm
CD8_modality_integrate_cluster <- FindClusters(CD8_modality_integrate_cluster, resolution = 0.6, 
                                               graph.name = "wsnn", algorithm = 4, verbose = TRUE)

pdf("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/UMAP/Ag_cluster_leiden.pdf")
DimPlot(CD8_modality_integrate_cluster, reduction = "wnn.umap")
dev.off()

# adding the celltypes to the object
clus_celltypes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Table/clus_celltypes", header = TRUE)

patterns = clus_celltypes$clus
replacements = clus_celltypes$celltypes
# names(replacement) <- patterns
CD8_modality_integrate_cluster@meta.data$celltypes <- CD8_modality_integrate_cluster@meta.data$seurat_clusters
for (i in seq_along(patterns)) {
  pattern <- paste0("\\b", patterns[i], "\\b")
  CD8_modality_integrate_cluster@meta.data$celltypes <- str_replace_all(CD8_modality_integrate_cluster@meta.data$celltypes, pattern, replacements[i])
}

pdf(paste(savedir,"UMAP/celltypes_wnnumap.pdf",sep = ""))
DimPlot(CD8_modality_integrate_cluster, reduction = "wnn.umap", group.by = "celltypes", label = TRUE, label.size = 3)
dev.off()

saveRDS(CD8_modality_integrate_cluster, paste(savedir, "saveRDS_obj/CD8_modality_integrate_cluster_0.6.RDS",sep = ""))

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/ADT_featureplot.R")
ADT_featureplot(obj = CD8_modality_integrate_cluster, savedir = savedir,objname = "CD8_modality_integrate_cluster",
                process = "featureplots",reduction = "wnn.umap",x = "wnnUMAP_1",y = "wnnUMAP_2")

CD8_samplewise_counts <- table(CD8_modality_integrate_cluster@meta.data$orig.ident) %>% as.data.frame()
write.table(CD8_samplewise_counts, paste(savedir,"Table/CD8_samplewise_counts.txt",sep = ""), sep = "\t", row.names = F, col.names = T, quote = F)

Exp01_02_03_Ag <- read.table( "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Table/Exp01_02_03_Ag_cutoff.txt",sep = "\t", header = TRUE)
Exp01_02_03_Ag$Run <- gsub(".*._","",Exp01_02_03_Ag$samplename)
Exp01_02_03_Ag$Exp_barcode <- paste(Exp01_02_03_Ag$Run, "_geno_",Exp01_02_03_Ag$barcode,sep = "")

CD8_modality_integrate_cluster@meta.data$Antigen <- Exp01_02_03_Ag[match(rownames(CD8_modality_integrate_cluster@meta.data), Exp01_02_03_Ag$Exp_barcode),"max_Ag"]
Ind_per_antigen <- table(CD8_modality_integrate@meta.data$orig.ident,CD8_modality_integrate_cluster@meta.data$Antigen)
colnames(Ind_per_antigen) <- gsub(", ","_",colnames(Ind_per_antigen))
write.table(Ind_per_antigen, paste(savedir,"Table/Ind_per_Ag.txt",sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)

CD8_modality_integrate@meta.data$sample_Ag <- paste(CD8_modality_integrate@meta.data$orig.ident,CD8_modality_integrate_cluster@meta.data$Antigen,sep = "_")
Ind_Ag_per_clus <- table(CD8_modality_integrate@meta.data$sample_Ag,CD8_modality_integrate_cluster@meta.data$seurat_clusters)
write.table(Ind_Ag_per_clus, paste(savedir,"Table/Ind_Ag_per_clus.txt",sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)

DefaultAssay(CD8_modality_integrate) <- "RNA"
Ag_markers <- FindAllMarkers(CD8_modality_integrate)

write.table(Ag_markers, paste(savedir,"Table/Ag_markers.txt",sep = ""), row.names = T, col.names = T, sep = "\t", quote = F)


### Performing the single cell markers young individuals versus all cells from old individuals, per antigen?
CD8_modality_integrate_cluster@meta.data$Age_Ag_clus <- paste(CD8_modality_integrate_cluster@meta.data$Age,CD8_modality_integrate_cluster@meta.data$Antigen, CD8_modality_integrate_cluster@meta.data$seurat_clusters, sep="_")

Ags <- grep(",",CD8_modality_integrate_cluster@meta.data$Antigen,invert=TRUE,value=TRUE) %>% unique()

DefaultAssay(CD8_modality_integrate_cluster) <- "RNA"
clus <- c(2,4)
dir.create(paste(savedir,"scdifferential",sep=""), showWarnings = FALSE)

for(i in 1:length(Ags)){
  group1 = paste("O",Ags[i],"4",sep="_")
  group2 = paste("Y",Ags[i],"4",sep="_")
  O_vs_Y_markers <- FindMarkers(CD8_modality_integrate_cluster, ident.1 = group1, 
                                ident.2 = group2, group.by = "Age_Ag_clus", assay = "RNA", 
                                logfc.threshold = 0)
  write.table(O_vs_Y_markers, paste(savedir,"scdifferential/",group1,"_VS_",group2,".txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
}

### Reference Mapping of each individuals in Ag to Bulk ####
### Using this method for the reference mapping https://satijalab.org/seurat/articles/multimodal_reference_mapping.html 
reference = Bulk
DefaultAssay(reference) = "integrated" 
reference <- RunUMAP(reference, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
              reduction.key = "wnnUMAP_", return.model = TRUE)
p <- DimPlot(reference, group.by = "celltypes", reduction = "wnn.umap") 

pdf("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/reference_celltypes.pdf")
p
dev.off()

# reference=readRDS(paste(savedir,"saveRDS_obj/reference.RDS",sep = ""))

# Computing an sPCA transformation 
# As described in our manuscript, we first compute a ‘supervised’ PCA. This identifies the transformation of the transcriptome data that best encapsulates
# the structure of the WNN graph. This allows a weighted combination of the protein and RNA measurements to ‘supervise’ the PCA, and highlight the
# most relevant sources of variation. After computing this transformation, we can project it onto a query dataset. We can also compute and project a 
# PCA projection, but recommend the use of sPCA when working with multimodal references that have been constructed with WNN analysis.
# The sPCA calculation is performed once, and then can be rapidly projected onto each query dataset.

reference <- ScaleData(reference, assay = 'RNA')
reference <- FindVariableFeatures(reference, nfeatures = 4000, assay = "RNA")
reference <- RunSPCA(reference, assay = 'RNA', graph = 'wsnn')

# Computing a cached neighbor index
# Since we will be mapping multiple query samples to the same reference, we can cache particular steps that only involve the reference. 
# This step is optional but will improve speed when mapping multiple samples.
# We compute the first 50 neighbors in the sPCA space of the reference. We store this information in the spca.annoy.neighbors object within the 
# reference Seurat object and also cache the annoy index data structure (via cache.index = TRUE).

reference <- FindNeighbors(
  object = reference,
  reduction = "spca",
  dims = 1:50,
  graph.name = "spca.annoy.neighbors", 
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

saveRDS(reference,paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Bulk/saveRDS_obj/reference.RDS",sep = ""))
### Optional 
# If you want to save and load a cached index for a Neighbor object generated with method = "annoy" and cache.index = TRUE, use the 
# SaveAnnoyIndex()/LoadAnnoyIndex() functions. Importantly, this index cannot be saved normally to an RDS or RDA file, so it will not persist 
# correctly across R session restarts or saveRDS/readRDS for the Seurat object containing it. Instead, use LoadAnnoyIndex() to add the Annoy index to
# the Neighbor object every time R restarts or you load the reference Seurat object from RDS. The file created by SaveAnnoyIndex() can be distributed
# along with a reference Seurat object, and added to the Neighbor object in the reference.

# bm[["spca.annoy.neighbors"]]

## A Neighbor object containing the 50 nearest neighbors for 30672 cells

# SaveAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "/brahms/shared/vignette-data/reftmp.idx")
# bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "/brahms/shared/vignette-data/reftmp.idx")

# Query dataset preprocessing

# Here we will demonstrate mapping multiple donor bone marrow samples to the multimodal bone marrow reference. These query datasets are derived from 
# the Human Cell Atlas (HCA) Immune Cell Atlas Bone marrow dataset and are available through SeuratData. This dataset is provided as a single merged 
# object with 8 donors. We first split the data back into 8 separate Seurat objects, one for each original donor to map individually.
library(Seurat)
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/"
CD8 <- readRDS(paste(maindir,"saveRDS_obj/CD8_modality_integrate_cluster_0.6.RDS",sep = ""))
query <- CD8
query.batches <- SplitObject(query, split.by = "orig.ident")

# We then normalize the query in the same manner as the reference. Here, the reference was normalized using log-normalization via NormalizeData().
# If the reference had been normalized using SCTransform(), the query must be normalized with SCTransform() as well.

query.batches <- lapply(X = query.batches, FUN = SCTransform, verbose = TRUE)

# Mapping
# We then find anchors between each donor query dataset and the multimodal reference. This command is optimized to minimize mapping time, by passing 
# in a pre-computed set of reference neighbors, and turning off anchor filtration.

anchors <- list()
for (i in 1:length(query.batches)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = reference,
    query = query.batches[[i]],
    normalization.method = "SCT",
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )
}

# We then individually map each of the datasets.
for (i in 1:length(query.batches)) {
  query.batches[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = query.batches[[i]],
    reference = reference, 
    refdata = list(
      celltype = "seurat_clusters2", 
      predicted_ADT = "ADT"),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
}

# Explore the mapping results
# Now that mapping is complete, we can visualize the results for individual objects

p1 <- DimPlot(query.batches[[1]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3) + ggtitle(names(query.batches)[1])
p2 <- DimPlot(query.batches[[2]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3) + ggtitle(names(query.batches)[2])
p3 <- DimPlot(query.batches[[3]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3) + ggtitle(names(query.batches)[3])
p4 <- DimPlot(query.batches[[4]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3) + ggtitle(names(query.batches)[4])
p5 <- DimPlot(query.batches[[5]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3) + ggtitle(names(query.batches)[5])
p6 <- DimPlot(query.batches[[6]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3) + ggtitle(names(query.batches)[6])
p7 <- DimPlot(query.batches[[7]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3) + ggtitle(names(query.batches)[7])
p8 <- DimPlot(query.batches[[8]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3) + ggtitle(names(query.batches)[8])
p9 <- DimPlot(query.batches[[9]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3) + ggtitle(names(query.batches)[9])
p10 <- DimPlot(query.batches[[10]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3) + ggtitle(names(query.batches)[10])
p11 <- DimPlot(query.batches[[11]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3) + ggtitle(names(query.batches)[11])
p12 <- DimPlot(query.batches[[12]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3) + ggtitle(names(query.batches)[12])
p13 <- DimPlot(query.batches[[13]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3) + ggtitle(names(query.batches)[13])
p14 <- DimPlot(query.batches[[14]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3) + ggtitle(names(query.batches)[14])

require(gridExtra)
pdf(paste(savedir,"UMAP/reference_mapping_Ag_to_Bulk.pdf",sep = ""), width = 15, height = 12)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14, ncol = 4, nrow = 4)
dev.off()

saveRDS(query.batches, paste(savedir,"saveRDS_obj/query.batches.RDS",sep = ""))
Ag_ind_clus <- dplyr::select(query@meta.data, seurat_clusters, Antigen, orig.ident)
Ag_req <- grep(",",unique(Ag_ind_clus$Antigen),invert=TRUE,value=TRUE)
Ag_ind_clus_req <- Ag_ind_clus[grep(paste("^",Ag_req,"$",sep="", collapse="|"), Ag_ind_clus$Antigen),]
Ag_ind_clus_req_table <- table(Ag_ind_clus_req$orig.ident, Ag_ind_clus_req$seurat_clusters, Ag_ind_clus_req$Antigen)

### Since this is three dimension table first need to convert to matrix 
Ags <- names(Ag_ind_clus_req_table[1,1,])

for (i in 1:length(Ags)) {
  Ag_table <- Ag_ind_clus_req_table[,,Ags[i]]
  assign(paste(Ags[i], "_table",sep = ""), Ag_table)
}

required_tables <- ls(pattern="[0-9]_table")

rbind(get(required_tables[1]),get(required_tables[2]),get(required_tables[3]),get(required_tables[4]),
      get(required_tables[5]),get(required_tables[6]),get(required_tables[7]),get(required_tables[8]),
      get(required_tables[9]),EBV_EBNA3C_table) -> combined_Ag_table

combined_Ag_table_df <- as.data.frame(combined_Ag_table)

Ag_vector <- c(rep(required_tables[1],14),rep(required_tables[2],14),rep(required_tables[3],14),rep(required_tables[4],14),
               rep(required_tables[5],14),rep(required_tables[6],14),rep(required_tables[7],14),rep(required_tables[8],14),rep(required_tables[9],14),
               rep("EBV_EBNA3C",14))
combined_Ag_table_df$Ags <- gsub("_table","",Ag_vector)

write.table(combined_Ag_table_df, paste(savedir,"Table/combined_Ag_table_df.txt",sep = ""),
            quote = F, col.names = T, row.names = T, sep = "\t")

individuals <- names(query.batches)

for (i in 1:length(individuals)) {
  obj <- query.batches[[i]]
  metadataa=(obj@meta.data)
  metadataa_Ag <- metadataa[grep(",",metadataa$Antigen,invert=TRUE),]
  metadataa_Ag$predicted.celltype <- factor(metadataa_Ag$predicted.celltype, levels=c(0:24))
  Ag_table <- table(metadataa_Ag$Antigen, metadataa_Ag$predicted.celltype)
  write.table(as.matrix(Ag_table), paste(savedir,"Table/",individuals[i],"_Ag_table.txt",sep = ""), 
              quote = F, row.names = T, col.names = T, sep = "\t")
}

### Ag scTCR ####
library(Seurat)
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/"
CD8 <- readRDS(paste(maindir,"saveRDS_obj/CD8_modality_integrate_cluster_0.6.RDS",sep = ""))
savedir <- paste(maindir,"scTCR/",sep = "")
dir.create(savedir, showWarnings = FALSE)

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Table/"
Exp01 <- read.table(paste(savedir,"Exp01_Ag_UMI_gt_50.txt",sep = ""), header = TRUE, sep = "\t")
Exp02 <- read.table(paste(savedir,"Exp02_Ag_UMI_gt_50.txt",sep = ""), header = TRUE, sep = "\t")
Exp03 <- read.table(paste(savedir,"Exp03_Ag_UMI_gt_50.txt",sep = ""), header = TRUE, sep = "\t")

combined <- rbind(Exp01, Exp02, Exp03)
combined$Exp <- gsub(".*._","",combined$samplename)
combined$Exp_barcode <- paste(combined$Exp,"_geno_",combined$barcode,sep="")
combined_noNA <- na.omit(combined)
combined_noNA_ordered <- combined_noNA[match(rownames(CD8@meta.data),combined_noNA$Exp_barcode),]
stopifnot(all(combined_noNA_ordered$Exp_barcode == rownames(CD8@meta.data)))
CD8@meta.data$max_specificity_score <- combined_noNA_ordered$max_specificity_score
CD8@meta.data$max_Ag <- combined_noNA_ordered$max_Ag
CD8@meta.data$Ag_UMI <- combined_noNA_ordered$Ag_UMI
CD8@meta.data$Con_UMI <- combined_noNA_ordered$Con_UMI

Ag_table <- select(CD8@meta.data, max_specificity_score, max_Ag, seurat_clusters, celltypes)
write.table(Ag_table, paste(savedir,"Table/Ag_score_clus_celltype.txt", sep = ""), quote = F, col.names = T, row.names = T, sep = "\t")

p <- VlnPlot(CD8, "max_specificity_score")
pdf("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/QC_Vln/max_specificity_score.pdf")
p
dev.off()

p <- VlnPlot(CD8, "max_specificity_score", group.by = "celltypes")
pdf("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/QC_Vln/max_specificity_score_celltypes.pdf")
p
dev.off()

### Performing the shell commands
# cd /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/preprocessing/Exp01_First_4/post_processing/genotype_demux/Exp01_geno/outs/per_sample_outs/Exp01_geno/vdj_t/
# awk '{print "Exp01_geno_"$0}' filtered_contig_annotations.csv > Exp01_filtered_contig_annotations.csv
# At this location. 1. Combined all the Experiment and removed the header and added it.
cells <- (rownames(CD8@meta.data))
combined_TCR <- read.csv(paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/All_Exp_filtered_contig_annotations_wid_header.csv"), header = TRUE)

CD8_combined_TCR <- combined_TCR[match(rownames(CD8@meta.data), combined_TCR$barcode, nomatch=0),]
dir.create(paste(savedir,"table",sep = ""), showWarnings = FALSE)
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/"
write.table(CD8_combined_TCR, paste(savedir,"Table/CD8_combined_TCR.txt",sep = ""), sep = "\t", row.names = T, col.names = T, quote = F)

CD8_combined_TCR_selected <- select(CD8_combined_TCR,barcode, cdr1,cdr1_nt, cdr2,cdr2_nt, cdr3, cdr3_nt, chain, v_gene, d_gene, j_gene, c_gene, reads, umis)
CD8_combined_TCR_selected$cdr_aa <- paste(CD8_combined_TCR_selected$cdr1, CD8_combined_TCR_selected$cdr2, CD8_combined_TCR_selected$cdr3, sep="_")
CD8_combined_TCR_selected$cdr_nt <- paste(CD8_combined_TCR_selected$cdr1_nt, CD8_combined_TCR_selected$cdr2_nt, CD8_combined_TCR_selected$cdr3_nt, sep="_")

CD8@meta.data$cdr_aa <- CD8_combined_TCR_selected[match(rownames(CD8@meta.data), CD8_combined_TCR_selected$barcode),"cdr_aa"]
CD8@meta.data$cdr_nt <- CD8_combined_TCR_selected[match(rownames(CD8@meta.data), CD8_combined_TCR_selected$barcode),"cdr_nt"]
CD8@meta.data$id_CDR3 <- paste(CD8@meta.data$CDR3,CD8@meta.data$orig.ident,sep="_")

CDR3 <- table(CD8@meta.data$id_CDR3) %>% as.data.frame()
# clonal <- CDR3[order(-CDR3$Freq),]
clonal_noNA <- CDR3[grep("^NA_",CDR3$Var1,invert=TRUE),]
clonal_noNA_gt_2 <- clonal_noNA[clonal_noNA$Freq > 1,]
clonal_noNA_gt_2$Var1 <- factor(clonal_noNA_gt_2$Var1, levels = clonal_noNA_gt_2[order(-clonal_noNA_gt_2$Freq),"Var1"])

library(ggplot2)
dir.create(paste(savedir,"clones",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"clones/CD8_individual_clonal_expanded_gt_1.pdf",sep = ""), width = 60, height = 6)
ggplot(clonal_noNA_gt_2, aes(x=Var1,y=Freq, fill = Var1)) + geom_bar(stat="identity", color="black") +  
  ggtitle("CD8 CDR3 aa gt 1") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 3),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) + NoLegend()
dev.off()

write.table(clonal_noNA_gt_2[order(-clonal_noNA_gt_2$Freq),], paste(savedir,"clones/CD8_individual_clonal_expanded_gt_1.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
# length(clonal_noNA_gt_2$Var1)

### Further we are adding the Antigen 
clonal_noNA_gt_2 <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scTCR/clones/CD8_individual_clonal_expanded_gt_1.txt", header = TRUE)
clonal_noNA_gt_2[,"Ag-num"] <- NA
total_clones <- (clonal_noNA_gt_2$Var1)
for (i in 1:length(total_clones)) {
  clone_spec <- CD8@meta.data[grep(total_clones[i],CD8@meta.data$id_CDR3),"Antigen"] %>% table() %>% as.data.frame()
  colnames(clone_spec) <- c("Antigen","cellnum")
  clone_spec_ordered <- clone_spec[order(-clone_spec$cellnum),]
  clone_spec_ordered$Antigen <- gsub(", ","_",clone_spec_ordered$Antigen)
  clonal_noNA_gt_2[i,"Ag-num"] <- paste(clone_spec_ordered$Antigen, clone_spec_ordered$cellnum,sep="-", collapse=", ")
}

write.table(clonal_noNA_gt_2,paste(savedir,"clones/CD8_individual_clonal_expanded_gt_1_wid_Ag.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

clonal_noNA_gt_2 <- clonal_noNA_gt_2[grep("\\*",clonal_noNA_gt_2$Var1,invert = TRUE),]

rm(plot_list)
plot_list <- list()
for(i in 1:length(clonal_noNA_gt_2$Var1)){
  clonal_expand <- rownames(CD8@meta.data[grep(clonal_noNA_gt_2$Var1[i],CD8@meta.data$id_CDR3),])
  # samplesname <- CD8@meta.data[grep(clonal_noNA_gt_2$Var1[i],CD8@meta.data$id_CDR3),"orig.ident"]  %>% unique()
  # samplesname <- paste(samplesname,collapse=" & ")
  p <- DimPlot(CD8,cells.highlight = clonal_expand, reduction = "wnn.umap", 
               label = TRUE, cols.highlight = "deeppink2",
               cols = "gray92") + 
    ggtitle(paste(clonal_noNA_gt_2$Var1[i])) + 
    NoLegend() + 
    theme(plot.title = element_text(hjust = 0.2))
  # pdf(paste(savedir,"CD8/CD8_",clonal_noNA_gt_2$Var1[i],".pdf",sep = ""))
  # print(p)
  # dev.off()
  plot_list[[i]] <- p
}

dir.create(paste(savedir,"clones/UMAP",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"clones/UMAP/CD8_all_UMAP_combined_1.pdf",sep = ""))
plot_list[1:200]
dev.off()

pdf(paste(savedir,"clones/UMAP/CD8_all_UMAP_combined_2.pdf",sep = ""))
plot_list[201:400]
dev.off()

pdf(paste(savedir,"clones/UMAP/CD8_all_UMAP_combined_3.pdf",sep = ""))
plot_list[401:600]
dev.off()

pdf(paste(savedir,"clones/UMAP/CD8_all_UMAP_combined_4.pdf",sep = ""))
plot_list[601:800]
dev.off()

pdf(paste(savedir,"clones/UMAP/CD8_all_UMAP_combined_5.pdf",sep = ""))
plot_list[801:1000]
dev.off()

pdf(paste(savedir,"clones/UMAP/CD8_all_UMAP_combined_6.pdf",sep = ""))
plot_list[1001:1200]
dev.off()

pdf(paste(savedir,"clones/UMAP/CD8_all_UMAP_combined_7.pdf",sep = ""))
plot_list[1201:1500]
dev.off()

pdf(paste(savedir,"clones/UMAP/CD8_all_UMAP_combined_8.pdf",sep = ""))
plot_list[1501:1767]
dev.off()

### Making a table
CD8_table <- CD8@meta.data[,c("orig.ident","CDR3","seurat_clusters")]

### Removing the naive cluster 1,6, and 12 since they are niave and we are focused on memory cells


paste(CD8_table$orig.ident, CD8_table$seurat_clusters,sep="_") -> CD8_table$id_cluster
id_clone_gt_1 <- gsub("_.*","",clonal_noNA_gt_2$Var1)
CD8_table_gt_1 <- CD8_table[grep(paste(id_clone_gt_1, collapse = "|"), CD8_table$CDR3),]
CDR3_cluster <- table(CD8_table_gt_1$CDR3,CD8_table_gt_1$seurat_clusters)
dir.create(paste(savedir,"CD8/Table/",sep = ""), showWarnings = FALSE)
write.table(CDR3_id_cluster, paste(savedir,"clones/CDR3_clone_cluster_matrix.txt",sep = ""),
            row.names = T, col.names = T, quote = F, sep = "\t")

### CD8 All TRB 
all_TRB <- rownames(CD8@meta.data[!is.na(CD8@meta.data$CDR3),])
p <- DimPlot(CD8,cells.highlight = all_TRB,
             reduction = "wnn.umap", sizes.highlight = 0.2, label = TRUE, cols.highlight = "deeppink2") + ggtitle("CD8_all") +
  NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5))

pdf(paste(savedir,"clones/UMAP/CD8_all.pdf",sep = ""))
p
dev.off()

### CD8 TRB gt 1
all_gt_1 <- rownames(CD8@meta.data[grep(paste(clonal_noNA_gt_2$Var1, collapse = "|"),CD8@meta.data$id_CDR3),])
p <- DimPlot(CD8, cells.highlight = all_gt_1, reduction = "wnn.umap", sizes.highlight = 0.2,
             label = TRUE, cols.highlight = "deeppink2") + ggtitle("CD8_gt_1") +
  NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5))

pdf(paste(savedir,"clones/UMAP/CD8_gt_1.pdf",sep = ""))
p
dev.off()

# ### Unique TCR ######
# #### Identifying the Unique TCR absolute number
# CD8_table <- CD8@meta.data[,c("orig.ident","CDR3","seurat_clusters")]
# CD8_table_mat <- as.matrix(table(CD8_table$CDR3,CD8_table$orig.ident))
# clone_table <- colSums(CD8_table_mat) %>% as.data.frame()
# colnames(clone_table) <- "Total_clones"
# CD8_table_mat[CD8_table_mat > 0] <- 1
# clone_table$unique_clone <- as.data.frame(colSums(CD8_table_mat))[,1]
# clone_table$total_cells <- as.data.frame(table(CD8@meta.data$orig.ident))[,2]
# clone_table$unique_clone_norm <- clone_table$unique_clone/clone_table$total_cells
# clone_table$Total_clones_norm <- clone_table$Total_clones/clone_table$total_cells
# clone_table$Total_clones_norm_with_unique_clones <- clone_table$Total_clones/clone_table$unique_clone
# write.table(clone_table, paste(savedir,"CD8/Table/clone_summary.txt",sep = ""),
#             sep = "\t", row.names = T, col.names = T, quote = F)
# 
# clone_table$sampleid <- rownames(clone_table)
# 
# split_to_dataframe <- function(x) {
#   split_elements <- strsplit(x, "_")[[1]]
#   data.frame(t(split_elements))
# }
# 
# # Convert split elements to dataframe
# sample_split <- do.call(rbind, lapply(clone_table$sampleid, split_to_dataframe))
# rownames(sample_split) <- rownames(clone_table)
# colnames(sample_split) <- c("vaccine","age","ID","Gender","Run")
# clone_table_2 <- merge(clone_table, sample_split, by =0, all=TRUE)
# 
# p <- ggplot(clone_table_2, aes_string(x="sampleid", y="Total_clones", fill="vaccine")) +
#   geom_bar(stat = "identity",color="black", position = "dodge") +
#   ggtitle("Total number of Clones each Sample") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         panel.background = element_rect(fill = 'white', colour = 'black'),
#         panel.grid.minor = element_line(colour = "white"),
#         panel.grid.major = element_line(colour = "white"))
# 
# 
# pdf(paste(savedir,"CD8/Table/total_clones.pdf",sep = ""), width = 12, height = 8)
# p
# dev.off()
# 
# p <- ggplot(clone_table_2, aes_string(x="sampleid", y="unique_clone", fill="vaccine")) +
#   geom_bar(stat = "identity",color="black", position = "dodge") +
#   ggtitle("Unique Clones each Sample") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         panel.background = element_rect(fill = 'white', colour = 'black'),
#         panel.grid.minor = element_line(colour = "white"),
#         panel.grid.major = element_line(colour = "white"))
# 
# 
# pdf(paste(savedir,"CD8/Table/Unique_clones.pdf",sep = ""), width = 12, height = 8)
# p
# dev.off()
# 
# p <- ggplot(clone_table_2, aes_string(x="sampleid", y="unique_clone_norm", fill="vaccine")) +
#   geom_bar(stat = "identity",color="black", position = "dodge") +
#   ggtitle("Unique Clones Normalized each Sample") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         panel.background = element_rect(fill = 'white', colour = 'black'),
#         panel.grid.minor = element_line(colour = "white"),
#         panel.grid.major = element_line(colour = "white"))
# 
# pdf(paste(savedir,"CD8/Table/Unique_clones_Normalized.pdf",sep = ""), width = 12, height = 8)
# p
# dev.off()
# 
# p <- ggplot(clone_table_2, aes_string(x="sampleid", y="Total_clones_norm", fill="vaccine")) +
#   geom_bar(stat = "identity",color="black", position = "dodge") +
#   ggtitle("Total_clones Normalized each Sample") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         panel.background = element_rect(fill = 'white', colour = 'black'),
#         panel.grid.minor = element_line(colour = "white"),
#         panel.grid.major = element_line(colour = "white"))
# 
# 
# pdf(paste(savedir,"CD8/Table/Total_clones_Normalized.pdf",sep = ""), width = 12, height = 8)
# p
# dev.off()
# 
# p <- ggplot(clone_table_2, aes_string(x="sampleid", y="Total_clones_norm_with_unique_clones", fill="vaccine")) +
#   geom_bar(stat = "identity",color="black", position = "dodge") +
#   ggtitle("Total_clones_norm_with_unique_clones Normalized each Sample") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         panel.background = element_rect(fill = 'white', colour = 'black'),
#         panel.grid.minor = element_line(colour = "white"),
#         panel.grid.major = element_line(colour = "white"))
# 
# pdf(paste(savedir,"CD8/Table/Total_clones_norm_with_unique_clones_Normalized.pdf",sep = ""), width = 12, height = 8)
# p
# dev.off()

#### Clonal Sharing #####
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scTCR/"
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/VZV_clonal_sharing.R")
CD8$Age_Ag <- paste(CD8@meta.data$Age ,CD8@meta.data$Antigen, sep="_")
Age_Ag = grep(",",unique(CD8$Age_Ag),value=TRUE,invert=TRUE)

for (i in 1:length(Age_Ag)) {
  TCR_AJ_VZV(object = CD8, savedir = paste(savedir,sep = ""), clus_col = "seurat_clusters", TCR_col = "CDR3",
             group_col = "Age_Ag", group_val = Age_Ag[i], split_col = "orig.ident",
             column_name = c("id","Ag","Age","Age_number","gender","Run"))
}

#### Cell types
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scTCR/celltypes/"
dir.create(savedir, showWarnings = FALSE)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/VZV_clonal_sharing_gt_1_celltypes.R")
CD8$Age_Ag <- paste(CD8@meta.data$Age ,CD8@meta.data$Antigen, sep="_")
Age_Ag = grep(",",unique(CD8$Age_Ag),value=TRUE,invert=TRUE)
clus_length = length(unique(CD8@meta.data$celltypes))
clus_factor = c("naïve","CM","EM1","EM2","EM3","TEMRA","MAIT","HELIOS_high")

for (i in 1:length(Age_Ag)) {
  try({TCR_AJ_VZV_gt_1_CT(object = CD8, savedir = savedir, 
                          clus_col = "celltypes", TCR_col = "CDR3", 
                          group_col = "Age_Ag", group_val = Age_Ag[i], 
                          split_col = "orig.ident", column_name = c("id","Ag","Age","Age_number","gender","Run"),
                          total_clusters = clus_length, clus_fact = clus_factor)
  })
}

### Diversity ####
### Shannon Diversity index or Shannon Entropy ####
rm(df2)
sample_id <- CD8@meta.data$orig.ident %>% unique()
cluster_num <- paste("C",unique(CD8@meta.data$seurat_clusters)[order(unique(CD8@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"diversity",sep = ""),showWarnings = FALSE)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Shannon_Diversity.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Shannon_AJ(object = CD8, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "CDR3",
                             group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = length(cluster_num))
  })
}

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scTCR/"
dir.create(paste(savedir, "clusters",sep = ""), showWarnings = FALSE)
dir.create(paste(savedir, "clusters/diversity/",sep = ""), showWarnings = FALSE)
write.table(df2, paste(savedir,"clusters/diversity/Shannon_sample_diversity_all_cluster.txt",sep = ""), quote = F, 
            row.names = T, col.names = T, sep = "\t")
# df2 <- read.table(paste(savedir,"diversity/Shannon_sample_diversity_all_cluster.txt",sep = ""), header = TRUE)
df2$sample <- row.names(df2)
df2_melted <- melt(df2)
colnames(df2_melted) <- c("Sample","Cluster","Shannon_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name

if("condition" %in% colnames(df2_melted_split) & "vaccine" %in% colnames(df2_melted_split)){
  df2_melted_split[grep("DMSO",df2_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Shannon_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Shannon Diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"clusters/diversity/Age_shannon_diversity_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### Simpson Diversity index #####
rm(df2)
sample_id <- CD8@meta.data$orig.ident %>% unique()
cluster_num <- paste("C",unique(CD8@meta.data$seurat_clusters)[order(unique(CD8@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"diversity",sep = ""),showWarnings = FALSE)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Simpson_Diversity.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Simpson_AJ(object = CD8, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "CDR3",
                             group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = length(cluster_num))
  })
}

write.table(df2, paste(savedir,"clusters/diversity/Simpson_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df2$sample <- row.names(df2)
df2_melted <- melt(df2, id = "sample")
colnames(df2_melted) <- c("Sample","Cluster","Simpson")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name


for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Simpson)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Simpson")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"clusters/diversity/Age_Simpson_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

#### Gini Simpson Diversity index ####
df2 <- read.table(paste(savedir,"clusters/diversity/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1-df2

write.table(df3, paste(savedir,"clusters/diversity/Gini_Simpson_diversity_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df3$sample <- row.names(df3)
df3_melted <- melt(df3, id="sample")
colnames(df3_melted) <- c("Sample","Cluster","Gini_Simpson_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df3_melted_split <- do.call(rbind, lapply(df3_melted[,"Sample"], split_to_dataframe))
colnames(df3_melted_split) <- column_name

if("condition" %in% colnames(df3_melted_split) & "vaccine" %in% colnames(df3_melted_split)){
  df3_melted_split[grep("DMSO",df3_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df3_melted[,column_name[i]] <- df3_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df3_melted), aes(x=Cluster,y=Gini_Simpson_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Gini_Simpson_diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"clusters/diversity/vaccine_Gini_Simpson_diveristy_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

#### Inverse Simpson index #####
df2 <- read.table(paste(savedir,"clusters/diversity/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1/df2

write.table(df3, paste(savedir,"clusters/diversity/inverse_Simpson_diversity_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df3$sample <- row.names(df3)
df3_melted <- melt(df3, id="sample")
colnames(df3_melted) <- c("Sample","Cluster","inverse_Simpson_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df3_melted_split <- do.call(rbind, lapply(df3_melted[,"Sample"], split_to_dataframe))
colnames(df3_melted_split) <- column_name

if("condition" %in% colnames(df3_melted_split) & "vaccine" %in% colnames(df3_melted_split)){
  df3_melted_split[grep("DMSO",df3_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df3_melted[,column_name[i]] <- df3_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df3_melted), aes(x=Cluster,y=inverse_Simpson_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("inverse_Simpson_diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"clusters/diversity/vaccine_inverse_Simpson_diveristy_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### Diversity per person per Antigen ####
### Shannon Diversity index
metadata <- CD8@meta.data
metadata$rows <- rownames(metadata)
metadata <- metadata[!is.na(metadata$CDR3),]
metadata_mem <- metadata[grep("naïve", metadata$celltypes,invert=TRUE),]
metadata_mem_df <- as.data.frame(table(metadata_mem$orig.ident))
metadata_mem_df[order(metadata_mem_df[,2]),2]
### Please choose what cell number you want to keep
cellnumber = 1000 ### The chosen number
metadata_mem_equal <- metadata_mem %>% group_by(orig.ident) %>% slice_sample(n=cellnumber) 

cellnames=metadata_mem_equal$rows
CD8_subset <- subset(CD8, cells = cellnames)


rm(df2)
Ag_sampleid <- grep(",",paste(CD8_subset@meta.data$Antigen,CD8_subset@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
CD8_subset$Ag_sampleid <- paste(CD8_subset@meta.data$Antigen,CD8_subset@meta.data$orig.ident,sep = "_")
cluster_num <- paste("C",unique(CD8_subset@meta.data$seurat_clusters)[order(unique(CD8_subset@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"clusters/diversity",sep = ""),showWarnings = FALSE)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Shannon_Diversity.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Shannon_AJ(object = CD8_subset, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "CDR3",
                             group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = length(cluster_num))
  })
}

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scTCR/"
dir.create(paste(savedir, "clusters",sep = ""), showWarnings = FALSE)
dir.create(paste(savedir, "clusters/diversity/",sep = ""), showWarnings = FALSE)
dir.create(paste(savedir, "clusters/diversity/Ag_individuals_clusters/",sep = ""), showWarnings = FALSE)
write.table(df2, paste(savedir,"clusters/diversity/Ag_individuals_clusters/Shannon_sample_diversity_all_cluster_equal_num.txt",sep = ""), quote = F,
            row.names = T, col.names = T, sep = "\t")

### Simpson Diversity index 
rm(df2)
Ag_sampleid <- grep(",",paste(CD8_subset@meta.data$Antigen,CD8_subset@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
cluster_num <- paste("C",unique(CD8_subset@meta.data$seurat_clusters)[order(unique(CD8_subset@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",cluster_num)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Simpson_Diversity.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Simpson_AJ(object = CD8_subset, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "CDR3",
                             group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = length(cluster_num))
  })
}
write.table(df2, paste(savedir,"clusters/diversity/Ag_individuals_clusters/Simpson_sample_all_cluster_equal_num.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Gini Simpson Diversity index 
# df2 <- read.table(paste(savedir,"clusters/diversity/Ag_individuals_clusters/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1-df2
write.table(df3, paste(savedir,"clusters/diversity/Ag_individuals_clusters/Gini_Simpson_diversity_sample_all_cluster_equal_num.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Inverse Simpson index 
# df2 <- read.table(paste(savedir,"clusters/diversity/Ag_individuals_clusters/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1/df2
write.table(df3, paste(savedir,"clusters/diversity/Ag_individuals_clusters/inverse_Simpson_diversity_sample_all_cluster_equal_num.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

### Celltypes ####
rm(df2)
Ag_sampleid <- grep(",",paste(CD8_subset@meta.data$Antigen,CD8_subset@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
CD8_subset$Ag_sampleid <- paste(CD8_subset@meta.data$Antigen,CD8_subset@meta.data$orig.ident,sep = "_")
celltypes_fact <- c("naïve","CM","EM1","EM2","EM3","TEMRA","MAIT","HELIOS_high")
df2 <- data.frame(matrix(ncol = length(celltypes_fact)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",celltypes_fact)
dir.create(paste(savedir,"celltypes/diversity",sep = ""),showWarnings = FALSE)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Shannon_celltypes.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Shannon_CT_AJ(object = CD8_subset, savedir = savedir, clus_col = "celltypes", TCR_col = "CDR3",
                                group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                                column_name = c("id","Ag","Age","Age_number","gender","Run"), 
                                total_clusters = length(celltypes_fact), 
                                clus_fac = celltypes_fact)
  })
}

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scTCR/"
dir.create(paste(savedir, "celltypes",sep = ""), showWarnings = FALSE)
dir.create(paste(savedir, "celltypes/diversity/",sep = ""), showWarnings = FALSE)
dir.create(paste(savedir, "celltypes/diversity/Ag_individuals_celltypes/",sep = ""), showWarnings = FALSE)
write.table(df2, paste(savedir,"celltypes/diversity/Ag_individuals_celltypes/Shannon_sample_diversity_all_cluster_equal_num.txt",sep = ""), quote = F, 
            row.names = T, col.names = T, sep = "\t")

### Simpson Diversity index 
rm(df2)
Ag_sampleid <- grep(",",paste(CD8_subset@meta.data$Antigen,CD8_subset@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
CD8_subset$Ag_sampleid <- paste(CD8_subset@meta.data$Antigen,CD8_subset@meta.data$orig.ident,sep = "_")
celltypes_fact <- c("naïve","CM","EM1","EM2","EM3","TEMRA","MAIT","HELIOS_high")
df2 <- data.frame(matrix(ncol = length(celltypes_fact)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",celltypes_fact)
dir.create(paste(savedir,"celltypes/diversity",sep = ""),showWarnings = FALSE)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Simpson_celltypes.R")
for (i in 1:1:length(Ag_sampleid)) {
  try({df2[i,] <- Simpson_CT_AJ(object = CD8_subset, savedir = savedir, clus_col = "celltypes", TCR_col = "CDR3",
                             group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = length(celltypes_fact), 
                             clus_fac = celltypes_fact)
  })
}
write.table(df2, paste(savedir,"celltypes/diversity/Ag_individuals_celltypes/Simpson_sample_all_cluster_equal_num.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Gini Simpson Diversity index 
# df2 <- read.table(paste(savedir,"celltypes/diversity/Ag_individuals_celltypes/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1-df2
write.table(df3, paste(savedir,"celltypes/diversity/Ag_individuals_celltypes/Gini_Simpson_diversity_sample_all_cluster_equal_num.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Inverse Simpson index 
# df2 <- read.table(paste(savedir,"celltypes/diversity/Ag_individuals_celltypes/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1/df2
write.table(df3, paste(savedir,"celltypes/diversity/Ag_individuals_celltypes/inverse_Simpson_diversity_sample_all_cluster_equal_num.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")


### Bulk Pseudobulk ####
memory_cells <- rownames(CD8@meta.data[grep("naïve",CD8@meta.data$celltypes,invert=TRUE),])
CD8_mem <- subset(CD8, cells = memory_cells)
mem_cells_specific_Ag <- rownames(CD8_mem@meta.data[grep(",",CD8_mem@meta.data$Age_Ag,invert=TRUE),])
CD8_mem_spec_Ag <- subset(CD8_mem, cells = mem_cells_specific_Ag)
CD8_mem_spec_Ag$Age_Ag_Run <- paste(CD8_mem_spec_Ag$Age_Ag, CD8_mem_spec_Ag$Run,sep="_")

#### PCA
obj <- CD8_mem_spec_Ag
DefaultAssay(CD8_mem_spec_Ag) <- "RNA"
cell_samples <- table(CD8_mem_spec_Ag@meta.data[,paste(grouping_by,sep = "_")]) %>% as.data.frame()
sample_low <- cell_samples[cell_samples$Freq < cell_freq,1]
sample_low <- gsub("_","_",sample_low)
sample_low_2 <- c(sample_low,remove_samples)
message(paste("Removing this sample",sample_low_2,"\n",sep = " "))
Treg_cts <- AggregateExpression(obj, group.by = paste(grouping_by,sep = "_"),
                                assays = 'RNA', slot = "counts", return.seurat = FALSE)
Treg_cts_2 <- as.data.frame(Treg_cts$RNA)
Treg_cts_2 <- Treg_cts_2[,grep("Y_EBV_BMRF1_Exp03|O_EBV_BALF4_Exp02",colnames(Treg_cts_2),invert = TRUE)]
if(length(sample_low_2) != 0){
  index <- grep(paste(sample_low_2,collapse="|"),colnames(Treg_cts_2), invert=TRUE)
  Treg_cts_3 <- Treg_cts_2[,index]
} else {
  Treg_cts_3 <- Treg_cts_2
}


Treg_metadata <- data.frame(samples = colnames(Treg_cts_3))
# Function to split and convert to dataframe
split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, splitting)[[1]]
  data.frame(t(split_elements))
}

# Convert split elements to dataframe
Treg_metadata_2 <- do.call(rbind, lapply(Treg_metadata$samples, split_to_dataframe))
colnames(Treg_metadata_2) <- c(column_name_split)
Treg_metadata_2$orig.ident <- Treg_metadata$samples
# Treg_metadata_2$Run <- Treg_metadata_2$orig.ident

obj@meta.data[,paste(sample_col,batch_col,sep = "_")] <- paste(obj@meta.data[,sample_col],obj@meta.data[,batch_col],sep = "_")
batch_group <- unique(obj@meta.data[,batch_col])
# for (i in 1:length(batch_group)) {
#   donor_change <- unique(obj@meta.data[,paste(sample_col,batch_col,sep = "_")]) %>% grep(batch_group[i],.,value=TRUE) %>% 
#     gsub(paste(batch_group[i],sep = ""),"",.) %>% gsub("_","-",.) %>% paste(.,collapse=".*.|")
#   donor_change <- paste(donor_change,".*.",sep="")
#   Treg_metadata_2$Run <- gsub(donor_change,batch_group[i],Treg_metadata_2$Run)
# }
# Treg_metadata_2$Run <- gsub("-.*.","",Treg_metadata_2$Run)

stopifnot(all(Treg_metadata_2$orig.ident == colnames(Treg_cts_3))) # if TRUE move forward
# Treg_metadata$samplename <- gsub("gE_|DMSO_|_DMSO|_gE|B22015","",Treg_metadata$samples)
## Run will be treated as a covariate in the regression model,

design_AJ <- as.formula(paste("~",paste(colnames(Treg_metadata_2[,c(batch_col,"Age")]),collapse = "+")))
Treg_dds <- DESeqDataSetFromMatrix(countData = Treg_cts_3,
                                   colData = Treg_metadata_2,
                                   design = design_AJ)

keep = filterByExpr(Treg_dds, group=colData(Treg_dds)[,"Age"], min.count=gene_min_counts)
table(keep)
Treg_dds_2 <- Treg_dds[keep,]

rld<-vst(Treg_dds_2)
PCA <- DESeq2::plotPCA(object = rld,
                       intgroup = colnames(Treg_metadata_2),
                       returnData=TRUE,ntop = 5000)
percentVar <- round(100 * attr(PCA, "percentVar"))

savedir4 <- paste(savedir3,"PCA/",sep = "")
dir.create(savedir4, showWarnings = FALSE)


PC1 <-ggplot(PCA, aes_string("PC1", "PC2", label="name", color=batch_col, shape="Age")) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir4,"PCA_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".pdf",sep = ""),
    width = 9, height = 8)
print(PC1)
dev.off()

assay(rld) <- limma::removeBatchEffect(assay(rld),
                                       batch=rld$Run)

PCA_batch <- DESeq2::plotPCA(object = rld,
                             intgroup = colnames(Treg_metadata_2),
                             returnData=TRUE,ntop = 5000)

percentVar <- round(100 * attr(PCA_batch, "percentVar"))

PC1 <-ggplot(PCA_batch, aes_string("PC1", "PC2", label="name", color="Ag", shape="Age")) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir4,"PCA_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".pdf",sep = ""),
    width = 9, height = 8)
print(PC1)
dev.off()

write.table(PCA_batch, paste(savedir4,"PCA_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".txt",sep = ""),
            quote = F, row.names = T, col.names = T, sep = "\t")

PC1 <-ggplot(PCA_batch_2, aes_string("PC1", "PC2", label="Age", color="Ag", shape=batch_col)) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir4,"PCA_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,"_remove_outliers.pdf",sep = ""),
    width = 9, height = 8)
print(PC1)
dev.off()

PC1 <-ggplot(PCA_batch_2, aes_string("PC1", "PC2", label="name", color="Ag", shape="Age")) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir4,"PCA_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,"_remove_outliers_no_Run_marked.pdf",sep = ""),
    width = 9, height = 8)
print(PC1)
dev.off()


#### RNA
memory_cells <- rownames(CD8@meta.data[grep("naïve",CD8@meta.data$celltypes,invert=TRUE),])
CD8_mem <- subset(CD8, cells = memory_cells)
mem_cells_specific_Ag <- rownames(CD8_mem@meta.data[grep(",",CD8_mem@meta.data$Age_Ag,invert=TRUE),])
CD8_mem_spec_Ag <- subset(CD8_mem, cells = mem_cells_specific_Ag)
CD8_mem_spec_Ag@meta.data$Antigen2 <- gsub(".*._","",CD8_mem_spec_Ag@meta.data$Antigen)
CD8_mem_spec_Ag@meta.data$AgeAntigen2 <- paste(CD8_mem_spec_Ag@meta.data$Age,CD8_mem_spec_Ag@meta.data$Antigen2,sep=".")

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)

Ag = unique(CD8_mem_spec_Ag@meta.data$Antigen2)
AgeAg <- unique(CD8_mem_spec_Ag@meta.data$AgeAntigen2)

for (i in 1:length(Ag)) {
  try({
    cluster = "all"
    group1 = paste("O.",Ag[i],sep = "")
    group2 = paste("Y.",Ag[i],sep = "")
    remove_samples = grep(Ag[i],AgeAg,value = TRUE, invert = TRUE)
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
    colnames(design0) <- c("O","Y",paste("Run",1:(ncol(design0)-2),sep = ""))
    cm <- makeContrasts(O_VS_Y = O-Y,levels = design0)
    desl_clus <- LimmaEdgeR_differential(dds = CD8_subset_dds,
                                         design0 = design0,
                                         cm = cm, 
                                         savedir = savedir2,
                                         logfc = 0.5,
                                         p_value_adj = 0.05)
  })
}

obj <- CD8_mem_spec_Ag
DefaultAssay(obj) <- "RNA"
obj@meta.data$Ag_sampleid <- paste(obj@meta.data$Antigen, obj@meta.data$orig.ident, sep = "_")
cell_samples <- table(obj@meta.data[,"Ag_sampleid"]) %>% as.data.frame()
sample_low <- cell_samples[cell_samples$Freq < 20,1]
sample_low2 <- gsub("_","_",sample_low)

Treg_cts <- AggregateExpression(obj, group.by = paste("Ag_sampleid"), assays = 'RNA', slot = "counts", return.seurat = FALSE)
Treg_cts_2 <- as.data.frame(Treg_cts$RNA)
if(length(sample_low2) != 0){
  index <- grep(paste(sample_low2,collapse="|"),colnames(Treg_cts_2), invert=TRUE)
  Treg_cts_3 <- Treg_cts_2[,index]
} else {
  Treg_cts_3 <- Treg_cts_2
}

Treg_metadata <- data.frame(samples = colnames(Treg_cts_3))
# Function to split and convert to dataframe
split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, splitting)[[1]]
  data.frame(t(split_elements))
}

# Convert split elements to dataframe
Treg_metadata_2 <- do.call(rbind, lapply(Treg_metadata$samples, split_to_dataframe))
colnames(Treg_metadata_2) <- c("virus_type","Ag",column_name_split)
Treg_metadata_2$orig.ident <- Treg_metadata$samples
# Treg_metadata_2$Run <- Treg_metadata_2$orig.ident

sample_col = "orig.ident"
obj@meta.data[,paste(sample_col,batch_col,sep = "_")] <- paste(obj@meta.data[,sample_col],obj@meta.data[,batch_col],sep = "_")
batch_group <- unique(obj@meta.data[,batch_col])

stopifnot(all(Treg_metadata_2$orig.ident == colnames(Treg_cts_3))) # if TRUE move forward
# Treg_metadata$samplename <- gsub("gE_|DMSO_|_DMSO|_gE|B22015","",Treg_metadata$samples)
## Run will be treated as a covariate in the regression model,

design_AJ <- as.formula(paste("~",paste(colnames(Treg_metadata_2[,c(batch_col,"Ag")]),collapse = "+")))
Treg_dds <- DESeqDataSetFromMatrix(countData = Treg_cts_3,
                                   colData = Treg_metadata_2,
                                   design = design_AJ)

keep = filterByExpr(Treg_dds, group=colData(Treg_dds)[,"Ag"], min.count=gene_min_counts)
table(keep)
Treg_dds_2 <- Treg_dds[keep,]
rld <- vst(Treg_dds_2)
PCA <- DESeq2::plotPCA(object = rld,
                       intgroup = colnames(Treg_metadata_2),
                       returnData=TRUE,ntop = 5000)
percentVar <- round(100 * attr(PCA, "percentVar"))

savedir3 <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen_range_10/pseudobulk/all_Ag_CD8_mem/"
dir.create(savedir3, showWarnings = FALSE)
savedir4 <- paste(savedir3,"PCA/",sep = "")
dir.create(savedir4, showWarnings = FALSE)

PC1 <-ggplot(PCA, aes_string("PC1", "PC2", label="name", color="Ag", shape="Run")) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf(paste(savedir4,"PCA_cell_gt_all_individual_",cell_freq,".pdf",sep = ""),
    width = 9, height = 8)
print(PC1)
dev.off()

write.table(PCA, paste(savedir4,"PCA_with_batch_effect_cell_gt_all_individuals.txt",sep = ""),
            quote = F, row.names = T, col.names = T, sep = "\t")

assay(rld) <- limma::removeBatchEffect(assay(rld),
                                       batch=rld$Run)


# PCA_batch <- DESeq2::plotPCA(object = rld,
#                              intgroup = colnames(Treg_metadata_2),
#                              returnData=TRUE,ntop = 5000)
# 
# percentVar <- round(100 * attr(PCA_batch, "percentVar"))
# 
# PC1 <-ggplot(PCA_batch, aes_string("PC1", "PC2", label="name", color="Ag", shape="Run")) +
#   geom_point(size=5, alpha=0.7) +
#   geom_text_repel(size=2) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.background = element_rect(fill = 'white', colour = 'black'),
#         panel.grid.minor = element_line(colour = "grey"),
#         panel.grid.major = element_line(colour = "grey"))
# 
# pdf(paste(savedir4,"PCA_removed_batch_effect_cell_gt_all_individual_",cell_freq,".pdf",sep = ""),
#     width = 9, height = 8)
# print(PC1)
# dev.off()
# write.table(PCA_batch, paste(savedir4,"PCA_removed_batch_effect_cell_gt_all_individuals.txt",sep = ""),
#             quote = F, row.names = T, col.names = T, sep = "\t")

#### Further we are interested in all the PCs PC1 to all
pc <- prcomp(t(assay(rld)))
# attributes(pc)
loading <- pc$rotation 
summary <- summary(pc)

rownames(Treg_metadata_2) <- Treg_metadata_2$orig.ident
pc_with_metadata <- merge(pc$x,Treg_metadata_2,by = 'row.names', all = TRUE)

write.table(pc_with_metadata, paste(savedir4,"PCA_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,"_all_PCs.txt",sep = ""),
            quote = F, row.names = T, col.names = T, sep = "\t")

rm(plot_list)
plot_list <- list()

plot_list[[1]] <-ggplot(pc_with_metadata, aes_string("PC1", "PC2", label="orig.ident", color="Ag", shape="age")) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(round(as.vector(summary$importance[2,][1]),4)*100),"% variance")) +
  ylab(paste0("PC2: ",(round(as.vector(summary$importance[2,][2]),4)*100),"% variance")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

plot_list[[2]] <-ggplot(pc_with_metadata, aes_string("PC3", "PC4", label="orig.ident", color="Ag", shape="age")) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(round(as.vector(summary$importance[2,][3]),4)*100),"% variance")) +
  ylab(paste0("PC2: ",(round(as.vector(summary$importance[2,][4]),4)*100),"% variance")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

plot_list[[3]] <-ggplot(pc_with_metadata, aes_string("PC5", "PC6", label="orig.ident", color="Ag", shape="age")) +
  geom_point(size=5, alpha=0.7) +
  geom_text_repel(size=2) +
  xlab(paste0("PC1: ",(round(as.vector(summary$importance[2,][5]),4)*100),"% variance")) +
  ylab(paste0("PC2: ",(round(as.vector(summary$importance[2,][6]),4)*100),"% variance")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))


pdf(paste(savedir4,"PC1_to_PC6_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,"_fullname.pdf",sep = ""),
    width = 9, height = 8)
print(plot_list)
dev.off()

saveRDS(pc,paste(savedir4,"pca_prcomp.RDS",sep = ""))

log_cpm_all <- cpm(assay(Treg_dds_2), prior.count=2, log=TRUE)

write.table(log_cpm_all, paste("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen_range_10/pseudobulk/all_Ag_CD8_mem/log_CPM_all_individuals.txt",sep = ""),
            quote = F, row.names = T, col.names = T, sep = "\t")

#### Using the dataset to make the files with the range of +/- 10
CD8_Ag_table_ind <- table(paste(CD8@meta.data$orig.ident,CD8@meta.data$Ag_range_10), CD8@meta.data$seurat_clusters)

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/"
write.table(CD8_Ag_table_ind, paste(savedir,"Table/CD8_Ag_table_ind.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")


#### Decided that Ag +/-  10 #####
## Pseudobulk of all memory cells per single-specific antigen (not multi-specific/cross reactive) per person --> differential expression in each antigen comparing O and Y groups, only taking samples >20 cells; showing gene expression table of results as well as CPM across all conditions, so that they can be compared across antigen
#### RNA
library(Seurat)
maindir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/"
# CD8 <- readRDS(paste(maindir, "saveRDS_obj/CD8_modality_integrate_cluster_0.6.RDS",sep = ""))
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

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Range_10_rem_TRB2_rem_batch_effect/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)

Ag = unique(CD8_mem_spec_Ag@meta.data$Antigen2)
AgeAg <- unique(CD8_mem_spec_Ag@meta.data$AgeAntigen2)

for (i in 1:length(Ag)) {
  try({
    cluster = "all"
    group1 = paste("O.",Ag[i],sep = "")
    group2 = paste("Y.",Ag[i],sep = "")
    remove_samples = grep(Ag[i],AgeAg,value = TRUE, invert = TRUE)
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
    colnames(design0) <- c("O","Y",paste("Run",1:(ncol(design0)-2),sep = ""))
    cm <- makeContrasts(O_VS_Y = O-Y,levels = design0)
    desl_clus <- LimmaEdgeR_differential(dds = CD8_subset_dds,
                                         design0 = design0,
                                         cm = cm,
                                         savedir = savedir2,
                                         logfc = 0.5,
                                         p_value_adj = 0.05)
  })
}

### No batch effect removed
for (i in 1:length(Ag)) {
  try({
    cluster = "all"
    group1 = paste("O.",Ag[i],sep = "")
    group2 = paste("Y.",Ag[i],sep = "")
    remove_samples = grep(Ag[i],AgeAg,value = TRUE, invert = TRUE)
    remove_samples_name <- paste(remove_samples,collapse="_and_")
    DefaultAssay(CD8_mem_spec_Ag) <- "RNA"
    # savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run01_to_10_analysis/vaccine/S_vs_Z/"
    source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/COVID_pseudobulk_PCA_within_AJ_2.R")
    CD8_subset_dds <- COVID_pseudobulk_within_cluster_AJ(obj = CD8_mem_spec_Ag, savedir = savedir, group1 = group1, group2 = group2,
                                                         grouping_by = "AgeAntigen2", cluster = cluster, cell_freq = 20, remove_samples = remove_samples,
                                                         cluster_group = "seurat_clusters",sample_col = "orig.ident", batch_col = "Run",
                                                         gene_min_counts =  5,
                                                         column_name_split = c("sampleid","virus","age","age_number","gender","Run"))
    
    cluster2 <- paste(cluster, sep="_", collapse="_")
    savedir2 <- paste(savedir,"pseudobulk/clus_",cluster2,"_",group1,"_vs_",group2,"/",sep = "")
    # savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run04/Analysis/pseudobulk/clus_13_removed__DMSO_vs_Sprotein/"
    source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/RNAseq/pipeline_function/RNAseq_limma_EdgeR.R")
    design0 <- model.matrix(~ 0 + AgeAntigen2, data = colData(CD8_subset_dds))
    colnames(design0) <- c("O","Y")
    cm <- makeContrasts(O_VS_Y = O-Y,levels = design0)
    desl_clus <- LimmaEdgeR_differential(dds = CD8_subset_dds,
                                         design0 = design0,
                                         cm = cm, 
                                         savedir = savedir2,
                                         logfc = 0.5,
                                         p_value_adj = 0.05)
  })
}

#### TCR counting 
CD8@meta.data$TRB1 <- gsub(".*.--","",CD8@meta.data$TRA1_TRA2_TRB1_TRB2_cdrnt) %>% gsub("-.*.","",.)
CD8@meta.data$TRA1 <- gsub("--.*.","",CD8@meta.data$TRA1_TRA2_TRB1_TRB2_cdrnt) %>% gsub("-.*.","",.)
CD8@meta.data$TRA2 <- gsub("--.*.","",CD8@meta.data$TRA1_TRA2_TRB1_TRB2_cdrnt) %>% gsub(".*.-","",.)

CD8@meta.data$TRB1 <- na_if(CD8@meta.data$TRB1, "NA")
CD8@meta.data$TRA1 <- na_if(CD8@meta.data$TRA1, "NA")
CD8@meta.data$TRA2 <- na_if(CD8@meta.data$TRA2, "NA")

### We need to flip the TRA1 and TRA2
CD8@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc <- na_if(CD8@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc,"NA")

CD8@meta.data$TRA1_vdjc <- gsub("--.*.","",CD8@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc) %>% gsub("TRAC*-.*.","",.)
CD8@meta.data$TRA2_vdjc <- gsub("--.*.","",CD8@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc) %>% gsub(".*.-TRA","",.)

CD8@meta.data$TRA1_vdjc <- na_if(CD8@meta.data$TRA1_vdjc,"NA")
CD8@meta.data$TRA2_vdjc <- na_if(CD8@meta.data$TRA2_vdjc,"NA")

CD8@meta.data$TRA1_vdjc <- gsub("_.*.","",CD8@meta.data$TRA1_vdjc) %>% gsub("TRAV","",.) %>% gsub("-","",.) 
CD8@meta.data$TRA2_vdjc <- gsub("_.*.","",CD8@meta.data$TRA2_vdjc) %>% gsub("V","",.) %>% gsub("-","",.)

CD8@meta.data$TRA1_vdjc <- gsub("/.*.","",CD8@meta.data$TRA1_vdjc)
CD8@meta.data$TRA2_vdjc <- gsub("/.*.","",CD8@meta.data$TRA2_vdjc)

CD8@meta.data$TRA1_j <- gsub("--.*.","",CD8@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc) %>% gsub("TRAC*-.*.","",.) %>% gsub(".*?(TRAJ[0-9]+).*", "\\1",.)
CD8@meta.data$TRA2_j <- gsub("--.*.","",CD8@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc) %>% gsub(".*.-TRA","",.) %>% gsub(".*?(TRAJ[0-9]+).*", "\\1",.)

CD8@meta.data$TRA1_j <- gsub("TRAJ","",CD8@meta.data$TRA1_j)
CD8@meta.data$TRA2_j <- gsub("TRAJ","",CD8@meta.data$TRA2_j)

CD8@meta.data$TRA1_ordered<- CD8@meta.data$TRA1
CD8@meta.data$TRA2_ordered<- CD8@meta.data$TRA2

### Finding the cells with both TRA1 and TRA2
TRA1_TRA2_barcode <- rownames(CD8@meta.data[!is.na(CD8@meta.data$TRA2),])

for (i in 1:length(TRA1_TRA2_barcode)) {
  cell_index <- grep(TRA1_TRA2_barcode[i],rownames(CD8@meta.data))
  TRV1 <- as.numeric(CD8@meta.data[cell_index,"TRA1_vdjc"])
  TRV2 <- as.numeric(CD8@meta.data[cell_index,"TRA2_vdjc"])
  TRJ1 <- as.numeric(CD8@meta.data[cell_index,"TRA1_j"])
  TRJ2 <- as.numeric(CD8@meta.data[cell_index,"TRA2_j"])
  
  if (TRV1 == TRV2) {
    if (TRJ1 == TRJ2){
      TCRA1_length = str_length(CD8@meta.data[cell_index,"TRA1"])
      TCRA2_length = str_length(CD8@meta.data[cell_index,"TRA2"])
      if (TCRA1_length < TCRA2_length) { ### Since we have to keep in descending order to get more nt
        CD8@meta.data[cell_index,"TRA1_ordered"] <- CD8@meta.data[cell_index,"TRA2"]
        CD8@meta.data[cell_index,"TRA2_ordered"] <- CD8@meta.data[cell_index,"TRA1"]
        # print(CD8@meta.data[cell_index,c("TRA1","TRA2","TRA1_vdjc","TRA2_vdjc","TRA1_j","TRA2_j","TRA1_ordered","TRA2_ordered")])
      }
    }
    else if(TRJ1 > TRJ2) { ### We are keeping it in the ascending order
      print("TRJ1 > TRJ2")
      CD8@meta.data[cell_index,"TRA1_ordered"] <- CD8@meta.data[cell_index,"TRA2"]
      CD8@meta.data[cell_index,"TRA2_ordered"] <- CD8@meta.data[cell_index,"TRA1"]
      # print(CD8@meta.data[cell_index,c("TRA1","TRA2","TRA1_vdjc","TRA2_vdjc","TRA1_j","TRA2_j","TRA1_ordered","TRA2_ordered")])
    }
  }
  if(TRV1 > TRV2){
    print("TRV1 > TRV2")
    CD8@meta.data[cell_index,"TRA1_ordered"] <- CD8@meta.data[cell_index,"TRA2"]
    CD8@meta.data[cell_index,"TRA2_ordered"] <- CD8@meta.data[cell_index,"TRA1"]
    # print(CD8@meta.data[cell_index,c("TRA1","TRA2","TRA1_vdjc","TRA2_vdjc","TRA1_j","TRA2_j","TRA1_ordered","TRA2_ordered")])
  }
}

CD8@meta.data$TRA1_TRA2_TRB1_TRB2_ordered <- paste(CD8@meta.data$TRA1_ordered,"-" ,CD8@meta.data$TRA2_ordered,"--", CD8@meta.data$TRB1,"-",CD8@meta.data$TRB2 ,sep="")

clones <- CD8@meta.data$TRA1_TRA2_TRB1_TRB2_ordered
clone_no_TRB2 <- grep("-NA$",clones,value=TRUE)

clone_no_TRB2_df <- table(clone_no_TRB2) %>% as.data.frame()
table(clone_no_TRB2_df$Freq > 1)
table(clone_no_TRB2_df$Freq > 5)
table(clone_no_TRB2_df$Freq > 50)

clone_no_TRB2_with_TRB1 <- table(grep("NA-NA$",clone_no_TRB2, value=TRUE, invert=TRUE)) %>% as.data.frame()
write.table(clone_no_TRB2_with_TRB1,"/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1.txt", sep="\t",quote=F, col.names=T, row.names=F)

clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2 <- grep("^NA",clone_no_TRB2, value=TRUE, invert=TRUE) %>% grep("--NA-NA$",.,invert=TRUE, value=TRUE) %>% table() %>% as.data.frame()

### Performing single specific and multi specifc clones with cross reactivity
clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2 <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2.txt", header = TRUE)
CD8@meta.data$TRA1_or_TRA2_TRB1_no_TRB2 <- clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2[match(CD8@meta.data$TRA1_TRA2_TRB1_TRB2_ordered,clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2[,1]),1]
TCR_freq_df <- table(CD8@meta.data$TRA1_or_TRA2_TRB1_no_TRB2) %>% as.data.frame()

### Start with the > 50 expanded clones
TCR_freq_gt_50 <- TCR_freq_df[TCR_freq_df$Freq > 50,]
TCR_freq_gt_50[,1] <- as.character(TCR_freq_gt_50[,1]) ## To remove the factor
single_spec_gt_50 = 0
multi_spec_gt_50 = 0

### Creating an empty dataframe
single_spec_gt_50_df <- data.frame(matrix(nrow=dim(TCR_freq_gt_50)[1],ncol=6))
multi_spec_gt_50_df <- data.frame(matrix(nrow=dim(TCR_freq_gt_50)[1],ncol=6))
colnames(single_spec_gt_50_df) <- c("Clone","Ag","Ag_number","Ag_dom_percent","cross_reactive_percent","no_cross_reactive_percent")
colnames(multi_spec_gt_50_df) <- c("Clone","Ag","Ag_number","Ag_dom_percent","cross_reactive_percent","no_cross_reactive_percent")

for (i in 1:dim(TCR_freq_gt_50)[1]) {
  clone_spec_Ag <- table(CD8@meta.data[grep(paste("^",TCR_freq_gt_50[i,1],"$",sep=""),CD8@meta.data$TRA1_or_TRA2_TRB1_no_TRB2),"Ag_range_10"]) %>% as.data.frame()
  
  ### Checking whether it is multi-specific or single specific
  ## for multi specific conditions are 
  # Multi specific clones = cells of a clone have 2+ antigen specificities within max score -10%, (with or without true cross-reactivity, within max score -10%), 
  # e.g. cell 1 Ag 1 only, cell 2 Ag 2 only. for large clones (5+ cells) a margin of 5% (rounded up) is applied in which multi-specificity is not called 
  # e.g. of 10 cells 9 are Ag 1, 1 cell is Ag2 --> designated single specific, not multi 
  
  clone_spec_Ag_ordered <- clone_spec_Ag[order(-clone_spec_Ag$Freq),]
  Ag_dom <- clone_spec_Ag_ordered[1,1]
  Ag_dom_size <- length(grep(Ag_dom, clone_spec_Ag_ordered[,1]))
  total_Ag <- dim(clone_spec_Ag_ordered)[1]
  
  ### Ag_dominate percent
  Ag_dom_sum <- sum(clone_spec_Ag_ordered[grep(Ag_dom,clone_spec_Ag_ordered$Var1),2])
  total_Ag_sum <- sum(clone_spec_Ag_ordered$Freq)
  Ag_dom_percent <- (Ag_dom_sum/total_Ag_sum)*100
  
  ## Calculating Ag and its frequency
  for (j in 1:length(total_Ag)) {
    Ag1 <- paste(clone_spec_Ag_ordered$Var1,collapse="-")
    Ag1_freq <- paste(clone_spec_Ag_ordered$Freq,collapse="-")
  }
  
  ## Cross Reactivity
  X_reactive_total <- sum(clone_spec_Ag_ordered[grep(", ",clone_spec_Ag_ordered$Var1),"Freq"])
  X_reactive_percent <- (X_reactive_total/total_Ag_sum)*100
  no_X_reactive_percent <- 100-X_reactive_percent
  
  # if (Ag_dom_size == total_Ag) {
  if (Ag_dom_percent >= 95) {
    single_spec_gt_50 = single_spec_gt_50+1
    single_spec_gt_50_df[i,1] <- TCR_freq_gt_50[i,1]
    single_spec_gt_50_df[i,2] <- Ag1
    single_spec_gt_50_df[i,3] <- Ag1_freq
    single_spec_gt_50_df[i,4] <- Ag_dom_percent
    single_spec_gt_50_df[i,5] <- X_reactive_percent
    single_spec_gt_50_df[i,6] <- no_X_reactive_percent
  } else if (Ag_dom_percent < 95){
    # else if (Ag_dom_size < total_Ag){
    multi_spec_gt_50 = multi_spec_gt_50+1
    multi_spec_gt_50_df[i,1] <- TCR_freq_gt_50[i,1]
    multi_spec_gt_50_df[i,2] <- Ag1
    multi_spec_gt_50_df[i,3] <- Ag1_freq
    multi_spec_gt_50_df[i,4] <- Ag_dom_percent
    multi_spec_gt_50_df[i,5] <- X_reactive_percent
    multi_spec_gt_50_df[i,6] <- no_X_reactive_percent
  }
}

single_spec_gt_50_df_noNA <- na.omit(single_spec_gt_50_df)
multi_spec_gt_50_df_noNA <- na.omit(multi_spec_gt_50_df)

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/Table/"
write.table(single_spec_gt_50_df_noNA, paste(savedir,"clone_gt_50/single_spec_gt_50_df.txt",sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(multi_spec_gt_50_df_noNA, paste(savedir,"clone_gt_50/multi_spec_gt_50_df.txt",sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")


### Doing the same for > 5
### Start with the > 50 expanded clones
TCR_freq_gt_5 <- TCR_freq_df[TCR_freq_df$Freq > 5,]
TCR_freq_gt_5[,1] <- as.character(TCR_freq_gt_5[,1]) ## To remove the factor
single_spec_gt_5 = 0
multi_spec_gt_5 = 0

### Creating an empty dataframe
single_spec_gt_5_df <- data.frame(matrix(nrow=dim(TCR_freq_gt_5)[1],ncol=6))
multi_spec_gt_5_df <- data.frame(matrix(nrow=dim(TCR_freq_gt_5)[1],ncol=6))
colnames(single_spec_gt_5_df) <- c("Clone","Ag","Ag_number","Ag_dom_percent","cross_reactive_percent","no_cross_reactive_percent")
colnames(multi_spec_gt_5_df) <- c("Clone","Ag","Ag_number","Ag_dom_percent","cross_reactive_percent","no_cross_reactive_percent")

for (i in 1:dim(TCR_freq_gt_5)[1]) {
  clone_spec_Ag <- table(CD8@meta.data[grep(paste("^",TCR_freq_gt_5[i,1],"$",sep=""),CD8@meta.data$TRA1_or_TRA2_TRB1_no_TRB2),"Ag_range_10"]) %>% as.data.frame()
  
  ### Checking whether it is multi-specific or single specific
  ## for multi specific conditions are 
  # Multi specific clones = cells of a clone have 2+ antigen specificities within max score -10%, (with or without true cross-reactivity, within max score -10%), 
  # e.g. cell 1 Ag 1 only, cell 2 Ag 2 only. for large clones (5+ cells) a margin of 5% (rounded up) is applied in which multi-specificity is not called 
  # e.g. of 10 cells 9 are Ag 1, 1 cell is Ag2 --> designated single specific, not multi 
  
  clone_spec_Ag_ordered <- clone_spec_Ag[order(-clone_spec_Ag$Freq),]
  Ag_dom <- clone_spec_Ag_ordered[1,1]
  Ag_dom_size <- length(grep(Ag_dom, clone_spec_Ag_ordered[,1]))
  total_Ag <- dim(clone_spec_Ag_ordered)[1]
  
  ### Ag_dominate percent
  Ag_dom_sum <- sum(clone_spec_Ag_ordered[grep(Ag_dom,clone_spec_Ag_ordered$Var1),2])
  total_Ag_sum <- sum(clone_spec_Ag_ordered$Freq)
  Ag_dom_percent <- (Ag_dom_sum/total_Ag_sum)*100
  
  ## Calculating Ag and its frequency
  for (j in 1:length(total_Ag)) {
    Ag1 <- paste(clone_spec_Ag_ordered$Var1,collapse="-")
    Ag1_freq <- paste(clone_spec_Ag_ordered$Freq,collapse="-")
  }
  
  ## Cross Reactivity
  X_reactive_total <- sum(clone_spec_Ag_ordered[grep(", ",clone_spec_Ag_ordered$Var1),"Freq"])
  X_reactive_percent <- (X_reactive_total/total_Ag_sum)*100
  no_X_reactive_percent <- 100-X_reactive_percent
  
  # if (Ag_dom_size == total_Ag) {
  if (Ag_dom_percent >= 95) {
    single_spec_gt_5 = single_spec_gt_5+1
    single_spec_gt_5_df[i,1] <- TCR_freq_gt_5[i,1]
    single_spec_gt_5_df[i,2] <- Ag1
    single_spec_gt_5_df[i,3] <- Ag1_freq
    single_spec_gt_5_df[i,4] <- Ag_dom_percent
    single_spec_gt_5_df[i,5] <- X_reactive_percent
    single_spec_gt_5_df[i,6] <- no_X_reactive_percent
  } else if (Ag_dom_percent < 95){
    # else if (Ag_dom_size < total_Ag){
    multi_spec_gt_5 = multi_spec_gt_5+1
    multi_spec_gt_5_df[i,1] <- TCR_freq_gt_5[i,1]
    multi_spec_gt_5_df[i,2] <- Ag1
    multi_spec_gt_5_df[i,3] <- Ag1_freq
    multi_spec_gt_5_df[i,4] <- Ag_dom_percent
    multi_spec_gt_5_df[i,5] <- X_reactive_percent
    multi_spec_gt_5_df[i,6] <- no_X_reactive_percent
  }
}

single_spec_gt_5_df_noNA <- na.omit(single_spec_gt_5_df)
multi_spec_gt_5_df_noNA <- na.omit(multi_spec_gt_5_df)

savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/Table/"
write.table(single_spec_gt_5_df_noNA, paste(savedir,"clone_gt_5/single_spec_gt_5_df.txt",sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(multi_spec_gt_5_df_noNA, paste(savedir,"clone_gt_5/multi_spec_gt_5_df.txt",sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")

#### TCR ####
### TCR Diversity
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/"
CD8@meta.data$TRA1_or_TRA2_TRB1_no_TRB2_samples <- paste(CD8@meta.data$TRA1_or_TRA2_TRB1_no_TRB2,CD8@meta.data$orig.ident,sep="_")
CDR3 <- table(CD8@meta.data$TRA1_or_TRA2_TRB1_no_TRB2_samples) %>% as.data.frame()
# clonal <- CDR3[order(-CDR3$Freq),]
clonal_noNA <- CDR3[grep("^NA_",CDR3$Var1,invert=TRUE),]
clonal_noNA_gt_2 <- clonal_noNA[clonal_noNA$Freq > 1,]
clonal_noNA_gt_2$Var1 <- factor(clonal_noNA_gt_2$Var1, levels = clonal_noNA_gt_2[order(-clonal_noNA_gt_2$Freq),"Var1"])

library(ggplot2)
dir.create(paste(savedir,"clones",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"clones/CD8_individual_clonal_expanded_gt_1.pdf",sep = ""), width = 60, height = 15)
ggplot(clonal_noNA_gt_2, aes(x=Var1,y=Freq, fill = Var1)) + geom_bar(stat="identity", color="black") +  
  ggtitle("CD8 CDR3 aa gt 1") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 3),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) + NoLegend()
dev.off()

write.table(clonal_noNA_gt_2[order(-clonal_noNA_gt_2$Freq),], paste(savedir,"clones/CD8_individual_clonal_expanded_gt_1.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
# length(clonal_noNA_gt_2$Var1)

### Further we are adding the Antigen 
clonal_noNA_gt_2 <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/clones/CD8_individual_clonal_expanded_gt_1.txt", header = TRUE)
clonal_noNA_gt_2[,"Ag-num"] <- NA
total_clones <- (clonal_noNA_gt_2$Var1)
for (i in 1:length(total_clones)) {
  clone_spec <- CD8@meta.data[grep(total_clones[i],CD8@meta.data$TRA1_or_TRA2_TRB1_no_TRB2_samples),"Ag_range_10"] %>% table() %>% as.data.frame()
  colnames(clone_spec) <- c("Antigen","cellnum")
  clone_spec_ordered <- clone_spec[order(-clone_spec$cellnum),]
  clone_spec_ordered$Antigen <- gsub(", ","_",clone_spec_ordered$Antigen)
  clonal_noNA_gt_2[i,"Ag-num"] <- paste(clone_spec_ordered$Antigen, clone_spec_ordered$cellnum,sep="-", collapse=", ")
}

write.table(clonal_noNA_gt_2,paste(savedir,"clones/CD8_individual_clonal_expanded_gt_1_wid_Ag.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")
clonal_noNA_gt_2 <- clonal_noNA_gt_2[grep("\\*",clonal_noNA_gt_2$Var1,invert = TRUE),]
CD8@meta.data$TRA1_or_TRA2_TRB1_no_TRB2_gt_2 <- clonal_noNA_gt_2[match(CD8@meta.data$TRA1_or_TRA2_TRB1_no_TRB2, clonal_noNA_gt_2$Var1),"Var1"]

CD8@meta.data$TRA1_or_TRA2_TRB1_no_TRB2_gt_2 <- na_if(CD8@meta.data$TRA1_or_TRA2_TRB1_no_TRB2_gt_2, "NA")
clone_gt_2 <- rownames(CD8@meta.data[!is.na(CD8@meta.data$TRA1_or_TRA2_TRB1_no_TRB2_gt_2),])

p <- DimPlot(CD8,cells.highlight = clone_gt_2, reduction = "wnn.umap",
             label = TRUE, cols.highlight = "deeppink2", cols = "gray92", sizes.highlight = 0.1) +
  ggtitle("clone_gt_2 wid ind") +
  NoLegend() +
  theme(plot.title = element_text(hjust = 0.5))

dir.create(paste(savedir,"UMAP",sep = ""), showWarnings = FALSE)
pdf(paste(savedir,"UMAP/clones_gt_2.pdf",sep = ""))
p
dev.off()


### Removing the naive cluster 1,6, and 12 since they are niave and we are focused on memory cells


#### Clonal Sharing #####
## Ag
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/"
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/VZV_clonal_sharing.R")
CD8$Age_Ag <- paste(CD8@meta.data$Age ,CD8@meta.data$Ag_range_10, sep="_")
Age_Ag = grep(",",unique(CD8$Age_Ag),value=TRUE,invert=TRUE)

for (i in 1:length(Age_Ag)) {
  TCR_AJ_VZV(object = CD8, savedir = paste(savedir,sep = ""), clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2_gt_2",
             group_col = "Age_Ag", group_val = Age_Ag[i], split_col = "orig.ident",
             column_name = c("id","Ag","Age","Age_number","gender","Run"))
}

#### Cell types
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/"
dir.create(savedir, showWarnings = FALSE)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/VZV_clonal_sharing_gt_1_celltypes.R")
CD8$Age_Ag <- paste(CD8@meta.data$Age ,CD8@meta.data$Ag_range_10, sep="_")
Age_Ag = grep(",",unique(CD8$Age_Ag),value=TRUE,invert=TRUE)
clus_length = length(unique(CD8@meta.data$celltypes))
clus_factor = c("naïve","CM","EM1","EM2","EM3","TEMRA","MAIT","HELIOS_high")

for (i in 1:length(Age_Ag)) {
  try({TCR_AJ_VZV_gt_1_CT(object = CD8, savedir = savedir, 
                          clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2_gt_2",
                          group_col = "Age_Ag", group_val = Age_Ag[i], 
                          split_col = "orig.ident", column_name = c("id","Ag","Age","Age_number","gender","Run"),
                          total_clusters = clus_length, clus_fact = clus_factor)
  })
}

#### Individuals
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/"
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/VZV_clonal_sharing.R")
CD8$sample_Ag <- paste(CD8@meta.data$orig.ident ,CD8@meta.data$Ag_range_10, sep="_")
sample_Ag = grep(",",unique(CD8$sample_Ag),value=TRUE,invert=TRUE)

for (i in 1:length(Age_Ag)) {
  TCR_AJ_VZV(object = CD8, savedir = paste(savedir,sep = ""), clus_col = "seurat_clusters", 
             TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2_gt_2",
             group_col = "sample_Ag", group_val = sample_Ag[i], split_col = "orig.ident",
             column_name = c("id","Ag","Age","Age_number","gender","Run"))
}

#### Cell types
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/"
dir.create(savedir, showWarnings = FALSE)
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/VZV_clonal_sharing_gt_1_celltypes.R")
CD8$sample_Ag <- paste(CD8@meta.data$orig.ident ,CD8@meta.data$Ag_range_10, sep="_")
sample_Ag = grep(",",unique(CD8$sample_Ag),value=TRUE,invert=TRUE)
clus_length = length(unique(CD8@meta.data$celltypes))
clus_factor = c("naïve","CM","EM1","EM2","EM3","TEMRA","MAIT","HELIOS_high")

for (i in 1:length(Age_Ag)) {
  try({TCR_AJ_VZV_gt_1_CT(object = CD8, savedir = savedir, 
                          clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2_gt_2",
                          group_col = "sample_Ag", group_val = sample_Ag[i], 
                          split_col = "orig.ident", column_name = c("id","Ag","Age","Age_number","gender","Run"),
                          total_clusters = clus_length, clus_fact = clus_factor)
  })
}

### All Diversity ####
### Shannon Diversity index or Shannon Entropy ####
rm(df2)
sample_id <- CD8@meta.data$orig.ident %>% unique()
cluster_num <- paste("C",unique(CD8@meta.data$seurat_clusters)[order(unique(CD8@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"diversity",sep = ""),showWarnings = FALSE)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Shannon_Diversity.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Shannon_AJ(object = CD8, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                             group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = length(cluster_num))
  })
}

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/diversity/"
dir.create(paste(savedir, "clusters",sep = ""), showWarnings = FALSE)
write.table(df2, paste(savedir,"clusters/Shannon_sample_diversity_all_cluster.txt",sep = ""), quote = F, 
            row.names = T, col.names = T, sep = "\t")
# df2 <- read.table(paste(savedir,"diversity/Shannon_sample_diversity_all_cluster.txt",sep = ""), header = TRUE)
df2$sample <- row.names(df2)
df2_melted <- melt(df2)
colnames(df2_melted) <- c("Sample","Cluster","Shannon_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

column_name = c("id","Ag","Age","Age_number","gender","Run")

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name

if("condition" %in% colnames(df2_melted_split) & "vaccine" %in% colnames(df2_melted_split)){
  df2_melted_split[grep("DMSO",df2_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Shannon_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Shannon Diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"clusters/Age_shannon_diversity_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### Simpson Diversity index #####
rm(df2)
sample_id <- CD8@meta.data$orig.ident %>% unique()
cluster_num <- paste("C",unique(CD8@meta.data$seurat_clusters)[order(unique(CD8@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"diversity",sep = ""),showWarnings = FALSE)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Simpson_Diversity.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Simpson_AJ(object = CD8, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                             group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = length(cluster_num))
  })
}

write.table(df2, paste(savedir,"clusters/Simpson_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df2$sample <- row.names(df2)
df2_melted <- melt(df2, id = "sample")
colnames(df2_melted) <- c("Sample","Cluster","Simpson")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name


for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Simpson)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Simpson")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"clusters/Age_Simpson_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

#### Gini Simpson Diversity index ####
df2 <- read.table(paste(savedir,"clusters/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1-df2

write.table(df3, paste(savedir,"clusters/Gini_Simpson_diversity_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df3$sample <- row.names(df3)
df3_melted <- melt(df3, id="sample")
colnames(df3_melted) <- c("Sample","Cluster","Gini_Simpson_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df3_melted_split <- do.call(rbind, lapply(df3_melted[,"Sample"], split_to_dataframe))
colnames(df3_melted_split) <- column_name

if("condition" %in% colnames(df3_melted_split) & "vaccine" %in% colnames(df3_melted_split)){
  df3_melted_split[grep("DMSO",df3_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df3_melted[,column_name[i]] <- df3_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df3_melted), aes(x=Cluster,y=Gini_Simpson_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Gini_Simpson_diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"clusters/vaccine_Gini_Simpson_diveristy_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

#### Inverse Simpson index #####
df2 <- read.table(paste(savedir,"clusters/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1/df2

write.table(df3, paste(savedir,"clusters/dinverse_Simpson_diversity_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df3$sample <- row.names(df3)
df3_melted <- melt(df3, id="sample")
colnames(df3_melted) <- c("Sample","Cluster","inverse_Simpson_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df3_melted_split <- do.call(rbind, lapply(df3_melted[,"Sample"], split_to_dataframe))
colnames(df3_melted_split) <- column_name

if("condition" %in% colnames(df3_melted_split) & "vaccine" %in% colnames(df3_melted_split)){
  df3_melted_split[grep("DMSO",df3_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df3_melted[,column_name[i]] <- df3_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df3_melted), aes(x=Cluster,y=inverse_Simpson_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("inverse_Simpson_diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"clusters/vaccine_inverse_Simpson_diveristy_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### All Celltypes ####
rm(df2)
sample_id <- CD8@meta.data$orig.ident %>% unique()
celltypes_num <- unique(CD8@meta.data$celltypes)[order(unique(CD8@meta.data$celltypes))]
df2 <- data.frame(matrix(ncol = length(celltypes_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",celltypes_num)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Shannon_Diversity.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Shannon_AJ(object = CD8, savedir = savedir, clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                             group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = celltypes_num)
  })
}

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/diversity/all_cells/"
dir.create(paste(savedir, "celltypes",sep = ""), showWarnings = FALSE)
write.table(df2, paste(savedir,"celltypes/Shannon_sample_diversity_all_celltypes_gt_20.txt",sep = ""), quote = F, 
            row.names = T, col.names = T, sep = "\t")

# df2 <- read.table(paste(savedir,"diversity/Shannon_sample_diversity_all_celltypes.txt",sep = ""), header = TRUE)
df2$sample <- row.names(df2)
df2_melted <- melt(df2)
colnames(df2_melted) <- c("Sample","celltypes","Shannon_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

column_name = c("id","Ag","Age","Age_number","gender","Run")

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name

if("condition" %in% colnames(df2_melted_split) & "vaccine" %in% colnames(df2_melted_split)){
  df2_melted_split[grep("DMSO",df2_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=celltypes,y=Shannon_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Shannon Diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"celltypes/Age_shannon_diversity_celltypess_noNA_gt_20.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### Simpson Diversity index #####
rm(df2)
sample_id <- CD8@meta.data$orig.ident %>% unique()
celltypes_num <- unique(CD8@meta.data$celltypes)[order(unique(CD8@meta.data$celltypes))]
df2 <- data.frame(matrix(ncol = length(celltypes_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",celltypes_num)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Simpson_Diversity.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Simpson_AJ(object = CD8, savedir = savedir, clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                             group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = celltypes_num)
  })
}

write.table(df2, paste(savedir,"celltypes/Simpson_sample_all_celltypes_gt_20.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df2$sample <- row.names(df2)
df2_melted <- melt(df2, id = "sample")
colnames(df2_melted) <- c("Sample","celltypes","Simpson")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name


for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=celltypes,y=Simpson)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Simpson")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"celltypes/Age_Simpson_celltypess_noNA_gt_20.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

#### Gini Simpson Diversity index ####
df2 <- read.table(paste(savedir,"celltypes/Simpson_sample_all_celltypes_gt_20.txt",sep = ""), header = TRUE)
df3=1-df2

write.table(df3, paste(savedir,"celltypes/Gini_Simpson_diversity_sample_all_celltypes_gt_20.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df3$sample <- row.names(df3)
df3_melted <- melt(df3, id="sample")
colnames(df3_melted) <- c("Sample","celltypes","Gini_Simpson_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df3_melted_split <- do.call(rbind, lapply(df3_melted[,"Sample"], split_to_dataframe))
colnames(df3_melted_split) <- column_name

if("condition" %in% colnames(df3_melted_split) & "vaccine" %in% colnames(df3_melted_split)){
  df3_melted_split[grep("DMSO",df3_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df3_melted[,column_name[i]] <- df3_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df3_melted), aes(x=celltypes,y=Gini_Simpson_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Gini_Simpson_diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"celltypes/vaccine_Gini_Simpson_diveristy_celltypess_noNA_gt_20.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

#### Inverse Simpson index #####
df2 <- read.table(paste(savedir,"celltypes/Simpson_sample_all_celltypes_gt_20.txt",sep = ""), header = TRUE)
df3=1/df2

write.table(df3, paste(savedir,"celltypes/dinverse_Simpson_diversity_sample_all_celltypes_gt_20.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df3$sample <- row.names(df3)
df3_melted <- melt(df3, id="sample")
colnames(df3_melted) <- c("Sample","celltypes","inverse_Simpson_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df3_melted_split <- do.call(rbind, lapply(df3_melted[,"Sample"], split_to_dataframe))
colnames(df3_melted_split) <- column_name

if("condition" %in% colnames(df3_melted_split) & "vaccine" %in% colnames(df3_melted_split)){
  df3_melted_split[grep("DMSO",df3_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df3_melted[,column_name[i]] <- df3_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df3_melted), aes(x=celltypes,y=inverse_Simpson_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("inverse_Simpson_diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"celltypes/vaccine_inverse_Simpson_diveristy_celltypess_noNA_gt_20.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### mem Diversity ####
### Shannon Diversity index or Shannon Entropy ####
rm(df2)
sample_id <- CD8_mem@meta.data$orig.ident %>% unique()
cluster_num <- paste("C",unique(CD8_mem@meta.data$seurat_clusters)[order(unique(CD8_mem@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"diversity",sep = ""),showWarnings = FALSE)
rm(sha_ent)
cluster_name <- unique(CD8_mem@meta.data$seurat_clusters)[order(unique(CD8_mem@meta.data$seurat_clusters))]

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Shannon_Diversity.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Shannon_AJ(object = CD8_mem, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                             group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/diversity/mem_cells/"
dir.create(paste(savedir, "clusters",sep = ""), showWarnings = FALSE)
write.table(df2, paste(savedir,"clusters/Shannon_sample_diversity_all_cluster.txt",sep = ""), quote = F, 
            row.names = T, col.names = T, sep = "\t")
# df2 <- read.table(paste(savedir,"diversity/Shannon_sample_diversity_all_cluster.txt",sep = ""), header = TRUE)
df2$sample <- row.names(df2)
df2_melted <- melt(df2)
colnames(df2_melted) <- c("Sample","Cluster","Shannon_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

column_name = c("id","Ag","Age","Age_number","gender","Run")

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name

if("condition" %in% colnames(df2_melted_split) & "vaccine" %in% colnames(df2_melted_split)){
  df2_melted_split[grep("DMSO",df2_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Shannon_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Shannon Diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"clusters/Age_shannon_diversity_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### Simpson Diversity index #####
rm(df2)
sample_id <- CD8_mem@meta.data$orig.ident %>% unique()
cluster_num <- paste("C",unique(CD8_mem@meta.data$seurat_clusters)[order(unique(CD8_mem@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"diversity",sep = ""),showWarnings = FALSE)
rm(sha_ent)
cluster_name <- unique(CD8_mem@meta.data$seurat_clusters)[order(unique(CD8_mem@meta.data$seurat_clusters))]

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Simpson_Diversity.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Simpson_AJ(object = CD8_mem, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                             group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

write.table(df2, paste(savedir,"clusters/Simpson_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df2$sample <- row.names(df2)
df2_melted <- melt(df2, id = "sample")
colnames(df2_melted) <- c("Sample","Cluster","Simpson")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name


for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Simpson)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Simpson")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"clusters/Age_Simpson_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

#### Gini Simpson Diversity index ####
df2 <- read.table(paste(savedir,"clusters/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1-df2

write.table(df3, paste(savedir,"clusters/Gini_Simpson_diversity_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df3$sample <- row.names(df3)
df3_melted <- melt(df3, id="sample")
colnames(df3_melted) <- c("Sample","Cluster","Gini_Simpson_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df3_melted_split <- do.call(rbind, lapply(df3_melted[,"Sample"], split_to_dataframe))
colnames(df3_melted_split) <- column_name

if("condition" %in% colnames(df3_melted_split) & "vaccine" %in% colnames(df3_melted_split)){
  df3_melted_split[grep("DMSO",df3_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df3_melted[,column_name[i]] <- df3_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df3_melted), aes(x=Cluster,y=Gini_Simpson_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Gini_Simpson_diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"clusters/vaccine_Gini_Simpson_diveristy_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

#### Inverse Simpson index #####
df2 <- read.table(paste(savedir,"clusters/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1/df2

write.table(df3, paste(savedir,"clusters/dinverse_Simpson_diversity_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df3$sample <- row.names(df3)
df3_melted <- melt(df3, id="sample")
colnames(df3_melted) <- c("Sample","Cluster","inverse_Simpson_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df3_melted_split <- do.call(rbind, lapply(df3_melted[,"Sample"], split_to_dataframe))
colnames(df3_melted_split) <- column_name

if("condition" %in% colnames(df3_melted_split) & "vaccine" %in% colnames(df3_melted_split)){
  df3_melted_split[grep("DMSO",df3_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df3_melted[,column_name[i]] <- df3_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df3_melted), aes(x=Cluster,y=inverse_Simpson_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("inverse_Simpson_diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"clusters/vaccine_inverse_Simpson_diveristy_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### memory Celltypes ####
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/diversity/mem_cells/"
CD8_mem <- subset(CD8, idents = c(0,2,3,4,5,7,8,9,10,11,13,14,15,16,17)) ### removing the naive clusters
sample_id <- CD8_mem@meta.data$orig.ident %>% unique()
celltypes_num <- unique(CD8_mem@meta.data$celltypes)[order(unique(CD8_mem@meta.data$celltypes))]
df2 <- data.frame(matrix(ncol = length(celltypes_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",celltypes_num)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Shannon_Diversity.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Shannon_AJ(object = CD8_mem, savedir = savedir, clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                                group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                                column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = celltypes_num)
  })
}

dir.create(paste(savedir, "celltypes",sep = ""), showWarnings = FALSE)
write.table(df2, paste(savedir,"celltypes/Shannon_sample_diversity_all_celltypes_gt_20.txt",sep = ""), quote = F, 
            row.names = T, col.names = T, sep = "\t")
# df2 <- read.table(paste(savedir,"diversity/Shannon_sample_diversity_all_celltypes.txt",sep = ""), header = TRUE)
df2$sample <- row.names(df2)
df2_melted <- melt(df2)
colnames(df2_melted) <- c("Sample","celltypes","Shannon_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

column_name = c("id","Ag","Age","Age_number","gender","Run")

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name

if("condition" %in% colnames(df2_melted_split) & "vaccine" %in% colnames(df2_melted_split)){
  df2_melted_split[grep("DMSO",df2_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=celltypes,y=Shannon_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Shannon Diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"celltypes/Age_shannon_diversity_celltypess_noNA_gt_20.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### Simpson Diversity index #####
rm(df2)
sample_id <- CD8_mem@meta.data$orig.ident %>% unique()
celltypes_num <- unique(CD8_mem@meta.data$celltypes)[order(unique(CD8_mem@meta.data$celltypes))]
df2 <- data.frame(matrix(ncol = length(celltypes_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",celltypes_num)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Simpson_Diversity.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Simpson_AJ(object = CD8_mem, savedir = savedir, clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                                group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                                column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = celltypes_num)
  })
}

write.table(df2, paste(savedir,"celltypes/Simpson_sample_all_celltypes_gt_20.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df2$sample <- row.names(df2)
df2_melted <- melt(df2, id = "sample")
colnames(df2_melted) <- c("Sample","celltypes","Simpson")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name


for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=celltypes,y=Simpson)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Simpson")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"celltypes/Age_Simpson_celltypess_noNA_gt_20.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

#### Gini Simpson Diversity index ####
df2 <- read.table(paste(savedir,"celltypes/Simpson_sample_all_celltypes_gt_20.txt",sep = ""), header = TRUE)
df3=1-df2

write.table(df3, paste(savedir,"celltypes/Gini_Simpson_diversity_sample_all_celltypes_gt_20.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df3$sample <- row.names(df3)
df3_melted <- melt(df3, id="sample")
colnames(df3_melted) <- c("Sample","celltypes","Gini_Simpson_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df3_melted_split <- do.call(rbind, lapply(df3_melted[,"Sample"], split_to_dataframe))
colnames(df3_melted_split) <- column_name

if("condition" %in% colnames(df3_melted_split) & "vaccine" %in% colnames(df3_melted_split)){
  df3_melted_split[grep("DMSO",df3_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df3_melted[,column_name[i]] <- df3_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df3_melted), aes(x=celltypes,y=Gini_Simpson_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Gini_Simpson_diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"celltypes/vaccine_Gini_Simpson_diveristy_celltypess_noNA_gt_20.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

#### Inverse Simpson index #####
df2 <- read.table(paste(savedir,"celltypes/Simpson_sample_all_celltypes_gt_20.txt",sep = ""), header = TRUE)
df3=1/df2

write.table(df3, paste(savedir,"celltypes/dinverse_Simpson_diversity_sample_all_celltypes_gt_20.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df3$sample <- row.names(df3)
df3_melted <- melt(df3, id="sample")
colnames(df3_melted) <- c("Sample","celltypes","inverse_Simpson_diversity")

split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df3_melted_split <- do.call(rbind, lapply(df3_melted[,"Sample"], split_to_dataframe))
colnames(df3_melted_split) <- column_name

if("condition" %in% colnames(df3_melted_split) & "vaccine" %in% colnames(df3_melted_split)){
  df3_melted_split[grep("DMSO",df3_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df3_melted[,column_name[i]] <- df3_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df3_melted), aes(x=celltypes,y=inverse_Simpson_diversity)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("inverse_Simpson_diversity")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"celltypes/vaccine_inverse_Simpson_diveristy_celltypess_noNA_gt_20.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### Diversity per person per Antigen ####
### Shannon Diversity index
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/diversity/mem_cell_Ag/"
dir.create(savedir, showWarnings = FALSE)
rm(df2)
Ag_sampleid <- grep(",",paste(CD8_mem@meta.data$Ag_range_10,CD8_mem@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
CD8_mem$Ag_sampleid <- paste(CD8_mem@meta.data$Ag_range_10,CD8_mem@meta.data$orig.ident,sep = "_")
cluster_num <- paste("C",unique(CD8_mem@meta.data$seurat_clusters)[order(unique(CD8_mem@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"clusters",sep = ""),showWarnings = FALSE)
rm(sha_ent)
cluster_name <- unique(CD8_mem@meta.data$seurat_clusters)[order(unique(CD8_mem@meta.data$seurat_clusters))]

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Shannon_Diversity.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Shannon_AJ(object = CD8_mem, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                             group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/diversity/mem_cell_Ag/"
dir.create(paste(savedir, "clusters",sep = ""), showWarnings = FALSE)
write.table(df2, paste(savedir,"clusters/Shannon_sample_diversity.txt",sep = ""), quote = F,
            row.names = T, col.names = T, sep = "\t")

### Simpson Diversity index 
rm(df2)
Ag_sampleid <- grep(",",paste(CD8_mem@meta.data$Ag_range_10,CD8_mem@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
cluster_num <- paste("C",unique(CD8_mem@meta.data$seurat_clusters)[order(unique(CD8_mem@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",cluster_num)
rm(sha_ent)
cluster_name <- unique(CD8_mem@meta.data$seurat_clusters)[order(unique(CD8_mem@meta.data$seurat_clusters))]

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Simpson_Diversity.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Simpson_AJ(object = CD8_mem, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                             group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}
write.table(df2, paste(savedir,"clusters/Simpson_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Gini Simpson Diversity index 
# df2 <- read.table(paste(savedir,"clusters/diversity/Ag_individuals_clusters/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1-df2
write.table(df3, paste(savedir,"clusters/Gini_Simpson_diversity_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Inverse Simpson index 
# df2 <- read.table(paste(savedir,"clusters/diversity/Ag_individuals_clusters/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1/df2
write.table(df3, paste(savedir,"clusters/inverse_Simpson_diversity_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

### Celltypes ####
rm(df2)
Ag_sampleid <- grep(",",paste(CD8_mem@meta.data$Ag_range_10,CD8_mem@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
CD8_mem$Ag_sampleid <- paste(CD8_mem@meta.data$Ag_range_10,CD8_mem@meta.data$orig.ident,sep = "_")
celltypes_fact <- c("CM","EM1","EM2","EM3","TEMRA","MAIT","HELIOS_high")
df2 <- data.frame(matrix(ncol = length(celltypes_fact)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",celltypes_fact)
dir.create(paste(savedir,"celltypes",sep = ""),showWarnings = FALSE)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Shannon_Diversity.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Shannon_AJ(object = CD8_mem, savedir = savedir, clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                                group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                                column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = celltypes_fact)
  })
}

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/diversity/mem_cell_Ag/"
dir.create(paste(savedir, "celltypes",sep = ""), showWarnings = FALSE)
write.table(df2, paste(savedir,"celltypes/Shannon_sample_diversity_all_cluster_gt_20.txt",sep = ""), quote = F, 
            row.names = T, col.names = T, sep = "\t")

### Simpson Diversity index 
rm(df2)
Ag_sampleid <- grep(",",paste(CD8_mem@meta.data$Ag_range_10,CD8_mem@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
cluster_num <- c("CM","EM1","EM2","EM3","TEMRA","MAIT","HELIOS_high")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",cluster_num)
rm(sha_ent)
cluster_name <- cluster_num

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Simpson_Diversity.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Simpson_AJ(object = CD8_mem, savedir = savedir, clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                             group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

write.table(df2, paste(savedir,"celltypes/Simpson_sample_all_cluster_gt_20.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Gini Simpson Diversity index 
# df2 <- read.table(paste(savedir,"celltypes/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1-df2
write.table(df3, paste(savedir,"celltypes/Gini_Simpson_diversity_sample_all_celltypes_gt_20.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Inverse Simpson index 
# df2 <- read.table(paste(savedir,"celltypes/diversity/Ag_individuals_celltypes/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1/df2
write.table(df3, paste(savedir,"celltypes/inverse_Simpson_diversity_sample_gt_20.txt",sep = ""), 
            quote = F, row.names = T, col.names = T, sep = "\t")

#### Gini Clonality Index #####
### All Cells 
rm(df2)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/clonality/all_cells/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
sample_id <- CD8@meta.data$orig.ident %>% unique()
cluster_num <- paste("C",unique(CD8@meta.data$seurat_clusters)[order(unique(CD8@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"expansion",sep = ""),showWarnings = FALSE)
rm(sha_ent)
cluster_name <- gsub("C","",cluster_num)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Gini_coefficient_index.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Gini_coef_AJ(object = CD8, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                               group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                               column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

write.table(df2, paste(savedir,"expansion/Gini_coef_clonality_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df2$sample <- row.names(df2)
df2_melted <- melt(df2, id = "sample")
colnames(df2_melted) <- c("Sample","Cluster","Gini_coef")

column_name = c("id","Ag","Age","Age_number","gender","Run")
split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name

if("condition" %in% colnames(df2_melted_split) & "vaccine" %in% colnames(df2_melted_split)){
  df2_melted_split[grep("DMSO",df2_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Gini_coef)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Gini coefficient Clonality")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"expansion/vaccine_Gini_coef_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### Celltypes
rm(df2)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/clonality/all_cells/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
sample_id <- CD8@meta.data$orig.ident %>% unique()
CT_name <- unique(CD8@meta.data$celltypes)
cluster_num <- CT_name[order(CT_name)]
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"expansion",sep = ""),showWarnings = FALSE)
rm(sha_ent)
cluster_name <- cluster_num

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Gini_coefficient_index.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Gini_coef_AJ(object = CD8, savedir = savedir, clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                               group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                               column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

write.table(df2, paste(savedir,"expansion/Gini_coef_clonality_sample_all_celltypes.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df2$sample <- row.names(df2)
df2_melted <- melt(df2, id = "sample")
colnames(df2_melted) <- c("Sample","Cluster","Gini_coef")

column_name = c("id","Ag","Age","Age_number","gender","Run")
split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name

if("condition" %in% colnames(df2_melted_split) & "vaccine" %in% colnames(df2_melted_split)){
  df2_melted_split[grep("DMSO",df2_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Gini_coef)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Gini coefficient Clonality")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"expansion/vaccine_Gini_coef_celltypes_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### Naive Diversity per person per Antigen ####
### Shannon Diversity index
library(Seurat)
library(dplyr)
CD8 <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/saveRDS_obj/CD8_modality_integrate_cluster_0.6.RDS")
savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/diversity/naive_cell_Ag/"
dir.create(savedir, showWarnings = FALSE)
rm(df2)
naive_cells <- rownames(CD8@meta.data[grep("naïve",CD8@meta.data$celltypes),])
CD8_naive <- subset(CD8, cells = naive_cells)
Ag_sampleid <- grep(",",paste(CD8_naive@meta.data$Ag_range_10,CD8_naive@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
CD8_naive$Ag_sampleid <- paste(CD8_naive@meta.data$Ag_range_10,CD8_naive@meta.data$orig.ident,sep = "_")
cluster_num <- paste("C",unique(CD8_naive@meta.data$seurat_clusters)[order(unique(CD8_naive@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"clusters",sep = ""),showWarnings = FALSE)
rm(sha_ent)
cluster_name <- unique(CD8_naive@meta.data$seurat_clusters)[order(unique(CD8_naive@meta.data$seurat_clusters))]

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Shannon_Diversity.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Shannon_AJ(object = CD8_naive, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                             group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/diversity/naive_cell_Ag/"
dir.create(paste(savedir, "clusters",sep = ""), showWarnings = FALSE)
write.table(df2, paste(savedir,"clusters/Shannon_sample_diversity.txt",sep = ""), quote = F,
            row.names = T, col.names = T, sep = "\t")

### Simpson Diversity index 
rm(df2)
Ag_sampleid <- grep(",",paste(CD8_naive@meta.data$Ag_range_10,CD8_naive@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
cluster_num <- paste("C",unique(CD8_naive@meta.data$seurat_clusters)[order(unique(CD8_naive@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",cluster_num)
rm(sha_ent)
cluster_name <- unique(CD8_naive@meta.data$seurat_clusters)[order(unique(CD8_naive@meta.data$seurat_clusters))]

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Simpson_Diversity.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Simpson_AJ(object = CD8_naive, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                             group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}
write.table(df2, paste(savedir,"clusters/Simpson_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Gini Simpson Diversity index 
# df2 <- read.table(paste(savedir,"clusters/diversity/Ag_individuals_clusters/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1-df2
write.table(df3, paste(savedir,"clusters/Gini_Simpson_diversity_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Inverse Simpson index 
# df2 <- read.table(paste(savedir,"clusters/diversity/Ag_individuals_clusters/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1/df2
write.table(df3, paste(savedir,"clusters/inverse_Simpson_diversity_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

### Celltypes ####
rm(df2)
Ag_sampleid <- grep(",",paste(CD8_naive@meta.data$Ag_range_10,CD8_naive@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
CD8_naive$Ag_sampleid <- paste(CD8_naive@meta.data$Ag_range_10,CD8_naive@meta.data$orig.ident,sep = "_")
celltypes_fact <- c("naïve")
df2 <- data.frame(matrix(ncol = length(celltypes_fact)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",celltypes_fact)
dir.create(paste(savedir,"celltypes",sep = ""),showWarnings = FALSE)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Shannon_Diversity.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Shannon_AJ(object = CD8_naive, savedir = savedir, clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                             group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = celltypes_fact)
  })
}

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/diversity/naive_cell_Ag/"
dir.create(paste(savedir, "celltypes",sep = ""), showWarnings = FALSE)
write.table(df2, paste(savedir,"celltypes/Shannon_sample_diversity_all_cluster_gt_20.txt",sep = ""), quote = F, 
            row.names = T, col.names = T, sep = "\t")

### Simpson Diversity index 
rm(df2)
Ag_sampleid <- grep(",",paste(CD8_naive@meta.data$Ag_range_10,CD8_naive@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
celltypes_fact <- c("naïve")
df2 <- data.frame(matrix(ncol = length(celltypes_fact)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",celltypes_fact)
rm(sha_ent)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Simpson_Diversity.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Simpson_AJ(object = CD8_naive, savedir = savedir, clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                             group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = celltypes_fact)
  })
}

write.table(df2, paste(savedir,"celltypes/Simpson_sample_all_cluster_gt_20.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Gini Simpson Diversity index 
# df2 <- read.table(paste(savedir,"celltypes/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1-df2
write.table(df3, paste(savedir,"celltypes/Gini_Simpson_diversity_sample_all_celltypes_gt_20.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Inverse Simpson index 
# df2 <- read.table(paste(savedir,"celltypes/diversity/Ag_individuals_celltypes/Simpson_sample_all_cluster.txt",sep = ""), header = TRUE)
df3=1/df2
write.table(df3, paste(savedir,"celltypes/inverse_Simpson_diversity_sample_gt_20.txt",sep = ""), 
            quote = F, row.names = T, col.names = T, sep = "\t")

#### Gini Clonality Index #####
### Naive cells
rm(df2)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/clonality/all_cells/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
sample_id <- CD8@meta.data$orig.ident %>% unique()
cluster_num <- paste("C",unique(CD8@meta.data$seurat_clusters)[order(unique(CD8@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"expansion",sep = ""),showWarnings = FALSE)
rm(sha_ent)
cluster_name <- gsub("C","",cluster_num)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Gini_coefficient_index.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Gini_coef_AJ(object = CD8, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                               group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                               column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

write.table(df2, paste(savedir,"expansion/Gini_coef_clonality_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df2$sample <- row.names(df2)
df2_melted <- melt(df2, id = "sample")
colnames(df2_melted) <- c("Sample","Cluster","Gini_coef")

column_name = c("id","Ag","Age","Age_number","gender","Run")
split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name

if("condition" %in% colnames(df2_melted_split) & "vaccine" %in% colnames(df2_melted_split)){
  df2_melted_split[grep("DMSO",df2_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Gini_coef)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Gini coefficient Clonality")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"expansion/vaccine_Gini_coef_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### Celltypes
rm(df2)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/clonality/all_cells/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
sample_id <- CD8@meta.data$orig.ident %>% unique()
CT_name <- unique(CD8@meta.data$celltypes)
cluster_num <- CT_name[order(CT_name)]
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"expansion",sep = ""),showWarnings = FALSE)
rm(sha_ent)
cluster_name <- cluster_num

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Gini_coefficient_index.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Gini_coef_AJ(object = CD8, savedir = savedir, clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                               group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                               column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

write.table(df2, paste(savedir,"expansion/Gini_coef_clonality_sample_all_celltypes.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df2$sample <- row.names(df2)
df2_melted <- melt(df2, id = "sample")
colnames(df2_melted) <- c("Sample","Cluster","Gini_coef")

column_name = c("id","Ag","Age","Age_number","gender","Run")
split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name

if("condition" %in% colnames(df2_melted_split) & "vaccine" %in% colnames(df2_melted_split)){
  df2_melted_split[grep("DMSO",df2_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Gini_coef)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Gini coefficient Clonality")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"expansion/vaccine_Gini_coef_celltypes_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

#### CD8 memory cells
rm(df2)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/clonality/mem_cells/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
sample_id <- CD8_mem@meta.data$orig.ident %>% unique()
cluster_num <- paste("C",unique(CD8_mem@meta.data$seurat_clusters)[order(unique(CD8_mem@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"expansion",sep = ""),showWarnings = FALSE)
rm(sha_ent)
cluster_name <- gsub("C","",cluster_num)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Gini_coefficient_index.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Gini_coef_AJ(object = CD8_mem, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                               group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                               column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

write.table(df2, paste(savedir,"expansion/Gini_coef_clonality_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df2 <- read.table(paste(savedir,"expansion/Gini_coef_clonality_sample_all_cluster.txt",sep = ""), header = TRUE)
df2$sample <- row.names(df2)
df2_melted <- melt(df2, id = "sample")
colnames(df2_melted) <- c("Sample","Cluster","Gini_coef")

column_name = c("id","Ag","Age","Age_number","gender","Run")
split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name

if("condition" %in% colnames(df2_melted_split) & "vaccine" %in% colnames(df2_melted_split)){
  df2_melted_split[grep("DMSO",df2_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Gini_coef)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Gini coefficient Clonality")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"expansion/vaccine_Gini_coef_clusters_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### Celltypes
rm(df2)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/clonality/mem_Ag/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
sample_id <- CD8_mem@meta.data$orig.ident %>% unique()
CT_name <- unique(CD8_mem@meta.data$celltypes)
cluster_num <- CT_name[order(CT_name)]
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(sample_id)))
rownames(df2) <- sample_id
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"expansion",sep = ""),showWarnings = FALSE)
rm(sha_ent)
cluster_name <- cluster_num

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Gini_coefficient_index.R")
for (i in 1:length(sample_id)) {
  try({df2[i,] <- Gini_coef_AJ(object = CD8_mem, savedir = savedir, clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                               group_col = "orig.ident", group_val = sample_id[i], split_col = "orig.ident",
                               column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

write.table(df2, paste(savedir,"expansion/Gini_coef_clonality_sample_all_celltypes.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

df2$sample <- row.names(df2)
df2_melted <- melt(df2, id = "sample")
colnames(df2_melted) <- c("Sample","Cluster","Gini_coef")

column_name = c("id","Ag","Age","Age_number","gender","Run")
split_to_dataframe <- function(x) {
  split_elements <- strsplit(x, "_")[[1]]
  data.frame(t(split_elements))
}

df2_melted_split <- do.call(rbind, lapply(df2_melted[,"Sample"], split_to_dataframe))
colnames(df2_melted_split) <- column_name

if("condition" %in% colnames(df2_melted_split) & "vaccine" %in% colnames(df2_melted_split)){
  df2_melted_split[grep("DMSO",df2_melted_split$condition),"vaccine"] <- "DMSO"
}

for (i in 1:length(column_name)) {
  df2_melted[,column_name[i]] <- df2_melted_split[,column_name[i]]
}

p <- ggplot(na.omit(df2_melted), aes(x=Cluster,y=Gini_coef)) + 
  geom_boxplot(aes(fill = Age)) + geom_point(position=position_dodge(width=0.75),aes(group=Age)) +
  ggtitle(paste("Gini coefficient Clonality")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

pdf(paste(savedir,"expansion/vaccine_Gini_coef_celltypes_noNA.pdf",sep = ""), width = 21, height = 6)
p
dev.off()

### Memory Ag
#### CD8 memory cells
rm(df2)
Ag_sampleid <- grep(",",paste(CD8_mem@meta.data$Ag_range_10,CD8_mem@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
cluster_num <- paste("C",unique(CD8_mem@meta.data$seurat_clusters)[order(unique(CD8_mem@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",cluster_num)
rm(sha_ent)
cluster_name <- unique(CD8_mem@meta.data$seurat_clusters)[order(unique(CD8_mem@meta.data$seurat_clusters))]

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Simpson_Diversity.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Simpson_AJ(object = CD8_mem, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                             group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                             column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}
write.table(df2, paste(savedir,"clusters/Simpson_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")


### memory Ag
rm(df2)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/clonality/mem_Ag/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
CD8_mem@meta.data$Ag_sampleid <- paste(CD8_mem@meta.data$Ag_range_10,CD8_mem@meta.data$orig.ident,sep = "_")

Ag_sampleid <- grep(",",paste(CD8_mem@meta.data$Ag_range_10,CD8_mem@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
cluster_num <- paste("C",unique(CD8_mem@meta.data$seurat_clusters)[order(unique(CD8_mem@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"expansion",sep = ""),showWarnings = FALSE)
rm(sha_ent)
cluster_name <- gsub("C","",cluster_num)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Gini_coefficient_index.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Gini_coef_AJ(object = CD8_mem, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                               group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                               column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

write.table(df2, paste(savedir,"expansion/Gini_coef_clonality_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")


### Celltypes
rm(df2)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/clonality/mem_Ag/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
CD8_mem@meta.data$Ag_sampleid <- paste(CD8_mem@meta.data$Ag_range_10,CD8_mem@meta.data$orig.ident,sep = "_")

Ag_sampleid <- grep(",",paste(CD8_mem@meta.data$Ag_range_10,CD8_mem@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
CT_name <- unique(CD8_mem@meta.data$celltypes)
cluster_num <- CT_name[order(CT_name)]

df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"expansion",sep = ""),showWarnings = FALSE)
rm(sha_ent)
cluster_name <- cluster_num

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Gini_coefficient_index.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Gini_coef_AJ(object = CD8_mem, savedir = savedir, clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                               group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                               column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

write.table(df2, paste(savedir,"expansion/Gini_coef_clonality_sample_all_celltypes.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

### Naive Ag
rm(df2)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/clonality/naive_Ag/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
CD8_naive@meta.data$Ag_sampleid <- paste(CD8_naive@meta.data$Ag_range_10,CD8_naive@meta.data$orig.ident,sep = "_")

Ag_sampleid <- grep(",",paste(CD8_naive@meta.data$Ag_range_10,CD8_naive@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
cluster_num <- paste("C",unique(CD8_naive@meta.data$seurat_clusters)[order(unique(CD8_naive@meta.data$seurat_clusters))],sep = "")
df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"expansion",sep = ""),showWarnings = FALSE)
rm(sha_ent)
cluster_name <- gsub("C","",cluster_num)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Gini_coefficient_index.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Gini_coef_AJ(object = CD8_naive, savedir = savedir, clus_col = "seurat_clusters", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                               group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                               column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

write.table(df2, paste(savedir,"expansion/Gini_coef_clonality_sample_all_cluster.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")


### Celltypes
rm(df2)
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/clonality/naive_Ag/"
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
CD8_naive@meta.data$Ag_sampleid <- paste(CD8_naive@meta.data$Ag_range_10,CD8_naive@meta.data$orig.ident,sep = "_")

Ag_sampleid <- grep(",",paste(CD8_naive@meta.data$Ag_range_10,CD8_naive@meta.data$orig.ident,sep = "_"),invert=TRUE, value=TRUE) %>% unique()
CT_name <- unique(CD8_naive@meta.data$celltypes)
cluster_num <- CT_name[order(CT_name)]

df2 <- data.frame(matrix(ncol = length(cluster_num)+1, nrow = length(Ag_sampleid)))
rownames(df2) <- Ag_sampleid
colnames(df2) <- c("all",cluster_num)
dir.create(paste(savedir,"expansion",sep = ""),showWarnings = FALSE)
rm(sha_ent)
cluster_name <- cluster_num

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/Gini_coefficient_index.R")
for (i in 1:length(Ag_sampleid)) {
  try({df2[i,] <- Gini_coef_AJ(object = CD8_naive, savedir = savedir, clus_col = "celltypes", TCR_col = "TRA1_or_TRA2_TRB1_no_TRB2",
                               group_col = "Ag_sampleid", group_val = Ag_sampleid[i], split_col = "orig.ident",
                               column_name = c("id","Ag","Age","Age_number","gender","Run"), total_clusters = cluster_name)
  })
}

write.table(df2, paste(savedir,"expansion/Gini_coef_clonality_sample_all_celltypes.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

#### Performing the scdifferential for the same clones different Ag ####
### Clone1 
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/"
clone <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/clone_1/clone1", header=TRUE,sep="\t")
CD8@meta.data$same_clone_Ag <- clone[match(rownames(CD8@meta.data), clone$barcodes),"Ag"]
CD8@meta.data$same_clone_Ag <- gsub(", ","-",CD8@meta.data$same_clone_Ag)
clone_celltypes <- table(paste(CD8@meta.data$celltypes,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_celltypes_req <- clone_celltypes[grep("NA$",clone_celltypes$Var1, invert=TRUE),]
write.table(clone_celltypes_req, paste(savedir,"clone_1/clone_celltypes_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

clone_cluster <- table(paste(CD8@meta.data$seurat_clusters,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_cluster_req <- clone_cluster[grep("NA$",clone_cluster$Var1, invert=TRUE),]
write.table(clone_cluster_req, paste(savedir,"clone_1/clone_clus_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

### scdifferential
VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "VZV_IE62", ident.2 = "EBV_EBNA1", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_1/VZV_IE62_vs_EBV_EBNA1.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_cells <- rownames(CD8@meta.data[grep("^VZV_IE62$", CD8@meta.data$same_clone_Ag),])
VZV_IE62_EBV_EBNA1_cells <- rownames(CD8@meta.data[grep("-", CD8@meta.data$same_clone_Ag),])

VZV_IE62_vs_VZV_IE62_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = VZV_IE62_cells, ident.2 = VZV_IE62_EBV_EBNA1_cells, logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_VZV_IE62_EBV_EBNA1, paste(savedir,"clone_1/VZV_IE62_vs_VZV_IE62_EBV_EBNA1.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/clone_UMAP.R")
CD8@meta.data$same_clone_Ag2 <- CD8@meta.data$same_clone_Ag
CD8@meta.data[grep(",",CD8@meta.data$same_clone_Ag2),"same_clone_Ag2"] <- "VZV_IE62-EBV_EBNA1"
clone_UMAP_AJ(object=CD8, reductions = "wnn.umap",umap1="wnnUMAP_1", umap2="wnnUMAP_2", highlight="same_clone_Ag2", savedir=paste(savedir,"clone_1/",sep = ""), colour = c("red","green","blue"))

### Clone4
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/"
clone <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/clone_4/clone4", header=TRUE,sep="\t")
CD8@meta.data$same_clone_Ag <- clone[match(rownames(CD8@meta.data), clone$barcodes),"Ag"]
CD8@meta.data$same_clone_Ag <- gsub(", ","-",CD8@meta.data$same_clone_Ag)
clone_celltypes <- table(paste(CD8@meta.data$celltypes,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_celltypes_req <- clone_celltypes[grep("NA$",clone_celltypes$Var1, invert=TRUE),]
write.table(clone_celltypes_req, paste(savedir,"clone_4/clone_celltypes_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

clone_cluster <- table(paste(CD8@meta.data$seurat_clusters,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_cluster_req <- clone_cluster[grep("NA$",clone_cluster$Var1, invert=TRUE),]
write.table(clone_cluster_req, paste(savedir,"clone_4/clone_clus_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

### scdifferential
VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "EBV_EBNA1", ident.2 = "VZV_IE62", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_4/EBV_EBNA1_vs_VZV_IE62.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/clone_UMAP.R")
clone_UMAP_AJ(object=CD8, reductions = "wnn.umap",umap1="wnnUMAP_1", umap2="wnnUMAP_2", highlight="same_clone_Ag", savedir=paste(savedir,"clone_4/",sep = ""), colour = c("red","green"))

### clone7
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/"
clone <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/clone_7/clone7", header=TRUE,sep="\t")
CD8@meta.data$same_clone_Ag <- clone[match(rownames(CD8@meta.data), clone$barcodes),"Ag"]
CD8@meta.data$same_clone_Ag <- gsub(", ","-",CD8@meta.data$same_clone_Ag)
clone_celltypes <- table(paste(CD8@meta.data$celltypes,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_celltypes_req <- clone_celltypes[grep("NA$",clone_celltypes$Var1, invert=TRUE),]
write.table(clone_celltypes_req, paste(savedir,"clone_7/clone_celltypes_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

clone_cluster <- table(paste(CD8@meta.data$seurat_clusters,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_cluster_req <- clone_cluster[grep("NA$",clone_cluster$Var1, invert=TRUE),]
write.table(clone_cluster_req, paste(savedir,"clone_7/clone_clus_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

### scdifferential
VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "EBV_EBNA1", ident.2 = "EBV_EBNA3C", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_7/EBV_EBNA1_vs_EBV_EBNA3C.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "EBV_EBNA1", ident.2 = "VZV_IE62", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_7/EBV_EBNA1_vs_VZV_IE62.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "EBV_EBNA3C", ident.2 = "VZV_IE62", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_7/EBV_EBNA3C_vs_VZV_IE62.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/clone_UMAP.R")
clone_UMAP_AJ(object=CD8, reductions = "wnn.umap",umap1="wnnUMAP_1", umap2="wnnUMAP_2", highlight="same_clone_Ag", savedir=paste(savedir,"clone_7/",sep = ""), colour = c("red","green","blue"))

### Clone9
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/"
clone <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/clone_9/clone9", header=TRUE,sep="\t")
CD8@meta.data$same_clone_Ag <- clone[match(rownames(CD8@meta.data), clone$barcodes),"Ag"]
CD8@meta.data$same_clone_Ag <- gsub(", ","-",CD8@meta.data$same_clone_Ag)
clone_celltypes <- table(paste(CD8@meta.data$celltypes,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_celltypes_req <- clone_celltypes[grep("NA$",clone_celltypes$Var1, invert=TRUE),]
write.table(clone_celltypes_req, paste(savedir,"clone_9/clone_celltypes_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

clone_cluster <- table(paste(CD8@meta.data$seurat_clusters,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_cluster_req <- clone_cluster[grep("NA$",clone_cluster$Var1, invert=TRUE),]
write.table(clone_cluster_req, paste(savedir,"clone_9/clone_clus_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

### scdifferential
VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "EBV_BRLF1", ident.2 = "VZV_IE62", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_9/EBV_BRLF1_vs_VZV_IE62.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

EBV_BRLF1_cells <- rownames(CD8@meta.data[grep("^EBV_BRLF1$", CD8@meta.data$same_clone_Ag),])
VZV_IE62_cells <- rownames(CD8@meta.data[grep("^VZV_IE62$", CD8@meta.data$same_clone_Ag),])
EBV_BRLF1_VZV_IE62_cells <- rownames(CD8@meta.data[grep("-", CD8@meta.data$same_clone_Ag),])

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = EBV_BRLF1_cells, ident.2 = EBV_BRLF1_VZV_IE62_cells, group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_9/EBV_BRLF1_vs_EBV_BRLF1_VZV_IE62_cells.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_vs_VZV_IE62_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = VZV_IE62_cells, ident.2 = EBV_BRLF1_VZV_IE62_cells, logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_VZV_IE62_EBV_EBNA1, paste(savedir,"clone_9/VZV_IE62_vs_EBV_BRLF1_VZV_IE62.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/clone_UMAP.R")
CD8@meta.data$same_clone_Ag2 <- CD8@meta.data$same_clone_Ag
CD8@meta.data[grep(",",CD8@meta.data$same_clone_Ag2),"same_clone_Ag2"] <- "EBV_BRLF1-VZV_IE62"
clone_UMAP_AJ(object=CD8, reductions = "wnn.umap",umap1="wnnUMAP_1", umap2="wnnUMAP_2", highlight="same_clone_Ag2", savedir=paste(savedir,"clone_9/",sep = ""), colour = c("red","green","blue"))


### clone10
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/"
clone <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/clone_10/clone10", header=TRUE,sep="\t")
CD8@meta.data$same_clone_Ag <- clone[match(rownames(CD8@meta.data), clone$barcodes),"Ag"]
CD8@meta.data$same_clone_Ag <- gsub(", ","-",CD8@meta.data$same_clone_Ag)
clone_celltypes <- table(paste(CD8@meta.data$celltypes,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_celltypes_req <- clone_celltypes[grep("NA$",clone_celltypes$Var1, invert=TRUE),]
write.table(clone_celltypes_req, paste(savedir,"clone_10/clone_celltypes_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

clone_cluster <- table(paste(CD8@meta.data$seurat_clusters,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_cluster_req <- clone_cluster[grep("NA$",clone_cluster$Var1, invert=TRUE),]
write.table(clone_cluster_req, paste(savedir,"clone_10/clone_clus_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

### scdifferential
VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "VZV_IE62", ident.2 = "EBV_EBNA1", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_10/VZV_IE62_vs_EBV_EBNA1.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "VZV_IE62", ident.2 = "EBV_BMLF1", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_10/VZV_IE62_vs_EBV_BMLF1.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "EBV_EBNA1", ident.2 = "EBV_BMLF1", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_10/EBV_EBNA1_vs_EBV_BMLF1.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/clone_UMAP.R")
clone_UMAP_AJ(object=CD8, reductions = "wnn.umap",umap1="wnnUMAP_1", umap2="wnnUMAP_2", highlight="same_clone_Ag", savedir=paste(savedir,"clone_10/",sep = ""), colour = c("red","green","blue"))


### Clone 13
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/"
clone <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/clone_13/clone13", header=TRUE,sep="\t")
CD8@meta.data$same_clone_Ag <- clone[match(rownames(CD8@meta.data), clone$barcodes),"Ag"]
CD8@meta.data$same_clone_Ag <- gsub(", ","-",CD8@meta.data$same_clone_Ag)
clone_celltypes <- table(paste(CD8@meta.data$celltypes,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_celltypes_req <- clone_celltypes[grep("NA$",clone_celltypes$Var1, invert=TRUE),]
write.table(clone_celltypes_req, paste(savedir,"clone_13/clone_celltypes_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

clone_cluster <- table(paste(CD8@meta.data$seurat_clusters,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_cluster_req <- clone_cluster[grep("NA$",clone_cluster$Var1, invert=TRUE),]
write.table(clone_cluster_req, paste(savedir,"clone_13/clone_clus_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

### scdifferential
VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "EBV_EBNA3C", ident.2 = "VZV_IE62", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_13/EBV_EBNA3C_vs_VZV_IE62.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "EBV_EBNA3C", ident.2 = "EBV_EBNA1", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_13/EBV_EBNA3C_vs_EBV_EBNA1.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "VZV_IE62", ident.2 = "EBV_EBNA1", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_13/VZV_IE62_vs_EBV_EBNA1.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

EBV_EBNA3C_cells <- rownames(CD8@meta.data[grep("^EBV_EBNA3C$", CD8@meta.data$same_clone_Ag),])
VZV_IE62_cells <- rownames(CD8@meta.data[grep("^VZV_IE62$", CD8@meta.data$same_clone_Ag),])
EBV_EBNA3C_EBV_BRLF1_cells <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/clone_13/EBV_EBNA3C_EBV_BRLF1_cellbarcodes.txt", header = FALSE)[,1]
EBV_EBNA3C_VZV_IE62_cells <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/clone_13/EBV_EBNA3C_VZV_IE62_cellbarodes.txt", header = FALSE)[,1]

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = EBV_EBNA3C_cells, ident.2 = EBV_EBNA3C_EBV_BRLF1_cells, logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_13/EBV_EBNA3C_vs_EBV_EBNA3C_EBV_BRLF1.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = EBV_EBNA3C_cells, ident.2 = EBV_EBNA3C_VZV_IE62_cells, logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_13/EBV_EBNA3C_vs_EBV_EBNA3C_VZV_IE62.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = EBV_EBNA3C_EBV_BRLF1_cells, ident.2 = EBV_EBNA3C_VZV_IE62_cells, logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_13/EBV_EBNA3C_EBV_BRLF1_vs_EBV_EBNA3C_VZV_IE62.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = VZV_IE62_cells, ident.2 = EBV_EBNA3C_VZV_IE62_cells, logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_13/VZV_IE62_vs_EBV_EBNA3C_VZV_IE62.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/clone_UMAP.R")
CD8@meta.data$same_clone_Ag2 <- CD8@meta.data$same_clone_Ag
CD8@meta.data[grep("EBV_EBNA3C-EBV_BRLF1|EBV_BRLF1-EBV_EBNA3C",CD8@meta.data$same_clone_Ag2),"same_clone_Ag2"] <- "EBV_EBNA3C-EBV_BRLF1"
CD8@meta.data[grep("EBV_EBNA3C-VZV_IE62|VZV_IE62-EBV_EBNA3C",CD8@meta.data$same_clone_Ag2),"same_clone_Ag2"] <- "EBV_EBNA3C-VZV_IE62"
clone_UMAP_AJ(object=CD8, reductions = "wnn.umap",umap1="wnnUMAP_1", umap2="wnnUMAP_2", highlight="same_clone_Ag2", savedir=paste(savedir,"clone_13/",sep = ""), colour = c("red","green","blue","magenta","black"))


### Clone 18
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/"
clone <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/clone_18/clone18", header=TRUE,sep="\t")
CD8@meta.data$same_clone_Ag <- clone[match(rownames(CD8@meta.data), clone$barcodes),"Ag"]
CD8@meta.data$same_clone_Ag <- gsub(", ","-",CD8@meta.data$same_clone_Ag)
clone_celltypes <- table(paste(CD8@meta.data$celltypes,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_celltypes_req <- clone_celltypes[grep("NA$",clone_celltypes$Var1, invert=TRUE),]
write.table(clone_celltypes_req, paste(savedir,"clone_18/clone_celltypes_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

clone_cluster <- table(paste(CD8@meta.data$seurat_clusters,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_cluster_req <- clone_cluster[grep("NA$",clone_cluster$Var1, invert=TRUE),]
write.table(clone_cluster_req, paste(savedir,"clone_18/clone_clus_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

### scdifferential
VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "EBV_EBNA1", ident.2 = "EBV_EBNA3C", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_18/EBV_EBNA1_vs_EBV_EBNA3C.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "EBV_EBNA1", ident.2 = "VZV_IE62", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_18/EBV_EBNA1_vs_VZV_IE62.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "EBV_EBNA3C", ident.2 = "VZV_IE62", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_18/EBV_EBNA3C_vs_VZV_IE62.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/clone_UMAP.R")
clone_UMAP_AJ(object=CD8, reductions = "wnn.umap",umap1="wnnUMAP_1", umap2="wnnUMAP_2", highlight="same_clone_Ag", savedir=paste(savedir,"clone_18/",sep = ""), colour = c("red","green","blue"))


### Clone 21
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/"
clone <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scdifferential/diff_Ag_same_clones/clone_21/clone21", header=TRUE,sep="\t")
CD8@meta.data$same_clone_Ag <- clone[match(rownames(CD8@meta.data), clone$barcodes),"Ag"]
CD8@meta.data$same_clone_Ag <- gsub(", ","-",CD8@meta.data$same_clone_Ag)
clone_celltypes <- table(paste(CD8@meta.data$celltypes,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_celltypes_req <- clone_celltypes[grep("NA$",clone_celltypes$Var1, invert=TRUE),]
write.table(clone_celltypes_req, paste(savedir,"clone_21/clone_celltypes_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

clone_cluster <- table(paste(CD8@meta.data$seurat_clusters,CD8@meta.data$same_clone_Ag)) %>% as.data.frame()
clone_cluster_req <- clone_cluster[grep("NA$",clone_cluster$Var1, invert=TRUE),]
write.table(clone_cluster_req, paste(savedir,"clone_21/clone_clus_req.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

### scdifferential
VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "EBV_EBNA1", ident.2 = "EBV_EBNA3C", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_21/EBV_EBNA1_vs_EBV_EBNA3C.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "EBV_EBNA1", ident.2 = "VZV_IE62", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_21/EBV_EBNA1_vs_VZV_IE62.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

VZV_IE62_vs_EBV_EBNA1 <- FindMarkers(CD8, ident.1 = "EBV_EBNA3C", ident.2 = "VZV_IE62", group.by="same_clone_Ag", logfc.threshold = 0, min.pct = 0)
write.table(VZV_IE62_vs_EBV_EBNA1, paste(savedir,"clone_21/EBV_EBNA3C_vs_VZV_IE62.txt",sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scTCRseq/clone_UMAP.R")
clone_UMAP_AJ(object=CD8, reductions = "wnn.umap",umap1="wnnUMAP_1", umap2="wnnUMAP_2", highlight="same_clone_Ag", savedir=paste(savedir,"clone_21/",sep = ""), colour = c("red","green","blue"))


#### TRAV MAIT
CD8_all <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/saveRDS_obj/CD8_modality_integrate_cluster_0.6.RDS")
TRAV1_2_TRAJ12_cells <- CD8_all@meta.data[grep("TRAV1-2__TRAJ12",CD8_all@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc),] %>% rownames()
TRAV1_2__TRAJ20_cells <- CD8_all@meta.data[grep("TRAV1-2__TRAJ20_TRAC",CD8_all@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc),] %>% rownames()
TRAV1_2__TRAJ33_cells <- CD8_all@meta.data[grep("TRAV1-2__TRAJ33_TRAC",CD8_all@meta.data$TRA1_TRA2_TRB1_TRB2_vdjc),] %>% rownames()

p1 <- DimPlot(CD8_all,cells.highlight = TRAV1_TRAJ12_cells, reduction = "wnn.umap", 
             label = TRUE, cols.highlight = "deeppink2",
             cols = "gray92") + 
  ggtitle("TRAV1_TRAJ12_cells") + 
  NoLegend() + 
  theme(plot.title = element_text(hjust = 0.2))

p2 <- DimPlot(CD8_all,cells.highlight = TRAV1_2__TRAJ20_cells, reduction = "wnn.umap", 
              label = TRUE, cols.highlight = "deeppink2",
              cols = "gray92") + 
  ggtitle("TRAV1_2__TRAJ20_cells") + 
  NoLegend() + 
  theme(plot.title = element_text(hjust = 0.2))

p3 <- DimPlot(CD8_all,cells.highlight = TRAV1_2__TRAJ33_cells, reduction = "wnn.umap", 
              label = TRUE, cols.highlight = "deeppink2",
              cols = "gray92") + 
  ggtitle("TRAV1_2__TRAJ33_cells") + 
  NoLegend() + 
  theme(plot.title = element_text(hjust = 0.2))

pdf("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/VDJ/clone_no_TRB2_with_TRB1_with_TRA1_or_TRA2/UMAP/TRAV1_TRAJ12_cells.pdf", width = 15, height = 5)
p1+p2+p3
dev.off()

### Making an Upset plot ####
diff_files <- list.files("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/Range_10_rem_TRB2_rem_batch_effect/pseudobulk/",
                         recursive = TRUE, full.names = TRUE, pattern = "O_VS_Y_sig_up.txt")

Ags <- gsub("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/pseudobulk//clus_all_O.","",diff_files) %>% gsub("/Differential/O_VS_Y_sig_up.txt","",.) %>% gsub("_vs_.*.","",.)

for (i in 1:length(diff_files)) {
  spec_Ag <- rownames(read.table(diff_files[i]))
  assign(Ags[i],spec_Ag)
}

### Please check how many Antigens are there then move forward in Ags based on that move ahead and change the numbering based on that
print(Ags)
library(data.table)
BMLF1 <- get(Ags[1])
# BMRF1 <- get(Ags[2])
BRLF1 <- get(Ags[2])
EBNA1 <- get(Ags[3])
EBNA3C <- get(Ags[4])
IE62 <- get(Ags[5])
LMP1 <- get(Ags[6])
LMP2 <- get(Ags[7])

n <- max(length(BMLF1), length(BRLF1), length(EBNA1), 
         length(EBNA3C), length(IE62), length(LMP1), length(LMP2))

length(BMLF1) = n                      
# length(BMRF1) = n
length(BRLF1) = n
length(EBNA1) = n
length(EBNA3C) = n
length(IE62) = n
length(LMP1) = n
length(LMP2) = n 

df = as.data.frame(cbind(BMLF1, BRLF1, EBNA1, EBNA3C, IE62, LMP1, LMP2))

df$V1 <- NULL

fld <- fromList(as.list(df))

# savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/pseudobulk/"
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
write.table(as.data.frame(final), paste(savedir,"Upregulated_gene_number_2.txt",sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")

### Downregulated Genes
diff_files <- list.files("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/pseudobulk/",
                         recursive = TRUE, full.names = TRUE, pattern = "O_VS_Y_sig_down.txt")

Ags <- gsub("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/pseudobulk//clus_all_O.","",diff_files) %>% gsub("/Differential/O_VS_Y_sig_down.txt","",.) %>% gsub("_vs_.*.","",.)

for (i in 1:length(diff_files)) {
  spec_Ag <- rownames(read.table(diff_files[i]))
  assign(Ags[i],spec_Ag)
}

library(data.table)
library(UpSetR)
BMLF1 <- get(Ags[1])
BMRF1 <- get(Ags[2])
BRLF1 <- get(Ags[3])
EBNA1 <- get(Ags[4])
EBNA3C <- get(Ags[5])
IE62 <- get(Ags[6])
LMP1 <- get(Ags[7])
LMP2 <- get(Ags[8])

n <- max(length(BMLF1), length(BMRF1), length(BRLF1), length(EBNA1),
         length(EBNA3C), length(IE62), length(LMP1), length(LMP2))

length(BMLF1) = n
length(BMRF1) = n
length(BRLF1) = n
length(EBNA1) = n
length(EBNA3C) = n
length(IE62) = n
length(LMP1) = n
length(LMP2) = n 

df = as.data.frame(cbind(BMLF1,BMRF1, BRLF1, EBNA1, EBNA3C, IE62, LMP1, LMP2))

df$V1 <- NULL

fld <- fromList(as.list(df))

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/pseudobulk/"
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

#### Clones Kuchroo Paper Fig 2E ####
required_col <- select(CD8_mem_spec_Ag@meta.data, "TRA1_or_TRA2_TRB1_no_TRB2", "AgeAntigen2", "celltypes", "seurat_clusters")
required_col_noNA <- required_col[!is.na(required_col$TRA1_or_TRA2_TRB1_no_TRB2),]
Age_Antigens <- unique(required_col_noNA$AgeAntigen2)

write.table(required_col_noNA, "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scTCR/celltypes/Khuchroo_paper/mem/clones_required_col_noNA.txt", sep="\t", col.names=T, row.names=T, quote=F)

### Running it in RStudio
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scTCR/celltypes/Khuchroo_paper/mem/"
required_col_noNA <- read.table(paste(savedir,"clones_required_col_noNA.txt",sep = ""), header = TRUE)
Age_Antigens <- unique(required_col_noNA$AgeAntigen2)[order(unique(required_col_noNA$AgeAntigen2))]

for (k in 1:length(Age_Antigens)) {
  try({
  required_col_noNA_Ag <- required_col_noNA[grep(Age_Antigens[k], required_col_noNA$AgeAntigen2),]
  required_col_noNA_Ag_table <- table(required_col_noNA_Ag$TRA1_or_TRA2_TRB1_no_TRB2, required_col_noNA_Ag$celltypes)
  rownames(required_col_noNA_Ag_table) <- NULL
  
  index <- vector()
  for (i in 1:nrow(required_col_noNA_Ag_table)) {
    ### removing those column that are not sharing with other clusters
    clone_count <- as.data.frame(table(required_col_noNA_Ag_table[i,]))
    total_col <- ncol(required_col_noNA_Ag_table)
    ### If total_col - 1 ==0 column has 0 we remove that column
    number = total_col-1
    if (clone_count[1,1] == 0) {
      if (clone_count[1,2] == number) {
        print(paste(i,"not Sharing",sep = " "))
        index[i] = NA
      }
      else{
        print(paste(i,"is Sharing",sep = " "))
        index[i] <- i
      }
    }
    else{
      index[i] <- i
    }
  }
  
  required_col_noNA_Ag_table_sharing <- required_col_noNA_Ag_table[index[!is.na(index)],]
  
  ### Making an empty dataframe
  Ag_sharing_percent <- as.data.frame(matrix(nrow=nrow(required_col_noNA_Ag_table_sharing), ncol=ncol(required_col_noNA_Ag_table_sharing)))
  colnames(Ag_sharing_percent) <- colnames(required_col_noNA_Ag_table_sharing)
  
  for (j in 1:nrow(required_col_noNA_Ag_table_sharing)) {
    Ag_sharing_percent[j,] <- (required_col_noNA_Ag_table_sharing[j,]/rowSums(required_col_noNA_Ag_table_sharing)[j])*100
  }
  
  assign(paste(Age_Antigens[k],"_percent",sep = ""),Ag_sharing_percent)
  })
}


Ag_vec <- ls(pattern = "_percent")
Ag_vec <- Ag_vec[-1]

### Please find the correct order of the celltypes
celltypes_order <- colnames(get(Ag_vec[2]))

### Making the equal columns with the same order
for (i in 1:length(Ag_vec)) {
  if (ncol(get(Ag_vec[i])) < length(unique(required_col_noNA$celltypes))) {
    df <- get(Ag_vec[i])
    missing_index <- grep(paste(colnames(df),collapse = "|"),celltypes_order,invert = TRUE,value = TRUE)
    df[,missing_index] <- 0
    
    ## ordering it
    df <- df[,match(celltypes_order,colnames(df))]
    assign(Ag_vec[i],df)
  }
}

## Combining all
all_combine <- rbind(get(Ag_vec[1]),get(Ag_vec[2]),get(Ag_vec[3]),get(Ag_vec[4]),get(Ag_vec[5]),get(Ag_vec[6]),get(Ag_vec[7]),get(Ag_vec[8]),get(Ag_vec[9]),get(Ag_vec[10]),
      get(Ag_vec[11]),get(Ag_vec[12]),get(Ag_vec[13]),get(Ag_vec[14]),get(Ag_vec[15]),get(Ag_vec[16]),get(Ag_vec[17]))

Ag_name <- gsub("_percent","",Ag_vec)

heatmap_rownames <- list()
for (i in 1:17) {
  heatmap_rownames[[i]] <- rep(Ag_name[i],nrow(get(Ag_vec[i])))
}

column_name <- rowAnnotation(Ag_group = unlist(heatmap_rownames))

library(ComplexHeatmap)
p <- Heatmap(all_combine,
             name = "Percentage",
             column_title = "Celltypes",
             row_title = "Clones",
             row_names_gp = gpar(fontsize = 12),
             column_names_gp = gpar(fontsize = 12),
        left_annotation = column_name,
             cluster_columns = FALSE,
             cluster_rows = FALSE,
        col=colorRamp2(c(0, 100), c("white", "maroon"))
)

pdf(paste(savedir,"clones_heatmap_mem_celltypes.pdf",sep = ""), width = 7, height = 9)
p
dev.off()

all_combine$Age_Antigen <- unlist(heatmap_rownames)
write.table(all_combine, paste(savedir,"clones_heatmap_mem_celltypes.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

###  Clusters ####
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scTCR/celltypes/Khuchroo_paper/mem/"
required_col_noNA <- read.table(paste(savedir,"clones_required_col_noNA.txt",sep = ""), header = TRUE)
Age_Antigens <- unique(required_col_noNA$AgeAntigen2)[order(unique(required_col_noNA$AgeAntigen2))]

for (k in 1:length(Age_Antigens)) {
  try({
    required_col_noNA_Ag <- required_col_noNA[grep(Age_Antigens[k], required_col_noNA$AgeAntigen2),]
    required_col_noNA_Ag_table <- table(required_col_noNA_Ag$TRA1_or_TRA2_TRB1_no_TRB2, required_col_noNA_Ag$seurat_clusters)
    rownames(required_col_noNA_Ag_table) <- NULL
    
    index <- vector()
    for (i in 1:nrow(required_col_noNA_Ag_table)) {
      ### removing those column that are not sharing with other clusters
      clone_count <- as.data.frame(table(required_col_noNA_Ag_table[i,]))
      total_col <- ncol(required_col_noNA_Ag_table)
      ### If total_col - 1 ==0 column has 0 we remove that column
      number = total_col-1
      if (clone_count[1,1] == 0) {
        if (clone_count[1,2] == number) {
          print(paste(i,"not Sharing",sep = " "))
          index[i] = NA
        }
        else{
          print(paste(i,"is Sharing",sep = " "))
          index[i] <- i
        }
      }
      else{
        index[i] <- i
      }
    }
    
    required_col_noNA_Ag_table_sharing <- required_col_noNA_Ag_table[index[!is.na(index)],]
    
    ### Making an empty dataframe
    Ag_sharing_percent <- as.data.frame(matrix(nrow=nrow(required_col_noNA_Ag_table_sharing), ncol=ncol(required_col_noNA_Ag_table_sharing)))
    colnames(Ag_sharing_percent) <- colnames(required_col_noNA_Ag_table_sharing)
    
    for (j in 1:nrow(required_col_noNA_Ag_table_sharing)) {
      Ag_sharing_percent[j,] <- (required_col_noNA_Ag_table_sharing[j,]/rowSums(required_col_noNA_Ag_table_sharing)[j])*100
    }
    
    assign(paste(Age_Antigens[k],"_percent",sep = ""),Ag_sharing_percent)
  })
}


Ag_vec <- ls(pattern = "_percent")
Ag_vec <- Ag_vec[-1]

### Please find the correct order of the celltypes
celltypes_order <- as.data.frame(table(required_col_noNA$seurat_clusters))[,1]

### Making the equal columns with the same order
for (i in 1:length(Ag_vec)) {
  if (ncol(get(Ag_vec[i])) < length(unique(required_col_noNA$seurat_clusters))) {
    df <- get(Ag_vec[i])
    missing_index <- grep(paste("^",colnames(df),"$",collapse = "|",sep = ""),celltypes_order,invert = TRUE,value = TRUE)
    df[,missing_index] <- 0
    
    ## ordering it
    df <- df[,match(celltypes_order,colnames(df))]
    assign(Ag_vec[i],df)
  }
}

## Combining all
all_combine <- rbind(get(Ag_vec[1]),get(Ag_vec[2]),get(Ag_vec[3]),get(Ag_vec[4]),get(Ag_vec[5]),get(Ag_vec[6]),get(Ag_vec[7]),get(Ag_vec[8]),get(Ag_vec[9]),get(Ag_vec[10]),
                     get(Ag_vec[11]),get(Ag_vec[12]),get(Ag_vec[13]),get(Ag_vec[14]),get(Ag_vec[15]),get(Ag_vec[16]),get(Ag_vec[17]))

Ag_name <- gsub("_percent","",Ag_vec)

heatmap_rownames <- list()
for (i in 1:17) {
  heatmap_rownames[[i]] <- rep(Ag_name[i],nrow(get(Ag_vec[i])))
}

column_name <- rowAnnotation(Ag_group = unlist(heatmap_rownames))

library(ComplexHeatmap)
p <- Heatmap(all_combine,
             name = "Percentage",
             column_title = "Clusters",
             row_title = "Clones",
             row_names_gp = gpar(fontsize = 12),
             column_names_gp = gpar(fontsize = 12),
             left_annotation = column_name,
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             col=colorRamp2(c(0, 100), c("white", "maroon"))
)

pdf(paste(savedir,"clones_heatmap_mem_clusters.pdf",sep = ""), width = 7, height = 9)
p
dev.off()

all_combine$Age_Antigen <- unlist(heatmap_rownames)
write.table(all_combine, paste(savedir,"clones_heatmap_mem_clusters.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

### For all Including naive #####
required_col <- select(CD8_all@meta.data, "TRA1_or_TRA2_TRB1_no_TRB2", "Age_Ag", "celltypes", "seurat_clusters")
required_col_noNA <- required_col[!is.na(required_col$TRA1_or_TRA2_TRB1_no_TRB2),]
Age_Antigens <- unique(required_col_noNA$Age_Ag)
Age_Antigens <- grep(",",Age_Antigens, invert = TRUE, value = TRUE)

savedir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scTCR/celltypes/Khuchroo_paper/naive/"
dir.create(savedir, showWarnings = FALSE)
write.table(required_col_noNA, paste(savedir,"clones_required_col_noNA.txt",sep = ""), sep="\t", col.names=T, row.names=T, quote=F)

### Running it in RStudio
rm(list = ls())
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scTCR/celltypes/Khuchroo_paper/naive/"
required_col_noNA <- read.table(paste(savedir,"clones_required_col_noNA.txt",sep = ""), header = TRUE, sep = "\t")
Age_Antigens <- grep(",",unique(required_col_noNA$Age_Ag)[order(unique(required_col_noNA$Age_Ag))],value = TRUE,invert = TRUE)

for (k in 1:length(Age_Antigens)) {
  try({
    required_col_noNA_Ag <- required_col_noNA[grep(Age_Antigens[k], required_col_noNA$Age_Ag),]
    required_col_noNA_Ag_table <- table(required_col_noNA_Ag$TRA1_or_TRA2_TRB1_no_TRB2, required_col_noNA_Ag$celltypes)
    rownames(required_col_noNA_Ag_table) <- NULL
    
    index <- vector()
    for (i in 1:nrow(required_col_noNA_Ag_table)) {
      ### removing those column that are not sharing with other clusters
      clone_count <- as.data.frame(table(required_col_noNA_Ag_table[i,]))
      total_col <- ncol(required_col_noNA_Ag_table)
      ### If total_col - 1 ==0 column has 0 we remove that column
      number = total_col-1
      if (clone_count[1,1] == 0) {
        if (clone_count[1,2] == number) {
          print(paste(i,"not Sharing",sep = " "))
          index[i] = NA
        }
        else{
          print(paste(i,"is Sharing",sep = " "))
          index[i] <- i
        }
      }
      else{
        index[i] <- i
      }
    }
    
    required_col_noNA_Ag_table_sharing <- required_col_noNA_Ag_table[index[!is.na(index)],]
    
    ### Making an empty dataframe
    Ag_sharing_percent <- as.data.frame(matrix(nrow=nrow(required_col_noNA_Ag_table_sharing), ncol=ncol(required_col_noNA_Ag_table_sharing)))
    colnames(Ag_sharing_percent) <- colnames(required_col_noNA_Ag_table_sharing)
    
    for (j in 1:nrow(required_col_noNA_Ag_table_sharing)) {
      Ag_sharing_percent[j,] <- (required_col_noNA_Ag_table_sharing[j,]/rowSums(required_col_noNA_Ag_table_sharing)[j])*100
    }
    
    assign(paste(Age_Antigens[k],"_percent",sep = ""),Ag_sharing_percent)
  })
}

Ag_vec <- ls(pattern = "_percent")
Ag_vec <- Ag_vec[-1]

### Please find the correct order of the celltypes
celltypes_order <- colnames(get(Ag_vec[2]))

### Making the equal columns with the same order
for (i in 1:length(Ag_vec)) {
  if (ncol(get(Ag_vec[i])) < length(unique(required_col_noNA$celltypes))) {
    df <- get(Ag_vec[i])
    missing_index <- grep(paste(colnames(df),collapse = "|"),celltypes_order,invert = TRUE,value = TRUE)
    df[,missing_index] <- 0
    
    ## ordering it
    df <- df[,match(celltypes_order,colnames(df))]
    assign(Ag_vec[i],df)
  }
}

## Combining all
all_combine <- rbind(get(Ag_vec[1]),get(Ag_vec[2]),get(Ag_vec[3]),get(Ag_vec[4]),get(Ag_vec[5]),get(Ag_vec[6]),get(Ag_vec[7]),get(Ag_vec[8]),get(Ag_vec[9]),get(Ag_vec[10]),
                     get(Ag_vec[11]),get(Ag_vec[12]),get(Ag_vec[13]),get(Ag_vec[14]),get(Ag_vec[15]),get(Ag_vec[16]),get(Ag_vec[17]))

Ag_name <- gsub("_percent","",Ag_vec)

heatmap_rownames <- list()
for (i in 1:17) {
  heatmap_rownames[[i]] <- rep(Ag_name[i],nrow(get(Ag_vec[i])))
}

column_name <- rowAnnotation(Ag_group = unlist(heatmap_rownames))

library(ComplexHeatmap)
p <- Heatmap(all_combine,
             name = "Percentage",
             column_title = "Celltypes",
             row_title = "Clones",
             row_names_gp = gpar(fontsize = 12),
             column_names_gp = gpar(fontsize = 12),
             left_annotation = column_name,
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             col=colorRamp2(c(0, 100), c("white", "maroon"))
)

pdf(paste(savedir,"clones_heatmap_naive_celltypes.pdf",sep = ""), width = 7, height = 9)
p
dev.off()

all_combine$Age_Antigen <- unlist(heatmap_rownames)
write.table(all_combine, paste(savedir,"clones_heatmap_naive_celltypes.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")

###  Clusters ####
rm(list = ls())
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/BEAM_T/downstream_analysis/after_geno_demux/Antigen/scTCR/celltypes/Khuchroo_paper/naive/"
required_col_noNA <- read.table(paste(savedir,"clones_required_col_noNA.txt",sep = ""), header = TRUE, sep = "\t")
Age_Antigens <- grep(",",unique(required_col_noNA$Age_Ag)[order(unique(required_col_noNA$Age_Ag))],invert = TRUE, value = TRUE)

for (k in 1:length(Age_Antigens)) {
  try({
    required_col_noNA_Ag <- required_col_noNA[grep(Age_Antigens[k], required_col_noNA$Age_Ag),]
    required_col_noNA_Ag_table <- table(required_col_noNA_Ag$TRA1_or_TRA2_TRB1_no_TRB2, required_col_noNA_Ag$seurat_clusters)
    rownames(required_col_noNA_Ag_table) <- NULL
    
    index <- vector()
    for (i in 1:nrow(required_col_noNA_Ag_table)) {
      ### removing those column that are not sharing with other clusters
      clone_count <- as.data.frame(table(required_col_noNA_Ag_table[i,]))
      total_col <- ncol(required_col_noNA_Ag_table)
      ### If total_col - 1 ==0 column has 0 we remove that column
      number = total_col-1
      if (clone_count[1,1] == 0) {
        if (clone_count[1,2] == number) {
          print(paste(i,"not Sharing",sep = " "))
          index[i] = NA
        }
        else{
          print(paste(i,"is Sharing",sep = " "))
          index[i] <- i
        }
      }
      else{
        index[i] <- i
      }
    }
    
    required_col_noNA_Ag_table_sharing <- required_col_noNA_Ag_table[index[!is.na(index)],]
    
    ### Making an empty dataframe
    Ag_sharing_percent <- as.data.frame(matrix(nrow=nrow(required_col_noNA_Ag_table_sharing), ncol=ncol(required_col_noNA_Ag_table_sharing)))
    colnames(Ag_sharing_percent) <- colnames(required_col_noNA_Ag_table_sharing)
    
    for (j in 1:nrow(required_col_noNA_Ag_table_sharing)) {
      Ag_sharing_percent[j,] <- (required_col_noNA_Ag_table_sharing[j,]/rowSums(required_col_noNA_Ag_table_sharing)[j])*100
    }
    
    assign(paste(Age_Antigens[k],"_percent",sep = ""),Ag_sharing_percent)
  })
}


Ag_vec <- ls(pattern = "_percent")
Ag_vec <- Ag_vec[-1]

### Please find the correct order of the celltypes
celltypes_order <- as.data.frame(table(required_col_noNA$seurat_clusters))[,1]

### Making the equal columns with the same order
for (i in 1:length(Ag_vec)) {
  if (ncol(get(Ag_vec[i])) < length(unique(required_col_noNA$seurat_clusters))) {
    df <- get(Ag_vec[i])
    missing_index <- grep(paste("^",colnames(df),"$",collapse = "|",sep = ""),celltypes_order,invert = TRUE,value = TRUE)
    df[,missing_index] <- 0
    
    ## ordering it
    df <- df[,match(celltypes_order,colnames(df))]
    assign(Ag_vec[i],df)
  }
}

## Combining all
all_combine <- rbind(get(Ag_vec[1]),get(Ag_vec[2]),get(Ag_vec[3]),get(Ag_vec[4]),get(Ag_vec[5]),get(Ag_vec[6]),get(Ag_vec[7]),get(Ag_vec[8]),get(Ag_vec[9]),get(Ag_vec[10]),
                     get(Ag_vec[11]),get(Ag_vec[12]),get(Ag_vec[13]),get(Ag_vec[14]),get(Ag_vec[15]),get(Ag_vec[16]),get(Ag_vec[17]))

Ag_name <- gsub("_percent","",Ag_vec)

heatmap_rownames <- list()
for (i in 1:17) {
  heatmap_rownames[[i]] <- rep(Ag_name[i],nrow(get(Ag_vec[i])))
}

column_name <- rowAnnotation(Ag_group = unlist(heatmap_rownames))
library(ComplexHeatmap)
p <- Heatmap(all_combine,
             name = "Percentage",
             column_title = "Clusters",
             row_title = "Clones",
             row_names_gp = gpar(fontsize = 12),
             column_names_gp = gpar(fontsize = 12),
             left_annotation = column_name,
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             col=colorRamp2(c(0, 100), c("white", "maroon"))
)

pdf(paste(savedir,"clones_heatmap_naive_clusters.pdf",sep = ""), width = 7, height = 9)
p
dev.off()

all_combine$Age_Antigen <- unlist(heatmap_rownames)
write.table(all_combine, paste(savedir,"clones_heatmap_naive_clusters.txt",sep = ""), quote = F, row.names = F, col.names = T, sep = "\t")







