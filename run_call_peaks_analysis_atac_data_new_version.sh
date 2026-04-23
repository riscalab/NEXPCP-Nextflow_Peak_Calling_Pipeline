#!/bin/env bash

#SBATCH --mem=20GB
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=nextflow_chip
#SBATCH --partition=hpc_l40s_b
#hpc_l40s_b
#hpc_a10_a 
#hpc_l40s
#hpc_l40_b


source /lustre/fs4/home/rjohnson/.bashrc_rj_test.sh
# source /ru-auth/local/home/rjohnson/.bashrc_rj_test.sh # or use this, should be the same thing 


conda activate nextflow_three

# plot_idr = 0.05 // the default is 0.05 to report in the png plots which peaks passed this threshold
# return_idr = 1  // the default is all peaks will be returned even if the report plots the ones that pass a certain number. 1 as default here will give back all peaks like idr has set


##### parameters documentation ###########

# This pipeline only takes experiments that have 3 replicates per condition due to my implementation of IDR
# I will allow it to take 2 replicates soon 

# This peak calling workflow takes bam files with the two parameters below
# please note that you should provide the path to the bam files with a glob pattern that lets nextflow know which are your control and which are your wild-type
# also you need to have the bam index files in that directory. and make the glob pattern take both the bam and the index file as shown in the example below
# --control_bams : this is actually the treatment bams you have. I wll update this to be called treatment. The workflow does not take controls like (igg)
# --wt_bams : this would be your wt bams. 

# the pipeline uses idr. IDR creates plots and in those plots you want the colors to represent the correct data that passed the threshold so make the next two parameters the same value
# --plot_idr : takes a value between 0 and 1 to show on a plot how much of your peaks passed the threshold. this value does not affect the actual data recieved in your final peak file
# --return_idr: this value actually returns only the peaks that passed this idr threshold.

# the if you ran your pipeline through the labs fastq2bam pipeline, you might have a few log files that contain duplicate information. 
# you can provide any log files that hopefully contain duplicate information or any other quality metrics you would like to have multiqc aggrigate for you. this includes any fastqc output files contatin information about your fastq files.
# --make_html_report (default = false) : make this true if you want to have the workflow make an html report of your duplicate info or any quality check files you might have that works with multiqc (ex: files from fastqc)
# --dups_log : give a path to the correct dups.log file for all of the bam files you created. Make a glob pattern for them or copy them to a directory and provide the path to that dir with all the dups.log files in there

# IMPORTANT for making sure deseq2 puts your condition of interest (treatment condition) first in the experiment design formula
# --treatment_name : for your conditions please enter the contition treatment name with the same capitalization as what is seen as your input file name. This will let deseq2 have the treatment condition as the target factor
# --masterPeak100kb false : add this parameter but make it true if you want the peaks to ber merged by this amount; default is false, but one needs to be selected
# --masterPeak10kb false : add this parameter but make it true if you want the peaks to ber merged by this amount; default is false, but one needs to be selected
# --masterPeak30kb false : add this parameter but make it true if you want the peaks to ber merged by this amount; default is false, but one needs to be selected
##########################################

nextflow call_peaks_analysis_pipeline.nf -profile peak_calling_analysis -resume \
--control_bams '../atac_bams/atac_shifted_bams/dH1_*{bam,bam.bai}' \
--wt_bams '../atac_bams/atac_shifted_bams/scr_*{bam,bam.bai}' \
--make_html_report true \
--dups_log './dup_info/*_dups.log' \
--plot_idr 0.1 \
--return_idr 0.1 \
--PE \
--treatment_name 'dH1' \
--macs3 \
--narrowPeak_data true \
--atac_analysis true \
--tss_enrich_atac true \
--rna_treatment_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/proseq_merged_bams/H1low_proseq_mergedbams_dedup*{bam.bai,bam}' \
--rna_control_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/proseq_merged_bams/Scr_proseq_mergedbams_dedup*{bam.bai,bam}' \
--cpm_bigwig true \
--enrichTSS_bed_region '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/NEXPCP-Nextflow_Peak_Calling_Pipeline/bin/hg38_TSS_export_bedrearrange_gene.bed'
#  \
# --hpc_account '-A risc_condo_bank' \
# --hpc_partition 'risc_a'



# --control_histone_signal '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/merged_bams/Scr_H3K27me3*{bam,bam.bai}' \
# --treatment_histone_signal '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/merged_bams/H1lo_H3K27me3*{bam,bam.bai}' \
# --control_ATAC_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/atac_bams/atac_shifted_bams/merged_shifted_bams/dH1_ATAC_allmerged*{bam,bam.bai}' \
# --treatment_ATAC_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/atac_bams/atac_shifted_bams/merged_shifted_bams/scr_ATAC_allmerged*{bam,bam.bai}' \
# --cpm_bigwig true

# need to try 0.1 idr and try 0.05 idr, I normally use 0.4
# keeping the macs3 pvalue at 0.05
############ Let's try using sicer2 with narrow peak parameters ###############
# I would need to find the igg for the atac seq bam files to do this so, not trying now










###############################################################################################

# nextflow call_peaks_analysis_pipeline.nf -profile peak_calling_analysis -resume \
# --plot_idr 0.4 \
# --return_idr 0.4


########################################## end seq data ##########################################
# nextflow call_peaks_analysis_pipeline.nf -profile peak_calling_analysis -resume \
# --control_bams 'bam_files/H1low_*{bam,bam.bai}' \
# --wt_bams 'bam_files/Scrm_*{bam,bam.bai}' \
# --make_html_report true \
# --dups_log './dup_info/*_dups.log' \
# --plot_idr 0.4 \
# --return_idr 0.4 \
# --treatment_name 'H1low'\
# --masterPeak100kb

####################################################################################