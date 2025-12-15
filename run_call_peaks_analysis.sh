#!/bin/env bash

#SBATCH --mem=20GB
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=nextflow_chip
#SBATCH --partition=hpc_l40s_b


# --mail-type=FAIL,END

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

# nextflow call_peaks_analysis_pipeline.nf -profile peak_calling_analysis -resume \
# --control_bams 'bam_files/H1low_*{bam,bam.bai}' \
# --wt_bams 'bam_files/Scrm_*{bam,bam.bai}' \
# --make_html_report true \
# --dups_log './dup_info/*_dups.log' \
# --plot_idr 0.4 \
# --return_idr 0.4 \
# --treatment_name 'H1low' \
# --masterPeak100kb

# nextflow call_peaks_analysis_pipeline.nf -profile peak_calling_analysis -resume \
# --plot_idr 0.4 \
# --return_idr 0.4

########################################## sicer2 path ##########################################

# but only running it on H3k9me3
nextflow call_peaks_analysis_pipeline.nf -profile peak_calling_analysis -resume \
--control_bams '/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/bam_files/sicer2_h3k9me3_run/Scrm_H3k9me3*{bam,bam.bai}' \
--wt_bams '/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/bam_files/sicer2_h3k9me3_run/H1low_H3k9me3*{bam,bam.bai}' \
--make_html_report true \
--dups_log './dup_info/*_dups.log' \
--PE true \
--sicer2 true \
--sicer2_fdr '0.01' \
--sicer2_gap_size '8000' \
--sicer2_window_size '500' \
--make_html_report true \
--dups_log './dup_info/*_dups.log' \
--treatment_name 'H1low' \
--atac_analysis true \
--nochange_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/unchanging_ATAC_dH1scr_regulated_peaks.bed' \
--up_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/up_ATAC_dH1scr_regulated_peaks.bed' \
--down_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/down_ATAC_dH1scr_regulated_peaks.bed' \
--control_ATAC_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/atac_bams/atac_shifted_bams/merged_shifted_bams/dH1_ATAC_allmerged*{bam,bam.bai}' \
--treatment_ATAC_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/atac_bams/atac_shifted_bams/merged_shifted_bams/scr_ATAC_allmerged*{bam,bam.bai}' \
--rpgc_bigwig true \
--rpgc_num_effectiveGenomeSize '2864785220' \
--control_igg '/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/bam_files/sicer2_h3k9me3_run/igg_bams/Scrm*Igg*{bam,bam.bai}' \
--wt_igg '/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/bam_files/sicer2_h3k9me3_run/igg_bams/H1low*Igg*{bam,bam.bai}'\


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