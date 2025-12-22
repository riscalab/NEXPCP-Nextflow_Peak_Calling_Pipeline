#!/bin/env bash

#SBATCH --mem=20GB
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=nextflow_chip
#SBATCH --partition=hpc_l40_a


# --account=risc_condo_bank

#hpc_l40s_b
#hpc_a10_a 
#hpc_l40_a
#hpc_l40s
#hpc_l40_b
#hpc_bigmem_a
#risc_a
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
# --control_bams '../merged_bams/Scr_*{bam,bam.bai}' \
# --wt_bams '../merged_bams/H1lo_*{bam,bam.bai}' \
# --PE true \
# --sicer2 true \
# --sicer2_fdr '0.3' \
# --sicer2_gap_size '2000' \
# --make_html_report true \
# --plot_idr 0.4 \
# --return_idr 1 \
# --treatment_name 'H1lo' \
# --masterPeak100kb true \
# --atac_analysis true \
# --nochange_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/unchanging_ATAC_dH1scr_regulated_peaks.bed' \
# --up_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/up_ATAC_dH1scr_regulated_peaks.bed' \
# --down_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/down_ATAC_dH1scr_regulated_peaks.bed' \
# --control_ATAC_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/atac_bams/atac_shifted_bams/merged_shifted_bams/dH1_ATAC_allmerged*{bam,bam.bai}' \
# --treatment_ATAC_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/atac_bams/atac_shifted_bams/merged_shifted_bams/scr_ATAC_allmerged*{bam,bam.bai}' \
# --rpgc_num_effectiveGenomeSize '2864785220' \
# --control_igg '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/merged_bams/igg_bams/Scr_*AO.arnold*IgG*{bam,bam.bai}' \
# --wt_igg '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/merged_bams/igg_bams/H1lo*AO.arnold*IgG*{bam,bam.bai}'



#### try the h3k4me1 ######

# nextflow call_peaks_analysis_pipeline.nf -profile peak_calling_analysis -resume \
# --control_bams '/lustre/fs4/risc_lab/scratch/hcanaj/CutnTag_aligned_bams_hg38_B3B3_042524/Scrm*H3K4me1*_*{bam,bam.bai}' \
# --wt_bams '/lustre/fs4/risc_lab/scratch/hcanaj/CutnTag_aligned_bams_hg38_B3B3_042524/dH1*H3K4me1*_*{bam,bam.bai}' \
# --PE true \
# --narrowPeak_data true \
# --macs2 true \
# --make_html_report true \
# --plot_idr 0.4 \
# --return_idr 1 \
# --treatment_name 'dH1' \
# --atac_analysis true \
# --nochange_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/unchanging_ATAC_dH1scr_regulated_peaks.bed' \
# --up_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/up_ATAC_dH1scr_regulated_peaks.bed' \
# --down_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/down_ATAC_dH1scr_regulated_peaks.bed' \
# --control_ATAC_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/atac_bams/atac_shifted_bams/merged_shifted_bams/dH1_ATAC_allmerged*{bam,bam.bai}' \
# --treatment_ATAC_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/atac_bams/atac_shifted_bams/merged_shifted_bams/scr_ATAC_allmerged*{bam,bam.bai}' \
# --rpgc_num_effectiveGenomeSize '2864785220' \
# --control_igg '/lustre/fs4/risc_lab/scratch/hcanaj/CutnTag_aligned_bams_hg38_B3B3_042524/Scrm_IGG_r1_S1_001*{bam,bam.bai}' \
# --wt_igg '/lustre/fs4/risc_lab/scratch/hcanaj/CutnTag_aligned_bams_hg38_B3B3_042524/dH1_IGG_r1_S9_001*{bam,bam.bai}' \
# --rpgc_bigwig true


## try H3K4me3 ########

# nextflow call_peaks_analysis_pipeline.nf -profile peak_calling_analysis -resume \
# --control_bams '/lustre/fs4/risc_lab/scratch/hcanaj/CutnTag_aligned_bams_hg38_B3B3_042524/Scrm*H3K4me3*_*{bam,bam.bai}' \
# --wt_bams '/lustre/fs4/risc_lab/scratch/hcanaj/CutnTag_aligned_bams_hg38_B3B3_042524/dH1*H3K4me3*_*{bam,bam.bai}' \
# --PE true \
# --narrowPeak_data true \
# --macs2 true \
# --make_html_report true \
# --plot_idr 0.4 \
# --return_idr 1 \
# --treatment_name 'dH1' \
# --atac_analysis true \
# --nochange_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/unchanging_ATAC_dH1scr_regulated_peaks.bed' \
# --up_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/up_ATAC_dH1scr_regulated_peaks.bed' \
# --down_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/down_ATAC_dH1scr_regulated_peaks.bed' \
# --control_ATAC_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/atac_bams/atac_shifted_bams/merged_shifted_bams/dH1_ATAC_allmerged*{bam,bam.bai}' \
# --treatment_ATAC_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/atac_bams/atac_shifted_bams/merged_shifted_bams/scr_ATAC_allmerged*{bam,bam.bai}' \
# --rpgc_num_effectiveGenomeSize '2864785220' \
# --control_igg '/lustre/fs4/risc_lab/scratch/hcanaj/CutnTag_aligned_bams_hg38_B3B3_042524/Scrm_IGG_r1_S1_001*{bam,bam.bai}' \
# --wt_igg '/lustre/fs4/risc_lab/scratch/hcanaj/CutnTag_aligned_bams_hg38_B3B3_042524/dH1_IGG_r1_S9_001*{bam,bam.bai}' \
# --rpgc_bigwig true


#### try H3K27ac #######

# now i want to see

# nextflow call_peaks_analysis_pipeline.nf -profile peak_calling_analysis -resume \
# --control_bams '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_peak_call_with_new_sicer2/h3k27ac_bams/Scrm*H3K27ac*_*{bam,bam.bai}' \
# --wt_bams '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_peak_call_with_new_sicer2/h3k27ac_bams/dH1*H3K27ac*_*{bam,bam.bai}' \
# --PE true \
# --narrowPeak_data true \
# --macs2 true \
# --make_html_report true \
# --plot_idr 0.4 \
# --return_idr 1 \
# --treatment_name 'dH1' \
# --atac_analysis true \
# --nochange_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/unchanging_ATAC_dH1scr_regulated_peaks.bed' \
# --up_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/up_ATAC_dH1scr_regulated_peaks.bed' \
# --down_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/down_ATAC_dH1scr_regulated_peaks.bed' \
# --control_ATAC_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/atac_bams/atac_shifted_bams/merged_shifted_bams/dH1_ATAC_allmerged*{bam,bam.bai}' \
# --treatment_ATAC_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/atac_bams/atac_shifted_bams/merged_shifted_bams/scr_ATAC_allmerged*{bam,bam.bai}' \
# --rpgc_num_effectiveGenomeSize '2864785220' \
# --control_igg '/lustre/fs4/risc_lab/scratch/hcanaj/CutnTag_aligned_bams_hg38_B3B3_042524/Scrm_IGG_r1_S1_001*{bam,bam.bai}' \
# --wt_igg '/lustre/fs4/risc_lab/scratch/hcanaj/CutnTag_aligned_bams_hg38_B3B3_042524/dH1_IGG_r1_S9_001*{bam,bam.bai}' \
# --rpgc_bigwig true

#### try H3K27ac but with sicer2 #######

# now i want to see

nextflow call_peaks_analysis_pipeline.nf -profile peak_calling_analysis -resume \
--control_bams '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_peak_call_with_new_sicer2/h3k27ac_bams/Scrm*H3K27ac*_*{bam,bam.bai}' \
--wt_bams '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_peak_call_with_new_sicer2/h3k27ac_bams/dH1*H3K27ac*_*{bam,bam.bai}' \
--PE true \
--sicer2 true \
--sicer2_fdr '0.3' \
--sicer2_gap_size '0' \
--sicer2_window_size '50' \
--make_html_report true \
--plot_idr 0.4 \
--return_idr 1 \
--treatment_name 'dH1' \
--atac_analysis true \
--nochange_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/unchanging_ATAC_dH1scr_regulated_peaks.bed' \
--up_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/up_ATAC_dH1scr_regulated_peaks.bed' \
--down_ATAC_peaks_bed '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/nextflow_R_script_outputs/ATAC/down_ATAC_dH1scr_regulated_peaks.bed' \
--control_ATAC_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/atac_bams/atac_shifted_bams/merged_shifted_bams/dH1_ATAC_allmerged*{bam,bam.bai}' \
--treatment_ATAC_bam '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/test_atac_data_pipeline/atac_bams/atac_shifted_bams/merged_shifted_bams/scr_ATAC_allmerged*{bam,bam.bai}' \
--rpgc_num_effectiveGenomeSize '2864785220' \
--rpgc_bigwig true \
--control_igg '/lustre/fs4/risc_lab/scratch/hcanaj/CutnTag_aligned_bams_hg38_B3B3_042524/Scrm_Igg_neg_control*{bam,bam.bai}' \
--wt_igg '/lustre/fs4/risc_lab/scratch/hcanaj/CutnTag_aligned_bams_hg38_B3B3_042524/dH1_Igg_neg_control*{bam,bam.bai}' \
--hpc_partition 'risc_a' \
--hpc_account '-A risc_condo_bank'

# trying with the negative control igg instead of these below
# --control_igg '/lustre/fs4/risc_lab/scratch/hcanaj/CutnTag_aligned_bams_hg38_B3B3_042524/Scrm_IGG_r1_S1_001*{bam,bam.bai}' \
# --wt_igg '/lustre/fs4/risc_lab/scratch/hcanaj/CutnTag_aligned_bams_hg38_B3B3_042524/dH1_IGG_r1_S9_001*{bam,bam.bai}' 
# --control_igg '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/merged_bams/igg_bams/new_merged_igg/Scr*IgG*{bam,bam.bai}' \
# --wt_igg '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/merged_bams/igg_bams/new_merged_igg/H1lo*IgG*{bam,bam.bai}'


# not using the rna seq genes I called anymore
# using irene's 
# --wtvs_lowup '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/H1low_vs_Scrm_up_rna_seq_gr_obj.bed'  \
# --wtvs_lowdown_nochange '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/H1low_vs_Scrm_down_rna_seq_gr_obj.bed' \
# --wtvs_lownochanging '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/H1low_vs_Scrm_unchanging_rna_seq_gr_obj.bed' \

# \
# -with-dag peak_calling_nf_pipeline_flowchart.pdf \
# -preview

# --rpgc_effectiveGenomeSize true \


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