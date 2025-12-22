# NEXPCP-Nextflow_Peak_Calling_Pipeline

## This pipeline might need you to have a few files in the bin directory if you want to get back certain plots and heatmaps. However, if you would just like differential analysis then that should not be a problem

# Parameters that you may need to run the pipeline
### More coming soon


# Template files
```
Use these template files to show you how I choose parameters below to run my data

run_call_peaks_analysis_atac_data_version.sh
run_call_peaks_analysis_hera_active_marks.sh
run_call_peaks_analysis_merged_hera_arnold_version.sh

```

# Required parameters
```

--control_bams 
--wt_bams
--PE true or --SE true
--treatment_name 'H1low' or the first field name of your treatment file
--control_igg (if you are using --sicer2 parameter)
--wt_igg (if you are using --sicer2 parameter)

```

# Parameters for controlling the partitions, account and executor
```
hpc_partition = 'hpc_a10_a' // this will be the default, and if you want to use risc_a change this to risc_a and the account to risc_condo_bank and not risc_hotel_bank

hpc_account = '-A risc_condo_bank' // for lab members you will pair this with hpc_partition = 'risc_a'
hpc_account = '-A risc_hotel_bank' This is the default. This is what will allow people to choose a partition other than risc_a. if they want to use risc a, they will have to put in the command line '-A risc_condo_bank' 
hpc_executor = 'slurm' // for lab members, this is default and will not need to be changed. for others, I do not know how this will behave when running the pipeline on a different executor and if the other options in this section are even relevant with different executors


```

### The parameters here represent the default values used in the pipeline, you, as the user, will set your own
```
plot_idr = 0.05 // the default is 0.05 to report in the png plots which peaks passed this threshold
return_idr = 1  // the default is all peaks will be returned even if the report plots the ones that pass a certain number. 1 as default here will give back all peaks like idr has set
make_html_report = false
treatment_name = null
masterPeak100kb = false
masterPeak10kb = false
masterPeak30kb = false
dups_log = []
narrowPeak_data = false
broadPeak_data = false
atac_analysis = false // if the user wants to do the atac seq analysis, see section below to know how to use other parameters that go with this one
PE = false // this is for the pair end reads experiment so the pipeline knows what it is 
SE = false // this is for the single end reads experiment so the pipeline knows how to handle it
macs2 = false
raw_bigwig = false
cpm_bigwig = false
rpgc_bigwig = false
//rpgc_effectiveGenomeSize = false // this is needed for the rpgc normalized bigwig
rpgc_num_effectiveGenomeSize = '2864785220' // this will be the default we use. most likely the size for the hg38
sicer2 = false // not sure why this parameter was not here yet
sicer2_fdr = '0.3'
sicer2_gap_size = '2000'
sicer2_window_size = '200' // this is the defualt value in sicer2, now the user can change it


```

### pro-seq and rna-sseq parameters 

**use this parameter if you want to change the default knowngene file from hg38 to a different organism or version**
```
--knownGene_bed_file 'path to one knowngene bed file' 
```

**use this if you have pro-seq data and a bed file with the up regulated genes from an experiment of the same type as your peak calling experiment used. The up, down, unchanging genes should come from the experiment design where treatment vs control(wt)**
```
--proseq_up_genes 'path to one up genes bed file' 
--proseq_down_genes 
--proseq_unchnaging_genes 
```


**Use these parameters if you have pro-seq bam files and want to get a bar plot showing the log2fc over roadmap epigenetic states Must be a pair glob pattern, with example**
```
--rna_treatment_bam  '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/proseq_merged_bams/H1low_proseq_mergedbams_dedup*{bam.bai,bam}'
--rna_control_bam
```

**Use these parameters below to pass in roadmap epigenetic state bed files which can be a glob pattern that points to all of the type of peaks in the directory**

```
--histone_broad  '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/roadmap_peaks_to_use/*.broadPeak'
--histone_narrow  '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/roadmap_peaks_to_use/*.narrowPeak'
```

## This workflow in the pipeline that uses these two parameters is not finished yet (workflow: get_proximal_distal_atac_peaks_workflow)
### If running narrowPeak_data parameter (for atac-seq data or narrow histone marks) data

**This section should not be confused with the atac_analysis path below. This is for when you are running atac-seq or narrow histone mark data through the pipeline as the main data you're analyzing. Use these parameters when you have histone bam files and want to have a bar graph showing log2fc of this signal over any of your atac peaks that are proximal or distal from gene tss. But this workflow was not fully complete so these two parameters are not needed right now**
```
--control_histone_signal
--treatment_histone_signal

```

### ATAC-seq analysis parameters


**I mentioned this a few sections up, but for running atac analysis just set this and it will be true or you can write true next to it in the command line. Next, specify the locations of your atac-seq bam files and atac seq peak bed files so you can get heat maps and profile plots where signal from the current assay (example cut&tag ) will have  cut&tag signal plotted over peaks or atac-seq signal plotted over cut&tag signal. The atac analysis will also calculate the log2FC of the atac-seq signal over different epigenetic states**
```
--atac_analysis true or false (default false) 
```


**Use these two parameters to get your ATAC-seq bigwig into the pipeline to get heatmaps that involve atac-seq data. You might need to run this pipeline first through the atac seq path to get the atac seq bigwigs and peak files. Then run the pipeline again in a different directory where you specify the path that points to the bigwig outputs when using this**

```
--control_ATAC_bigwig ' path to one control atac bigwig, so probably merge if you have replicates '
--treatment_ATAC_bigwig
```

**These are the atac seq bam files, which are used to create a bargraph that shows log2fc of atac signal over different epigenetic states**

```
--control_ATAC_bam '/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/scr-allmerge*.{bai,bam}'
--treatment_ATAC_bam
```

**If you ran the atac seq path of this pipeline or have your own atac seq up, down, and unchanging peaks, then you can put a path that points to them here so you can get heatmaps of the cut&tag (or any peak based assay) signal over these peak**

```
--down_ATAC_peaks_bed 'path to one down atac seq bed file'
--up_ATAC_peaks_bed
--nochange_ATAC_peaks_bed

```


# More coming soon

