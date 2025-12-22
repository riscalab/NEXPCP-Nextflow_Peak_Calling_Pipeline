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
atac_analysis = false // if the user wants to do the atac seq analysis, they will need to add the path to bam files using the parameters (control_histone_signal, treatment_histone_signal)
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


# More coming soon

