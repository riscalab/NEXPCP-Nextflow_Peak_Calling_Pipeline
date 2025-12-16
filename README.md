# NEXPCP-Nextflow_Peak_Calling_Pipeline

## This pipeline might need you to have a few files in the bin directory if you want to get back certain plots and heatmaps. However, if you would just like differential analysis then that should not be a problem

# Parameters that you may need to run the pipeline
### More coming soon

# Required parameters
```

--control_bams 
--wt_bams
--PE true or --SE true
--treatment_name 'H1low' or the first field name of your treatment file
--control_igg (if you are using --sicer2 parameter)
--wt_igg (if you are using --sicer2 parameter)

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

