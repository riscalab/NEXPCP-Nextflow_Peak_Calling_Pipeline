

process bedtools_flank_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'

    publishDir "./methylation_results", mode: 'copy', pattern: '*'

    label 'normal_big_resources'


    input:

    path(up_genes)

    path(down_genes)

    path(ref_genome_size)

    path(methyl_bed)



    output:

    path("${methyl_signal_in_up_promoters}"), emit: up_methyl_promoter_signal
    path("${methyl_signal_in_down_promoters}"), emit: down_methyl_promoter_signal


    script:

    methyl_bed_sorted = "${methyl_bed.baseName}_sorted.bed"

    up_genes_promoter_out = "up_genes_3kb_promoter.bed"
    down_genes_promoter_out = "down_genes_3kb_promoter.bed"

    up_genes_promoter_out_sorted = "up_genes_3kb_promoter_sorted.bed"
    down_genes_promoter_out_sorted = "down_genes_3kb_promoter_sorted.bed"

    methyl_signal_in_up_promoters = "up_promoters_3kb_methylSignal.bed"
    methyl_signal_in_down_promoters = "down_promoters_3kb_methylSignal.bed"


    """
    #!/usr/bin/env bash

    # flanking the up genes first
    bedtools flank -i ${up_genes} \
    -g ${ref_genome_size} \
    -b 3000 \
    > ${up_genes_promoter_out}

    # flanking the down genes now

    bedtools flank -i ${down_genes} \
    -g ${ref_genome_size} \
    -b 3000 \
    > ${down_genes_promoter_out}


    # need to sort both up and down promoter genes out files

    sort -k1,1 -k2,2n ${up_genes_promoter_out} > ${up_genes_promoter_out_sorted}

    sort -k1,1 -k2,2n ${down_genes_promoter_out} > ${down_genes_promoter_out_sorted}

    # now sorting the methyl bed file
    sort -k1,1 -k2,2n ${methyl_bed} > ${methyl_bed_sorted}

    # now I should be able to use bedtools map to get how many methylated sites are in upgenes

    bedtools map \
    -a ${up_genes_promoter_out_sorted} \
    -b ${methyl_bed_sorted} \
    -c 4 \
    -o mean \
    > ${methyl_signal_in_up_promoters}

    # now for the down promoters methyl signal

    bedtools map \
    -a ${down_genes_promoter_out_sorted} \
    -b ${methyl_bed_sorted} \
    -c 4 \
    -o mean \
    > ${methyl_signal_in_down_promoters}

    # the very last column is what will have the methylation mean in the signal in promoters file. so for this it would be column 8

    """
}



process r_wilcox_test_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/new_r_lan_3_rj'

    publishDir ""

    label 'normal_big_resources'



    input:

    path(up_methyl_promoter_data)

    path(down_methyl_promoter_data)



    output:



    script:


    """
    #!/usr/bin/env Rscript

    up_vals   <- read.table("${up_methyl_promoter_data}")[,8]
    down_vals <- read.table("${down_methyl_promoter_data}")[,8]
    wilcox_data = wilcox.test(as.numeric(up_vals), as.numeric(down_vals))

    # since the median is 0 have to switch to fisher's exact test instead





    """
}