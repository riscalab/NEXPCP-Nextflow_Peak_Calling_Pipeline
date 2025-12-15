
include {

    bedtools_flank_process;
    r_wilcox_test_process


}from '../modules/meth_analysis_modules.nf'



workflow meth_enrichment_analysis_workflow {


    take:

    methyl_bed_ch
    up_genes
    down_genes
    ref_genome_size


    main:

    // now to use bedtools to find the overlap or map of methyl coordinates in up and down gene regions

    // first get the promoter regions of up and down genes using bedtools flank

    bedtools_flank_process(up_genes, down_genes, ref_genome_size, methyl_bed_ch)

    // now I should run a wilcoxon test on the last column of the promoter up and promoter down mean counts

    up_promoter_methylation_ch = bedtools_flank_process.out.up_methyl_promoter_signal
    down_promoter_methylation_ch = bedtools_flank_process.out.down_methyl_promoter_signal

    r_wilcox_test_process(up_promoter_methylation_ch, down_promoter_methylation_ch)



    //emit:


}