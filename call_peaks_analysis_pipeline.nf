nextflow.enable.dsl=2

// this is the path where i got hera's bam files from: /lustre/fs4/risc_lab/scratch/hcanaj/K562_cutandtag_H1low_alignedbam_012925_hg38/all_aligned_bams_bw_hg38

// the bams and bai files should have 3 replicates if you want it to work with idr
params.control_bams = 'bam_files/H1low_*{bam,bam.bai}'
params.other_control_bams = 'bam_files/dH1_*{bam,bam.bai}'
//params.published_bam = "bam_files/ENCFF915XIL_*{bam,bam.bai}"
//params.published_bam = "bam_files/ENCFF905CZD_*{bam,bam.bai}"

norm_control_bams_index_tuple_ch = Channel.fromFilePairs(params.control_bams)
other_control_bams_index_tuple_ch = Channel.fromFilePairs(params.other_control_bams)
//published_bam_index_tuple_ch = Channel.fromFilePairs(params.published_bam)

control_bams_index_tuple_ch = norm_control_bams_index_tuple_ch.concat(other_control_bams_index_tuple_ch )

//control_bams_index_tuple_ch.view()


params.wt_bams = 'bam_files/Scrm_*{bam,bam.bai}'
wt_bams_index_tuple_ch = Channel.fromFilePairs(params.wt_bams)

params.control_igg = '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/merged_bams/igg_bams/Scr*{bam,bam.bai}'
control_igg_bam_index_tuple_ch = Channel.fromFilePairs(params.control_igg)

params.wt_igg = '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/merged_bams/igg_bams/H1lo*{bam,bam.bai}'
wt_igg_bam_index_tuple_ch = Channel.fromFilePairs(params.wt_igg)


params.ref_genome = file('/lustre/fs4/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/WholeGenomeFasta/genome.fa')
ref_genome_ch = Channel.value(params.ref_genome)

params.ref_genome_size = file('/lustre/fs4/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/WholeGenomeFasta/genome.fa.fai')
ref_genome_size_ch = Channel.value(params.ref_genome_size)

params.gtf_file = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/RNA-seq/rep2/gencode.v38_ERCC.gtf')
gtf_ch = Channel.value(params.gtf_file)

// I need a list of gene regions. got this from irene but also included it in my motif analysis nextflow script
params.wtvs_lowup = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/RNA-seq/rep2/hg38-ERCC-UMI-alignment/DESeq2_results/WTvslowup-genebody.bed')
wtvslowup_genebody_ch = Channel.value(params.wtvs_lowup)

// including the nochange genes just to have something to compare to
// changing the file from WTvslowdown-basemeanmatchnochange-genebody.bed, to this WTvslowdown-genebody.bed
params.wtvs_lowdown_nochange = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/RNA-seq/rep2/hg38-ERCC-UMI-alignment/DESeq2_results/WTvslowdown-genebody.bed')
wtvslowdown_nochange_ch = Channel.value(params.wtvs_lowdown_nochange)

// now to get irene's nochange gene body file
params.wtvs_lownochanging = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/RNA-seq/rep2/hg38-ERCC-UMI-alignment/DESeq2_results/WTvslownochange-genebody.bed')
wtvslownochange = Channel.value(params.wtvs_lownochanging)

// I want a multiqc plot of the duplicate levels for the k27 data

if (params.dups_log) {
    params.dups_log = file('./dup_info/*_dups.log')
    dups_log_ch = Channel.fromPath(params.dups_log)
}
else {
    dups_log_ch = Channel.empty()
}


params.meth_bed_file = file('/lustre/fs4/home/rjohnson/pipelines/hera_pipeline/results_SE/bl_filt_bed/bed_graphs_deeptools/control_CpG_r1r2r3_fp_filt_filt_coor_sorted_BL_filt_sort2_normalized_cpm.bed')
meth_bed_ch = Channel.value(params.meth_bed_file)


params.up_peaks_file = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/test_published_data/up_regulated_peaks.bed')
up_peaks_ch = Channel.value(params.up_peaks_file)

params.down_peaks_file = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/test_published_data/down_regulated_peaks.bed')
down_peaks_ch = Channel.value(params.down_peaks_file)

params.unchanging_peaks_file = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/test_published_data/unchanging_regulated_peaks.bed')
unchanging_peaks_ch = Channel.value(params.unchanging_peaks_file)

params.master_peaks_file = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/test_published_data/10kb_merged_masterPeak_geo_control_h1low.bed')
master_peaks_ch = Channel.value(params.master_peaks_file)

params.knownGene_bed_file = file('/lustre/fs4/home/rjohnson/downloads/genomes/hg38/genes/knownGene.hg38.bed')
knownGene_ch = Channel.value(params.knownGene_bed_file)

// proseq_up_genesID_with_coord.bed
// params.proseq_up_genes = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/proseq_data/proseq_up_genes_with_coord.tsv')
params.proseq_up_genes = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/proseq_data/proseq_up_genesID_with_coord.bed')
proseq_up_gene_ch = Channel.value(params.proseq_up_genes)

// proseq_down_genesID_with_coord.bed
// params.proseq_down_genes = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/proseq_data/proseq_down_genes_with_coord.tsv')
params.proseq_down_genes = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/proseq_data/proseq_down_genesID_with_coord.bed')
proseq_down_gene_ch = Channel.value(params.proseq_down_genes)

// proseq_other_genesID_with_coord.bed
// params.proseq_unchanging_genes = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/proseq_data/proseq_unchanging_genes_with_coord.tsv')
params.proseq_unchanging_genes = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/proseq_data/proseq_other_genesID_with_coord.bed')
proseq_unchanging_gene_ch = Channel.value(params.proseq_unchanging_genes)


// now adding the ATAC-seq bigwig signal files to nextflow

params.control_ATAC_bigwig = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/tracks/scr-allmerge.bw')
control_atac_bigwig_ch = Channel.value(params.control_ATAC_bigwig)

params.treatment_ATAC_bigwig = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/tracks/dH1-allmerge.bw')
treatment_atac_bigwig_ch = Channel.value(params.treatment_ATAC_bigwig)


// lets try the atac-bam files

params.control_ATAC_bam = '/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/scr-allmerge*.{bai,bam}'
control_atac_bam_ch =  Channel.fromFilePairs(params.control_ATAC_bam)

params.treatment_ATAC_bam = '/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/dH1-allmerge*.{bai,bam}'
treatment_atac_bam_ch = Channel.fromFilePairs(params.treatment_ATAC_bam)

// now i want to add the CpG island bigwig file so I can plot the signal over the up and down peaks 
params.bisulfate_bigwig = file('/lustre/fs4/home/rjohnson/pipelines/hera_pipeline/results_SE/bl_filt_bed/bed_graphs_deeptools/control_CpG_r1r2r3_fp_filt_filt_coor_sorted_BL_filt_sort2_normalized_cpm.bigwig')
// this is the version that uses field 4 in the bedgraph generated by methyldackel
//params.bisulfate_bigwig = file('/lustre/fs4/home/rjohnson/pipelines/hera_pipeline/results_SE/methylDackel_CpG_data/control_CpG_r1r2r3_fp_filt_filt_coor_sorted_BL_filt_sort2_CpG.bedGraph_bedGraph_to.bigwig')
bisulfate_bigwig_ch = Channel.value(params.bisulfate_bigwig)

// there is geo data for actual cpg islands in bed format that I have
params.cpg_islands_unmasked_bed = file('bin/cpgIslands_unmasked_real_islands.bed')
cpg_island_unmasked_ch = Channel.value(params.cpg_islands_unmasked_bed)

// now adding the ATAC-seq peaks up and down, to plot the atac-seq signal and the k27me3 signal of scrm and h1low
params.down_ATAC_peaks_bed = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/DEseq/scrvslow-down.bed')
down_atac_peaks_ch = Channel.value(params.down_ATAC_peaks_bed)

params.up_ATAC_peaks_bed = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/DEseq/scrvslow-up.bed')
up_atac_peaks_ch = Channel.value(params.up_ATAC_peaks_bed)

// this would be the no change ATAC peaks 
params.nochange_ATAC_peaks_bed = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/DEseq/rep-peaks-scrvslow-nochange.bed')
nochange_atac_peaks_ch = Channel.value(params.nochange_ATAC_peaks_bed)

// get histone bam files from the user if they are running ATAC-seq analysis
params.control_histone_signal = '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/merged_bams/Scr_H3K27me3*{bam,bai}'
control_histone_bams = Channel.fromFilePairs(params.control_histone_signal)

params.treatment_histone_signal = '/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/merged_bams/H1lo_H3K27me3*{bam,bai}'
treatment_histone_bams = Channel.fromFilePairs(params.treatment_histone_signal)

// get the proseq bam files from the user if they are running the  ATAC-seq analysis
// first I will merge the pro-seq bam files and then use the merged files as input here
params.rna_treatment_bam = file('/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/proseq_merged_bams/H1low_proseq_mergedbams_dedup*{bam.bai,bam}')
treatment_proseq_bam = Channel.value(params.rna_treatment_bam)

params.rna_control_bam = file('/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/proseq_merged_bams/Scr_proseq_mergedbams_dedup*{bam.bai,bam}')
control_proseq_bam = Channel.value(params.rna_control_bam)

// control_proseq_bam.view{it -> "this is the input fromFilePairs for the control proseq: $it"}


/////////////////// Now for the roadmap bedfiles //////////////////////
// location:  /lustre/fs4/home/ascortea/store/ascortea/beds/k562
// params.k27me3_broad= file('/lustre/fs4/home/ascortea/store/ascortea/beds/k562/E123-H3K27me3.broadPeak')
// params.h2az_broad = file('/lustre/fs4/home/ascortea/store/ascortea/beds/k562/E123-H2A.Z.broadPeak')
// params.k27ac_broad = file('/lustre/fs4/home/ascortea/store/ascortea/beds/k562/E123-H3K27ac.broadPeak')
// params.k36me3_broad = file('/lustre/fs4/home/ascortea/store/ascortea/beds/k562/E123-H3K36me3.broadPeak')
// params.k4me1_broad = file('/lustre/fs4/home/ascortea/store/ascortea/beds/k562/E123-H3K4me1.broadPeak')
// params.k4me2_broad = file('/lustre/fs4/home/ascortea/store/ascortea/beds/k562/E123-H3K4me2.broadPeak')
// params.k4me3_broad = file('/lustre/fs4/home/ascortea/store/ascortea/beds/k562/E123-H3K4me3.broadPeak')
// params.k79me2_broad = file('/lustre/fs4/home/ascortea/store/ascortea/beds/k562/E123-H3K79me2.broadPeak')
// params.k9ac_broad = file('/lustre/fs4/home/ascortea/store/ascortea/beds/k562/E123-H3K9ac.broadPeak')
// params.k9me1_broad = file('/lustre/fs4/home/ascortea/store/ascortea/beds/k562/E123-H3K9me1.broadPeak')
// params.k9me2_broad = file('/lustre/fs4/home/ascortea/store/ascortea/beds/k562/E123-H3K9me3.broadPeak')
// params.k20me1_broad = file('/lustre/fs4/home/ascortea/store/ascortea/beds/k562/E123-H4K20me1.broadPeak')

// many_broad_histone_peaks = Channel.value{params.k27me3_broad; params.h2az_broad; params.k27ac_broad; params.k36me3_broad; params.k4me1_broad; params.k4me2_broad; params.k4me3_broad; params.k79me2_broad; params.k9ac_broad; params.k9me1_broad; params.k9me2_broad; params.k20me1_broad  }

// or just get all broadpeak files using a glob pattern

params.histone_broad = file('/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/roadmap_peaks_to_use/*.broadPeak')
params.histone_narrow = file('/lustre/fs4/home/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/roadmap_peaks_to_use/*.narrowPeak')

roadmap_broad_histones = Channel.value(params.histone_broad)
roadmap_narrow_histones = Channel.value(params.histone_narrow)

//////////////////////////////////////////////////////////////////////
// include {
//     make_alignment_bw_process_control

// }from './modules/peak_analysis_modules.nf'



include {
    mk_bw_call_peaks_workflow;
    plot_histone_data_workflow;
    plot_histone_calledpeak_workflow;
    plot_signal_up_down_peaks_workflow;
    plot_diff_peaks_over_diff_genes_workflow;
    plot_atac_signal_over_diff_peaks_workflow;
    find_then_plot_cpgIslands_in_peaks_workflow;
    get_roadmap_histone_enrichment_workflow;
    find_then_plot_atacseqPeaks_in_experiment_peaks_workflow;
    get_proximal_distal_atac_peaks_workflow;
    mk_bw_call_peaks_workflow_sicer2
    //get_diff_peaks_workflow

}from './workflows/call_peaks_workflow.nf'



include {

}from './workflows/find_diff_peaks_workflow.nf'

include {
    meth_enrichment_analysis_workflow



}from './workflows/methylation_enrichment_analysis_workflow.nf'




workflow {

    // now i will put the control bams and wt bams into the peak calling workflow
    // I will also put the reference genome in as the third entry

    // mk_bw_call_peaks_workflow(control_bams_index_tuple_ch, wt_bams_index_tuple_ch, ref_genome_ch, ref_genome_size_ch, dups_log_ch )

    // // get a channel with the final concat idr peaks
    // //final_idr_concat_peaks_ch = mk_bw_call_peaks_workflow.out.concat_idr_peaks
    
    // // take the emitted channels from the call peaks workflow
    // control_bw_meta_ch = mk_bw_call_peaks_workflow.out.control_meta_bw_ch

    // wt_bw_meta_ch = mk_bw_call_peaks_workflow.out.wt_meta_bw_ch

    // // getting the cpm normalized bigwigs
    // control_meta_cpm_bw = mk_bw_call_peaks_workflow.out.control_meta_cpm_bw_ch
    // wt_meta_cpm_bw = mk_bw_call_peaks_workflow.out.wt_meta_cpm_bw_ch

    // // the mk_bw_call_peaks_workflow workflow also emits the broad peaks that were called so I will grab those
    // all_broadpeaks_ch = mk_bw_call_peaks_workflow.out.macspeaks_ch

    // master_peak_list_true = mk_bw_call_peaks_workflow.out.master_peaks_list_ch
    // up_peaks_list_true = mk_bw_call_peaks_workflow.out.up_peaks_list_ch
    // down_peaks_list_true = mk_bw_call_peaks_workflow.out.down_peaks_list_ch
    // unchanging_peaks_list_true = mk_bw_call_peaks_workflow.out.unchanging_peaks_list_ch

    // idr_merged_peaks = mk_bw_call_peaks_workflow.out.group_concat_idr_peaks_ch

    // adding a workflow here to get the differential peaks

    //get_diff_peaks_workflow(all_broadpeaks_ch)


    // here I want to make a workflow that plots the chromatin features at a list of annotated sites (genes)
    // Also looking at the peaks around tss

    ////////// for testing do not run the plotting yet it takes too long //////////////////////////

    // NOT DOING THIS WORKFLOW ANYMORE. IT IS POINTLES TO SHOW THIS SIGNAL WITHOUT THE SCRM AND H1LOW ON THE SAME PLOT
    //plot_histone_data_workflow(control_bw_meta_ch, wt_bw_meta_ch, wtvslowup_genebody_ch, wtvslowdown_nochange_ch)

    // I will just make another workflow for plotting each histone-rep bigwig with its corresponding histone-rep broadPeak

    // NOT DOING THIS EITHER BECASUE THE SINGAL OVER THE LOOSE PARAMETER BROADPEAKS FROM MACS IS NOT USEFUL INFO BEFORE THE IDR AND MERGING FILTERING 
    //plot_histone_calledpeak_workflow(control_bw_meta_ch, wt_bw_meta_ch, all_broadpeaks_ch)



    // now putting the methylation bed file in the workflow with the up genes and down genes

    //meth_enrichment_analysis_workflow(meth_bed_ch, wtvslowup_genebody_ch, wtvslowdown_nochange_ch, ref_genome_size_ch)



    // WILL HAVE TO FIGURE OUT A NAMING FOR THE PEAK FILES SO I CAN AUTOMATE THIS (IF HISTONE IN BIGWIG IS SAME AS HISTONE_LABEL IN PEAKS THEN RUN. SOMETHING LIKE THAT)
    // now lets plot histone marks signal (cpm normalized bigwig) over the up and down regulated peaks

    //plot_signal_up_down_peaks_workflow(control_meta_cpm_bw, wt_meta_cpm_bw, up_peaks_ch, down_peaks_ch, bisulfate_bigwig_ch, master_peaks_ch, cpg_island_unmasked_ch)

    // emitting the combined bigwig meta channel grouped by both histone and replicate
    // it has this format tuple(condition, histone, replicate, file_name, file)
    //combined_bigwig_meta_2grouped_ch = plot_signal_up_down_peaks_workflow.out.experiment_group_meta_cpm_ch

    

    if (params.narrowPeak_data) {

        mk_bw_call_peaks_workflow(control_bams_index_tuple_ch, wt_bams_index_tuple_ch, ref_genome_ch, ref_genome_size_ch, dups_log_ch )

        // get a channel with the final concat idr peaks
        //final_idr_concat_peaks_ch = mk_bw_call_peaks_workflow.out.concat_idr_peaks
        
        // take the emitted channels from the call peaks workflow
        // control_bw_meta_ch = mk_bw_call_peaks_workflow.out.control_meta_bw_ch

        // wt_bw_meta_ch = mk_bw_call_peaks_workflow.out.wt_meta_bw_ch

        control_bw_meta_ch = mk_bw_call_peaks_workflow.out.control_bigwig

        wt_bw_meta_ch = mk_bw_call_peaks_workflow.out.wt_bigwig

        // getting the cpm normalized bigwigs
        // control_meta_cpm_bw = mk_bw_call_peaks_workflow.out.control_meta_cpm_bw_ch
        // wt_meta_cpm_bw = mk_bw_call_peaks_workflow.out.wt_meta_cpm_bw_ch

        // the mk_bw_call_peaks_workflow workflow also emits the broad peaks that were called so I will grab those
        all_broadpeaks_ch = mk_bw_call_peaks_workflow.out.macspeaks_ch

        master_peak_list_true = mk_bw_call_peaks_workflow.out.master_peaks_list_ch
        up_peaks_list_true = mk_bw_call_peaks_workflow.out.up_peaks_list_ch
        down_peaks_list_true = mk_bw_call_peaks_workflow.out.down_peaks_list_ch
        unchanging_peaks_list_true = mk_bw_call_peaks_workflow.out.unchanging_peaks_list_ch

        idr_merged_peaks = mk_bw_call_peaks_workflow.out.group_concat_idr_peaks_ch

        /// WILL THIS WORK FOR THE ATAC-SEQ DATA?  ///////////
        /////////////////////////////////////// For when creating master peaks generated in this pipeline /////////////
        // this is replicating the process that uses peaks so i can use the peaks that were made in the pipeline
        plot_signal_up_down_peaks_workflow(control_bw_meta_ch, wt_bw_meta_ch, up_peaks_list_true, down_peaks_list_true, bisulfate_bigwig_ch, master_peak_list_true, unchanging_peaks_list_true, cpg_island_unmasked_ch)

        // the grouped channel that uses only data from the pipeline not geo control
        combined_bigwig_peak_2grouped_ch = plot_signal_up_down_peaks_workflow.out.exper_rep_bigwig_peak_group

        // if i want to keep using this channel, it also has to be below the workflow (only for when using geo control data)
        combined_bigwig_meta_2grouped_ch = plot_signal_up_down_peaks_workflow.out.experiment_group_meta_cpm_ch

        bigwig_meta_ch_to_join = plot_signal_up_down_peaks_workflow.out.experiment_group_meta_to_join_ch
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // now lets see how to get differential peaks to be plot over the up and down genes TSS plus 5kb
        // I will also plot the signal over genes tss plus 20kb in this workflow, because I already made the 20kb changed from 5kb
        // will put the combined bigwig channel here emitted from the other workflow
        //plot_diff_peaks_over_diff_genes_workflow(up_peaks_ch, down_peaks_ch, unchanging_peaks_ch, master_peaks_ch, wtvslowup_genebody_ch, wtvslowdown_nochange_ch, wtvslownochange, gtf_ch, ref_genome_size_ch, knownGene_ch, proseq_up_gene_ch, proseq_down_gene_ch, proseq_unchanging_gene_ch, combined_bigwig_meta_2grouped_ch)

        // now recreating the above that will only use the true peak files generated from the pipeline itself
        plot_diff_peaks_over_diff_genes_workflow(up_peaks_list_true, down_peaks_list_true, unchanging_peaks_list_true, master_peak_list_true, wtvslowup_genebody_ch, wtvslowdown_nochange_ch, wtvslownochange, gtf_ch, ref_genome_size_ch, knownGene_ch, proseq_up_gene_ch, proseq_down_gene_ch, proseq_unchanging_gene_ch, combined_bigwig_peak_2grouped_ch)
        
        // now to add a workflow to plot the atac-seq signal over the cut&run peaks



        //////////////////////////////////////////////////////////////////////////////////////////////////
        // need to update the below workflow to take the actual peaks generated in the pipeline
        //plot_atac_signal_over_diff_peaks_workflow(control_atac_bigwig_ch, treatment_atac_bigwig_ch, up_peaks_ch, down_peaks_ch, unchanging_peaks_ch, down_atac_peaks_ch, up_atac_peaks_ch, combined_bigwig_meta_2grouped_ch, cpg_island_unmasked_ch )
        
        // this below will be the copy of the above workflow but using actual peaks generated in pipeline
        // FOR THE ATAC SEQ PART OF THE PIPELINE THIS WILL PLOT THE ATAC BIGWIG SIGNAL OVER THE ATAC PEAKS SO SIMILAR TO ANOTHER PROCESS THAT DOES THAT IN THE OTHER PATH BUT WITHOUT USER DEFINED ATAC PEAK FILE
        plot_atac_signal_over_diff_peaks_workflow(control_atac_bigwig_ch, treatment_atac_bigwig_ch, up_peaks_list_true, down_peaks_list_true, unchanging_peaks_list_true, down_atac_peaks_ch, up_atac_peaks_ch, nochange_atac_peaks_ch, combined_bigwig_peak_2grouped_ch, cpg_island_unmasked_ch )
        ///////////////////////////////////////////////////////////////////////////////////////////////////


        // this will be to find the atac-seq peaks that are near the proseq genes

        // i want to get the up peaks that are near proseq gene tss or distal from them
        get_proximal_distal_atac_peaks_workflow(up_peaks_list_true, down_peaks_list_true, unchanging_peaks_list_true, master_peak_list_true, proseq_up_gene_ch, proseq_down_gene_ch, proseq_unchanging_gene_ch, ref_genome_size_ch, control_histone_bams, treatment_histone_bams)


        

    }
    else { 

        
        if (params.macs2) {
            mk_bw_call_peaks_workflow(control_bams_index_tuple_ch, wt_bams_index_tuple_ch, ref_genome_ch, ref_genome_size_ch, dups_log_ch )

            // get a channel with the final concat idr peaks
            //final_idr_concat_peaks_ch = mk_bw_call_peaks_workflow.out.concat_idr_peaks
            
            // take the emitted channels from the call peaks workflow
            // control_bw_meta_ch = mk_bw_call_peaks_workflow.out.control_meta_bw_ch

            // wt_bw_meta_ch = mk_bw_call_peaks_workflow.out.wt_meta_bw_ch

            control_bw_meta_ch = mk_bw_call_peaks_workflow.out.control_bigwig

            wt_bw_meta_ch = mk_bw_call_peaks_workflow.out.wt_bigwig

            // getting the cpm normalized bigwigs
            // control_meta_cpm_bw = mk_bw_call_peaks_workflow.out.control_meta_cpm_bw_ch
            // wt_meta_cpm_bw = mk_bw_call_peaks_workflow.out.wt_meta_cpm_bw_ch

            // the mk_bw_call_peaks_workflow workflow also emits the broad peaks that were called so I will grab those
            all_broadpeaks_ch = mk_bw_call_peaks_workflow.out.macspeaks_ch

            master_peak_list_true = mk_bw_call_peaks_workflow.out.master_peaks_list_ch
            up_peaks_list_true = mk_bw_call_peaks_workflow.out.up_peaks_list_ch
            down_peaks_list_true = mk_bw_call_peaks_workflow.out.down_peaks_list_ch
            unchanging_peaks_list_true = mk_bw_call_peaks_workflow.out.unchanging_peaks_list_ch

            idr_merged_peaks = mk_bw_call_peaks_workflow.out.group_concat_idr_peaks_ch

            // when using narrowPeak_data, lets not do these analyses

            /////////////////////////////////////// For when creating master peaks generated in this pipeline /////////////
            // this is replicating the process that uses peaks so i can use the peaks that were made in the pipeline
            plot_signal_up_down_peaks_workflow(control_bw_meta_ch, wt_bw_meta_ch, up_peaks_list_true, down_peaks_list_true, bisulfate_bigwig_ch, master_peak_list_true, unchanging_peaks_list_true, cpg_island_unmasked_ch)

            // the grouped channel that uses only data from the pipeline not geo control
            combined_bigwig_peak_2grouped_ch = plot_signal_up_down_peaks_workflow.out.exper_rep_bigwig_peak_group

            // if i want to keep using this channel, it also has to be below the workflow (only for when using geo control data)
            combined_bigwig_meta_2grouped_ch = plot_signal_up_down_peaks_workflow.out.experiment_group_meta_cpm_ch

            bigwig_meta_ch_to_join = plot_signal_up_down_peaks_workflow.out.experiment_group_meta_to_join_ch
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // now lets see how to get differential peaks to be plot over the up and down genes TSS plus 5kb
            // I will also plot the signal over genes tss plus 20kb in this workflow, because I already made the 20kb changed from 5kb
            // will put the combined bigwig channel here emitted from the other workflow
            //plot_diff_peaks_over_diff_genes_workflow(up_peaks_ch, down_peaks_ch, unchanging_peaks_ch, master_peaks_ch, wtvslowup_genebody_ch, wtvslowdown_nochange_ch, wtvslownochange, gtf_ch, ref_genome_size_ch, knownGene_ch, proseq_up_gene_ch, proseq_down_gene_ch, proseq_unchanging_gene_ch, combined_bigwig_meta_2grouped_ch)

            // now recreating the above that will only use the true peak files generated from the pipeline itself
            plot_diff_peaks_over_diff_genes_workflow(up_peaks_list_true, down_peaks_list_true, unchanging_peaks_list_true, master_peak_list_true, wtvslowup_genebody_ch, wtvslowdown_nochange_ch, wtvslownochange, gtf_ch, ref_genome_size_ch, knownGene_ch, proseq_up_gene_ch, proseq_down_gene_ch, proseq_unchanging_gene_ch, combined_bigwig_peak_2grouped_ch)
            
            // now to add a workflow to plot the atac-seq signal over the cut&run peaks



            //////////////////////////////////////////////////////////////////////////////////////////////////
            // need to update the below workflow to take the actual peaks generated in the pipeline
            //plot_atac_signal_over_diff_peaks_workflow(control_atac_bigwig_ch, treatment_atac_bigwig_ch, up_peaks_ch, down_peaks_ch, unchanging_peaks_ch, down_atac_peaks_ch, up_atac_peaks_ch, combined_bigwig_meta_2grouped_ch, cpg_island_unmasked_ch )
            
            // this below will be the copy of the above workflow but using actual peaks generated in pipeline
            plot_atac_signal_over_diff_peaks_workflow(control_atac_bigwig_ch, treatment_atac_bigwig_ch, up_peaks_list_true, down_peaks_list_true, unchanging_peaks_list_true, down_atac_peaks_ch, up_atac_peaks_ch, nochange_atac_peaks_ch, combined_bigwig_peak_2grouped_ch, cpg_island_unmasked_ch )
            ///////////////////////////////////////////////////////////////////////////////////////////////////


            ///////////////////////////////////////////////////////////////////////////////////////////////////
            // now we want to first find the CpG islands that are overlapping up, down, unchanging experiment(histone) peaks
            // then we plot the h1low, scrm and bisulfate signal over the new "CpG_in_up_peaks, CpG_in_down_peaks, and CpG_in_unchanging_peaks"

            //find_then_plot_cpgIslands_in_peaks_workflow(up_peaks_ch, down_peaks_ch, unchanging_peaks_ch, cpg_island_unmasked_ch, combined_bigwig_meta_2grouped_ch )

            // this is the version of the above workflow that will only use the peaks generated within the pipeline

            
            find_then_plot_cpgIslands_in_peaks_workflow(up_peaks_list_true, down_peaks_list_true, unchanging_peaks_list_true, master_peak_list_true, cpg_island_unmasked_ch, combined_bigwig_meta_2grouped_ch, combined_bigwig_peak_2grouped_ch, bigwig_meta_ch_to_join )
            ///////////////////////////////////////////////////////////////////////////////////////////////////

            // now I want to use roadmap histone data to get and enrichment barplot of which histones have more ATAC-accessability

            //get_roadmap_histone_enrichment_workflow(roadmap_broad_histones, control_atac_bigwig_ch, treatment_atac_bigwig_ch)
            
            get_roadmap_histone_enrichment_workflow(roadmap_broad_histones, roadmap_narrow_histones, idr_merged_peaks, control_atac_bam_ch, treatment_atac_bam_ch, control_proseq_bam, treatment_proseq_bam, up_peaks_list_true, down_peaks_list_true)


            ///////////////////////////////////////////////////////////////////////////////////////////////////
            // recreating the plot cpgIslands workflow but for finding which ATAC-seq peaks are in histone mark peaks

            find_then_plot_atacseqPeaks_in_experiment_peaks_workflow(up_peaks_list_true, down_peaks_list_true, unchanging_peaks_list_true, master_peak_list_true, nochange_atac_peaks_ch, up_atac_peaks_ch, roadmap_broad_histones, roadmap_narrow_histones, idr_merged_peaks, control_atac_bam_ch, treatment_atac_bam_ch)
            
        }
        else if (params.sicer2) {

            mk_bw_call_peaks_workflow_sicer2(control_bams_index_tuple_ch, wt_bams_index_tuple_ch, ref_genome_ch, ref_genome_size_ch, dups_log_ch, control_igg_bam_index_tuple_ch, wt_igg_bam_index_tuple_ch )

            
            // control_bw_meta_ch = mk_bw_call_peaks_workflow_sicer2.out.control_meta_bw_ch

            // wt_bw_meta_ch = mk_bw_call_peaks_workflow_sicer2.out.wt_meta_bw_ch

            control_bw_meta_ch = mk_bw_call_peaks_workflow_sicer2.out.control_bigwig

            wt_bw_meta_ch = mk_bw_call_peaks_workflow_sicer2.out.wt_bigwig

            // getting the cpm normalized bigwigs
            // control_meta_cpm_bw = mk_bw_call_peaks_workflow_sicer2.out.control_meta_cpm_bw_ch
            // wt_meta_cpm_bw = mk_bw_call_peaks_workflow_sicer2.out.wt_meta_cpm_bw_ch

            // the mk_bw_call_peaks_workflow workflow also emits the broad peaks that were called so I will grab those
            // all_broadpeaks_ch = mk_bw_call_peaks_workflow.out.macspeaks_ch

            master_peak_list_true = mk_bw_call_peaks_workflow_sicer2.out.master_peaks_list_ch
            up_peaks_list_true = mk_bw_call_peaks_workflow_sicer2.out.up_peaks_list_ch
            down_peaks_list_true = mk_bw_call_peaks_workflow_sicer2.out.down_peaks_list_ch
            unchanging_peaks_list_true = mk_bw_call_peaks_workflow_sicer2.out.unchanging_peaks_list_ch

            sicer2_concat_meta_peak = mk_bw_call_peaks_workflow_sicer2.out.sicer2_concat_meta_peaks_ch


            /////////////////////////////////////// For when creating master peaks generated in this pipeline /////////////
            // this is replicating the process that uses peaks so i can use the peaks that were made in the pipeline
            plot_signal_up_down_peaks_workflow(control_bw_meta_ch, wt_bw_meta_ch, up_peaks_list_true, down_peaks_list_true, bisulfate_bigwig_ch, master_peak_list_true, unchanging_peaks_list_true, cpg_island_unmasked_ch)

            // the grouped channel that uses only data from the pipeline not geo control
            combined_bigwig_peak_2grouped_ch = plot_signal_up_down_peaks_workflow.out.exper_rep_bigwig_peak_group

            // if i want to keep using this channel, it also has to be below the workflow (only for when using geo control data)
            combined_bigwig_meta_2grouped_ch = plot_signal_up_down_peaks_workflow.out.experiment_group_meta_cpm_ch

            bigwig_meta_ch_to_join = plot_signal_up_down_peaks_workflow.out.experiment_group_meta_to_join_ch
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            plot_diff_peaks_over_diff_genes_workflow(up_peaks_list_true, down_peaks_list_true, unchanging_peaks_list_true, master_peak_list_true, wtvslowup_genebody_ch, wtvslowdown_nochange_ch, wtvslownochange, gtf_ch, ref_genome_size_ch, knownGene_ch, proseq_up_gene_ch, proseq_down_gene_ch, proseq_unchanging_gene_ch, combined_bigwig_peak_2grouped_ch)
            
            // now to add a workflow to plot the atac-seq signal over the cut&run peaks



            //////////////////////////////////////////////////////////////////////////////////////////////////
            // need to update the below workflow to take the actual peaks generated in the pipeline
            //plot_atac_signal_over_diff_peaks_workflow(control_atac_bigwig_ch, treatment_atac_bigwig_ch, up_peaks_ch, down_peaks_ch, unchanging_peaks_ch, down_atac_peaks_ch, up_atac_peaks_ch, combined_bigwig_meta_2grouped_ch, cpg_island_unmasked_ch )
            
            // this below will be the copy of the above workflow but using actual peaks generated in pipeline
            plot_atac_signal_over_diff_peaks_workflow(control_atac_bigwig_ch, treatment_atac_bigwig_ch, up_peaks_list_true, down_peaks_list_true, unchanging_peaks_list_true, down_atac_peaks_ch, up_atac_peaks_ch, nochange_atac_peaks_ch, combined_bigwig_peak_2grouped_ch, cpg_island_unmasked_ch )
            ///////////////////////////////////////////////////////////////////////////////////////////////////


            ///////////////////////////////////////////////////////////////////////////////////////////////////
            // now we want to first find the CpG islands that are overlapping up, down, unchanging experiment(histone) peaks
            // then we plot the h1low, scrm and bisulfate signal over the new "CpG_in_up_peaks, CpG_in_down_peaks, and CpG_in_unchanging_peaks"

            //find_then_plot_cpgIslands_in_peaks_workflow(up_peaks_ch, down_peaks_ch, unchanging_peaks_ch, cpg_island_unmasked_ch, combined_bigwig_meta_2grouped_ch )

            // this is the version of the above workflow that will only use the peaks generated within the pipeline

            
            find_then_plot_cpgIslands_in_peaks_workflow(up_peaks_list_true, down_peaks_list_true, unchanging_peaks_list_true, master_peak_list_true, cpg_island_unmasked_ch, combined_bigwig_meta_2grouped_ch, combined_bigwig_peak_2grouped_ch, bigwig_meta_ch_to_join )
            ///////////////////////////////////////////////////////////////////////////////////////////////////

            // now I want to use roadmap histone data to get and enrichment barplot of which histones have more ATAC-accessability

            //get_roadmap_histone_enrichment_workflow(roadmap_broad_histones, control_atac_bigwig_ch, treatment_atac_bigwig_ch)
            
            get_roadmap_histone_enrichment_workflow(roadmap_broad_histones, roadmap_narrow_histones, sicer2_concat_meta_peak, control_atac_bam_ch, treatment_atac_bam_ch, control_proseq_bam, treatment_proseq_bam, up_peaks_list_true, down_peaks_list_true)


            ///////////////////////////////////////////////////////////////////////////////////////////////////
            // recreating the plot cpgIslands workflow but for finding which ATAC-seq peaks are in histone mark peaks

            find_then_plot_atacseqPeaks_in_experiment_peaks_workflow(up_peaks_list_true, down_peaks_list_true, unchanging_peaks_list_true, master_peak_list_true, nochange_atac_peaks_ch, up_atac_peaks_ch, roadmap_broad_histones, roadmap_narrow_histones, sicer2_concat_meta_peak, control_atac_bam_ch, treatment_atac_bam_ch)

        }
    
    
    }
}