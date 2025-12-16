



include {
    make_alignment_bw_process_control;
    make_alignment_bw_process_wt;
    plot_histone_at_genes_process;
    merge_peaks_bedtools_process;
    macs2_call_peaks_process_both;
    plot_histones_at_peaks_process;
    find_idr_in_replicates_process;
    multiqc_process;
    mk_bedgraph_process;
    seacr_peakcalls_process;
    sicer2_peakcall_process;
    mk_bed_for_sicer2_process;
    get_pval_bedgraph;
    kenttools_get_bigwig;
    // find_diff_peaks_R_process;
    plot_at_up_down_peaks_process;
    signal_over_gene_tss_process;
    bedtools_stranded_create_process;
    merge_concat_peaks_process;
    diff_peaks_intersect_diff_genes_process;
    atac_signal_over_peaks_process;
    get_CpG_islands_in_peaks_process;
    plot_over_diff_cpg_regions_process;
    atac_enrich_counts_process;
    atac_enrich_counts_process as proseq_enrich_counts_process;
    r_atac_enrich_plot_process;
    get_atacPeaks_in_roadmapPeaks_process;
    get_atacPeaks_in_genetss_process;
    merge_bams_on_condition_process;
    get_merged_bigwig_process;
    atac_enrich_counts_2nd_version_process;
    r_atac_enrich_plot_2nd_version_process;
    sicer2_peakcall_process_noigg;
    concat_sicer2_peaks_process
    // macs2_call_peaks_process_wt
    

}from '../modules/peak_analysis_modules.nf'

include {
    find_diff_peaks_R_process;
    find_diff_peaks_R_process_SE
}from '../modules/r_differential_analysis_module.nf'




workflow mk_bw_call_peaks_workflow {


    take:
    control_bams_index_tuple_ch
    wt_bams_index_tuple_ch
    ref_genome_ch
    ref_genome_size_ch
    dups_log_ch

    


    main:

    // first I need to separate the control and wt by their histone marks

    control_bams_index_tuple_ch
        .map { unique_name, bam_index_tuple -> 
        
        bam_path = bam_index_tuple[0]
        bai_path = bam_index_tuple[1]

        basename = bam_path.baseName
        tokens = basename.tokenize("_")
        
        condition_label = tokens[0]

        histone_label = tokens[1]

        replicate_label = tokens[2]

        bio_label = tokens[3]

        bam_file_name = bam_path.name

        meta_name = "${tokens[0]}_${tokens[1]}_${tokens[2]}_${tokens[3]}"



        tuple(condition_label, histone_label, replicate_label, meta_name, bam_file_name, bam_path, bai_path)//.join(','))


        }
        .groupTuple(by:1, sort: true) // how the grouping would look  control h3k27me3 [[H1low, H1low, H1low], H3k27me3, [r2, r3, r1], [H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam, H1low_H3k27me3_r3_S27_001.trim.st.all.blft.qft.rmdup.sorted.bam, H1low_H3k27me3_r1_S25_001.trim.st.all.blft.qft.rmdup.sorted.bam],
        // .multiMap{row ->

        // histone_key = row[1] // hoping the second element in the grouping list is the histone key. look above at a snippet of how grouping looked

        // result = [:]
        // if (histone_key.contains('H3k27me3')) { result.h3k27me3 = row} 

        // if (histone_key.contains('H3k9me3')) { result.h3k9me3 = row}

        // return result

        // // h3k27me3: histone_key.contains('H3k27me3') ? row : null
        // // h3k9me3: histone_key.contains('H3k9me3') ? row : null

        // }
        .set{control_bam_meta_ch} // how one of the channels look [H1low,H3k27me3,r2,H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam,/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/bin/H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam,/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/bin/H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam.bai]
        //.view()
    //control_bam_meta_ch.view { v -> "control $v"} // not using multimap. here is how it looks.  [[H1low, H1low, H1low], H3k27me3, [r2, r3, r1], [H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam, H1low_H3k27me3_r3_S27_001.trim.st.all.blft.qft.rmdup.sorted.bam, H1low_H3k27me3_r1_S25_001.trim.st.all.blft.qft.rmdup.sorted.bam],

    //control_bam_meta_ch.view()

    wt_bams_index_tuple_ch
        .map { unique_name, bam_index_tuple -> 
        
        bam_path = bam_index_tuple[0]
        bai_path = bam_index_tuple[1]

        basename = bam_path.baseName
        tokens = basename.tokenize("_")
        
        condition_label = tokens[0]

        histone_label = tokens[1]

        replicate_label = tokens[2]

        bio_label = tokens[3]

        bam_file_name = bam_path.name

        meta_name = "${tokens[0]}_${tokens[1]}_${tokens[2]}_${tokens[3]}"

        

        tuple(condition_label, histone_label, replicate_label, meta_name, bam_file_name, bam_path, bai_path)//.join(', '))


        }
        .groupTuple(by:1, sort: true) // how the grouping would look [Scrm, Scrm, Scrm], H3k27me3, [r3, r2, r1], [Scrm_H3k27me3_r3_S3_001.trim.st.all.blft.qft.rmdup.sorted.bam, Scrm_H3k27me3_r2_S2_001.trim.st.all.blft.qft.rmdup.sorted.bam, Scrm_H3k27me3_r1_S1_001.trim.st.all.blft.qft.rmdup.sorted.bam],
        //.view()
        // .multiMap{row ->

        // histone_key = row[1] // hoping the second element in the grouping list is the histone key. look above at a snippet of how grouping looked

        // result = [:]
        // if (histone_key.contains('H3k27me3')) { result.h3k27me3 = row} 

        // if (histone_key.contains('H3k9me3')) { result.h3k9me3 = row}

        // return result

        // //h3k27me3: histone_key.contains('H3k27me3') ? row : null
        // //h3k9me3: histone_key.contains('H3k9me3') ? row : null

        // }
        .set{wt_bam_meta_ch} // how one of the wt channels look with out grouping, but i ended up grouping so look above [Scrm,H3k27me3,r3,Scrm_H3k27me3_r3_S3_001.trim.st.all.blft.qft.rmdup.sorted.bam,/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/bin/Scrm_H3k27me3_r3_S3_001.trim.st.all.blft.qft.rmdup.sorted.bam,/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/bin/Scrm_H3k27me3_r3_S3_001.trim.st.all.blft.qft.rmdup.sorted.bam.bai]
    // wt_bam_meta_ch.h3k27me3.view { v -> "h3k27me3 $v"}
    // wt_bam_meta_ch.h3k9me3.view { v -> "h3k9me3 $v"}


    // first I need to make a process to make bigwig files out of each bam in the histone marks groups.
    // these bigwig files will be made from the bam files 
    // all_control_bams = control_bam_meta_ch.h3k27me3.concat(control_bam_meta_ch.h3k9me3)
    // all_control_bams.view()
    
    //control_bam_meta_ch.view()

    make_alignment_bw_process_control(control_bam_meta_ch)

    make_alignment_bw_process_wt(wt_bam_meta_ch)


    if (params.raw_bigwig) {

        control_bigwig = make_alignment_bw_process_control.out.bigwig_meta_ch
        wt_bigwig = make_alignment_bw_process_wt.out.cpm_bigwig_meta_ch

    }
    else if (params.cpm_bigwig) {

        control_bigwig = make_alignment_bw_process_control.out.cpm_bigwig_meta_ch
        wt_bigwig = make_alignment_bw_process_wt.out.cpm_bigwig_meta_ch

    }
    else if (params.rpgc_bigwig) {

        control_bigwig = make_alignment_bw_process_control.out.rpgc_bigwig_meta_ch
        wt_bigwig = make_alignment_bw_process_wt.out.rpgc_bigwig_meta_ch

    }
    else {

        control_bigwig = make_alignment_bw_process_control.out.bigwig_meta_ch
        wt_bigwig = make_alignment_bw_process_wt.out.cpm_bigwig_meta_ch

    }


    // control_meta_bw_ch = make_alignment_bw_process_control.out.bigwig_meta_ch

    // wt_meta_bw_ch = make_alignment_bw_process_wt.out.bigwig_meta_ch

    // control_meta_cpm_bw_ch = make_alignment_bw_process_control.out.cpm_bigwig_meta_ch
    // wt_meta_cpm_bw_ch = make_alignment_bw_process_wt.out.cpm_bigwig_meta_ch

    // // now outputtting the rpgc bigwigs 
    // wt_meta_rpgc_bw_ch = make_alignment_bw_process_wt.out.rpgc_bigwig_meta_ch

    // control_meta_rpgc_bw_ch = make_alignment_bw_process_control.out.rpgc_bigwig_meta_ch


    // I want to view how the meta bigwig output channel looks.
    //control_meta_bw_ch.view()

    // now that i have the meta channels where I grouped them by the histone marks, I can put it into a process to call peaks on all the files for that histone mark and another process will spawn calling peaks for the other histone marks in parallel

    // i think if i concat the control_bam_meta_ch and the wt_bam_meta_ch I can parallelize the process
    wt_bam_meta_ch
        .concat(control_bam_meta_ch)
        .transpose()
        //.view()
        .set{concat_wt_control_bam_meta_ch}
    //concat_wt_control_bam_meta_ch.view()

    ////////////////// GIVE THE USER THE OPTION TO USE MACS2 OR SICER2 ////////////////////////////////////

    if (params.macs2) {

    
        macs2_call_peaks_process_both(concat_wt_control_bam_meta_ch, ref_genome_ch) // might need ref_genome
        //macs2_call_peaks_process_wt()

        // i want to get the ppois files from the macs2 process
        ppois_files_ch = macs2_call_peaks_process_both.out.ppois_macs2_file

        get_pval_bedgraph(ppois_files_ch, ref_genome_size_ch)

        pval_bedgraph_ch = get_pval_bedgraph.out.pvalue_bedgraph_file

        chrom_size_ch = get_pval_bedgraph.out.chrom_size_file

        kenttools_get_bigwig(pval_bedgraph_ch, ref_genome_size_ch)
        //chrom_size_ch
    
    
    }
    // if (params.sicer2) {

    //     sicer2_peakcall_process(concat_wt_control_bam_meta_ch, ref_genome_ch)
    //     sicer2_peaks = sicer2_peakcall_process.out.sicer2_peak_file

    //     sicer2_peaks.view{it -> "these are the sicer2 peaks $it"}
    // }
    



    if (params.narrowPeak_data ) {

        macspeaks_ch = macs2_call_peaks_process_both.out.narrow_peaks

        macspeaks_ch_backup = macs2_call_peaks_process_both.out.narrow_peaks


        macspeaks_ch_backup
            .map{ peakpath -> 
            
            basename = peakpath.baseName
            file_name = peakpath.name

            tokens = basename.tokenize("_")

            condition = tokens[0]
            histone = tokens[1]
            replicate = tokens[2]
            bio_rep = tokens[3]

            grouping_key = "${condition}_${histone}"

            tuple(grouping_key, condition, histone, replicate, bio_rep, file_name, basename, peakpath)
            
            }
            .groupTuple(by:0, sort:true)
            //.view()
            .set{macspeak_gtuple_meta_ch}

        // split the broadpeaks so i have them grouped by their condition label


        // now i want to get the idr peaks per each replicate combination
        //broadpeak_gtuple_meta_ch.view()
        find_idr_in_replicates_process(macspeak_gtuple_meta_ch)

        // these were the peaks without 10kb merged
        concat_idr_peaks = find_idr_in_replicates_process.out.final_concat_peaks.collect()

        // now group the concat_idr_peaks without merging
        // there are only two peak files that go below becasue of idr, and that will be one idr merged peak file for each condition
        concat_idr_peaks
            .flatten()
            .map {file -> 
            
            basename = file.baseName

            file_name = file.name

            // concat_IDR_H1low_H3k27me3_r1_vs_r2_vs_r3
            tokens = basename.tokenize("_")

            condition = tokens[2]
            histone = tokens[3]

            tuple(histone, condition, file_name, file)

            
            
            }
            .groupTuple(by:0, sort: true)
            .view{it -> "these are the peaks to send to merge_concat_peaks_process: $it"}
            // H1low is first and scrambled is next but if the file is named something else you should not hard code which is which
            // example: [H3k27me3, [conditions], [concat_IDR_H1low_H3k27me3_r1_vs_r2_vs_r3_0.4_pairs.broadPeak, concat_IDR_Scrm_H3k27me3_r1_vs_r2_vs_r3_0.4_pairs.broadPeak], [/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/work/32/d1583cfe217bf08b629e50f5d2c843/concat_IDR_H1low_H3k27me3_r1_vs_r2_vs_r3_0.4_pairs.broadPeak, /lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/work/72/7a61681f524cbe7de6d547ad41a86d/concat_IDR_Scrm_H3k27me3_r1_vs_r2_vs_r3_0.4_pairs.broadPeak]]
            .set{group_concat_idr_peaks_ch}

        merge_concat_peaks_process(group_concat_idr_peaks_ch)

        atac_concat_peaks_ch = merge_concat_peaks_process.out.concat_ATAC_peak
        
        // here i dont output the merge by 10, 30 ,100 kb becasue this is for ATAC peaks. I dont think we will ever merge
        // but i also dont think any narrow peaks assay will ever be merged

        atac_concat_peaks_ch
                //.flatten()
                .map { file -> 
                
                file_basename = file.baseName
                file_name = file.name

                tokens = file_basename.tokenize("_")

                condition_label = tokens[2]
                exper_type = tokens[3]
                //merge_dist = tokens[4]
                rep_label = tokens[4]

                tuple(condition_label, exper_type, rep_label, file_name, file)

                }
                //.groupTuple(by:1, sort:true)
                .view{it -> "these are the atac/narrowpeaks concat files $it" }
                .set{group_concat_meta_peaks_ch}

        

        control_bams_index_tuple_ch
            .concat(wt_bams_index_tuple_ch)
            .map { key, tuple ->
        
            
                bam = tuple[0]
                bai = tuple[1]
                basename = bam.baseName
                file_name = bam.name

                tokens = basename.tokenize("_")

                condition = tokens[0]
                histone = tokens[1]

                bam
                //tuple(histone, condition, file_name, bam)
            
            }
            .collect()
            .set{all_bams_paths}

        if (params.PE) {
            find_diff_peaks_R_process(group_concat_meta_peaks_ch, all_bams_paths)

            diff_peaks_tuple = find_diff_peaks_R_process.out.diff_peaks_ch
            //diff_peaks_tuple.view()

            // try keeping the peaks in separate channels
            // then when I put this into a process or workflow, just check to see if the histones match
            master_peaks_list_ch = find_diff_peaks_R_process.out.master_peak_emit.collect()
            up_peaks_list_ch = find_diff_peaks_R_process.out.up_peaks_emit.collect()
            down_peaks_list_ch = find_diff_peaks_R_process.out.down_peaks_emit.collect()
            // unchanging_peaks_list_ch = find_diff_peaks_R_process.out.unchanging_peaks_emit.collect()
            unchanging_peaks_list_ch = find_diff_peaks_R_process.out.other_peaks_emit.collect()
        }
        else if (params.SE) {
            find_diff_peaks_R_process_SE(group_concat_meta_peaks_ch, all_bams_paths)

            diff_peaks_tuple = find_diff_peaks_R_process_SE.out.diff_peaks_ch
            //diff_peaks_tuple.view()

            // try keeping the peaks in separate channels
            // then when I put this into a process or workflow, just check to see if the histones match
            master_peaks_list_ch = find_diff_peaks_R_process_SE.out.master_peak_emit.collect()
            up_peaks_list_ch = find_diff_peaks_R_process_SE.out.up_peaks_emit.collect()
            down_peaks_list_ch = find_diff_peaks_R_process_SE.out.down_peaks_emit.collect()
            // unchanging_peaks_list_ch = find_diff_peaks_R_process_SE.out.unchanging_peaks_emit.collect()
            unchanging_peaks_list_ch = find_diff_peaks_R_process_SE.out.other_peaks_emit.collect()
        }
        else {

            onError:
            throw new IllegalArgumentException("User Parameter Error: Please use the parameter --PE or --SE, if you have pair end reads or single end reads respectively!")
        }
        // using the 10kb merged idr peaks
        //find_diff_peaks_R_process(group_10kb_concat_idr_peaks_ch, all_bams_paths)

        // diff_peaks_tuple = find_diff_peaks_R_process.out.diff_peaks_ch
        // //diff_peaks_tuple.view()

        // // try keeping the peaks in separate channels
        // // then when I put this into a process or workflow, just check to see if the histones match
        // master_peaks_list_ch = find_diff_peaks_R_process.out.master_peak_emit.collect()
        // up_peaks_list_ch = find_diff_peaks_R_process.out.up_peaks_emit.collect()
        // down_peaks_list_ch = find_diff_peaks_R_process.out.down_peaks_emit.collect()
        // unchanging_peaks_list_ch = find_diff_peaks_R_process.out.unchanging_peaks_emit.collect()

    }
    else if (params.broadPeak_data) {
        // now getting the channel for the broadpeaks and will have to emit it from the workflow and put into a new workflow to plot the bigwig signal onto the called broad peaks
        
        // this is what i called the channels before changing it to macspeaks...
        //broadpeaks_ch
        //broadpeaks_ch_backup
        macspeaks_ch = macs2_call_peaks_process_both.out.broad_peaks
        macspeaks_ch_backup = macs2_call_peaks_process_both.out.broad_peaks
        //broadpeaks_ch.view() // now got the peaks so time to emit this channel.

        // okay, another thing to do is to use bedtools to merge peaks by 1kb, 2kb, and 5kb to see how they look

        merge_peaks_bedtools_process(macspeaks_ch)




        macspeaks_ch_backup
            .map{ peakpath -> 
            
            basename = peakpath.baseName
            file_name = peakpath.name

            tokens = basename.tokenize("_")

            condition = tokens[0]
            histone = tokens[1]
            replicate = tokens[2]
            bio_rep = tokens[3]

            grouping_key = "${condition}_${histone}"

            tuple(grouping_key, condition, histone, replicate, bio_rep, file_name, basename, peakpath)
            
            }
            .groupTuple(by:0, sort:true)
            //.view()
            .set{broadpeak_gtuple_meta_ch}

        // split the broadpeaks so i have them grouped by their condition label


        // now i want to get the idr peaks per each replicate combination
        //broadpeak_gtuple_meta_ch.view()
        find_idr_in_replicates_process(broadpeak_gtuple_meta_ch)

        // these were the peaks without 10kb merged
        concat_idr_peaks = find_idr_in_replicates_process.out.final_concat_peaks

        
        // now with the concat peaks, I need to put them in R and get the list of up peaks and down peaks

        concat_idr_peaks
            .map {file -> 
            
            basename = file.baseName

            file_name = file.name

            // concat_IDR_H1low_H3k27me3_r1_vs_r2_vs_r3
            tokens = basename.tokenize("_")

            condition = tokens[2]
            histone = tokens[3]

            tuple(histone, condition, file_name, file)

            
            
            }
            .groupTuple(by:0, sort: true)
            //.view()
            // H1low is first and scrambled is next but if the file is named something else you should not hard code which is which
            // example: [H3k27me3, [conditions], [concat_IDR_H1low_H3k27me3_r1_vs_r2_vs_r3_0.4_pairs.broadPeak, concat_IDR_Scrm_H3k27me3_r1_vs_r2_vs_r3_0.4_pairs.broadPeak], [/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/work/32/d1583cfe217bf08b629e50f5d2c843/concat_IDR_H1low_H3k27me3_r1_vs_r2_vs_r3_0.4_pairs.broadPeak, /lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/work/72/7a61681f524cbe7de6d547ad41a86d/concat_IDR_Scrm_H3k27me3_r1_vs_r2_vs_r3_0.4_pairs.broadPeak]]
            .set{group_concat_idr_peaks_ch}

        // now I should merge these concat peaks but keep them separate and output them in the same meta channel as above

        merge_concat_peaks_process(group_concat_idr_peaks_ch)

        // not using 10kb anymore, better to use 100kb merged peaks
        // first_10kb_peakfile = merge_concat_peaks_process.out.first_10kb_merged_peak
        // second_10kb_peakfile = merge_concat_peaks_process.out.second_10kb_merged_peak

        // NOT DOING THIS ANYMORE. I merged the condition peaks in the process and that keeps peaks from chromosomes that were missing peaks
        //first_10kb_peakfile = merge_concat_peaks_process.out.first_100kb_merged_peak
        //second_10kb_peakfile = merge_concat_peaks_process.out.second_100kb_merged_peak

        concat_master_peak_10kb = merge_concat_peaks_process.out.concat_10kb_merged_peak
        concat_master_peak_30kb = merge_concat_peaks_process.out.concat_30kb_merged_peak
        concat_master_peak_100kb = merge_concat_peaks_process.out.concat_100kb_merged_peak
        //first_10kb_peakfile.view()

        // this says 10kb merged but it is really whatever is used above
        //both_10kb_peakfiles = first_10kb_peakfile.concat(second_10kb_peakfile)

        // if (params.masterPeak100kb) {
        //     concat_master_peak_100kb
        //         .map { file -> 
                
        //         file_basename = file.baseName
        //         file_name = file.name

        //         tokens = file_basename.tokenize("_")

        //         condition_label = tokens[2]
        //         exper_type = tokens[3]

        //         tuple(condition_label, exper_type, file_name, file)

        //         }
        //         .groupTuple(by:1, sort:true)
        //         .set{group_10kb_concat_idr_peaks_ch}
        // }

        //group_10kb_concat_idr_peaks_ch = merge_concat_peaks_process.out.merged_10kb_concat_peaks
        //group_10kb_concat_idr_peaks_ch.view()


        // now I also need the bam files
        // load all bams in but use the histone mark to get the correct bams
        // the bam meta ch looks like this tuple(condition_label, histone_label, replicate_label, meta_name, bam_file_name, bam_path, bai_path)
        
        control_bams_index_tuple_ch
            .concat(wt_bams_index_tuple_ch)
            .map { key, tuple ->
        
            
                bam = tuple[0]
                bai = tuple[1]
                basename = bam.baseName
                file_name = bam.name

                tokens = basename.tokenize("_")

                condition = tokens[0]
                histone = tokens[1]

                bam
                //tuple(histone, condition, file_name, bam)
            
            }
            .collect()
            // just get all bams
            // . map { bam ->
            
            // basename = bam.baseName
            // file_name = bam.name

            // tokens = basename.tokenize("_")

            // condition = tokens[0]
            // histone = tokens[1]

            // tuple(histone, condition, file_name, bam)
            
            // }
            // .groupTuple(by:0, sort:true)
            //.view()
            //.set{meta_bam_histone_group_tuple_ch}
            .set{all_bams_paths}
        
        //all_bams_paths.view()
        //meta_bam_histone_group_tuple_ch.view()
        
        //group_concat_idr_peaks_ch.view()
        
        // filtering channels 
        // h3k27me3_idr_peaks_ch = group_concat_idr_peaks_ch.filter { it[0] == 'H3k27me3' }
        // h3k27me3_bams_ch = meta_bam_histone_group_tuple_ch.filter { it[0] == 'H3k27me3' }

        // h3k27me3_idr_peaks_ch.view()
        // h3k27me3_bams_ch.view()

        if (params.masterPeak100kb) {
            concat_master_peak_100kb
                .map { file -> 
                
                file_basename = file.baseName
                file_name = file.name

                tokens = file_basename.tokenize("_")

                condition_label = tokens[2]
                exper_type = tokens[3]
                merge_dist = tokens[4]

                tuple(condition_label, exper_type, merge_dist, file_name, file)

                }
                //.groupTuple(by:1, sort:true)
                .set{group_concat_idr_peaks_ch}

            
        }
        else if (params.masterPeak10kb) {
            concat_master_peak_10kb
                .map { file -> 
                
                file_basename = file.baseName
                file_name = file.name

                tokens = file_basename.tokenize("_")

                condition_label = tokens[2]
                exper_type = tokens[3]
                merge_dist = tokens[4]

                tuple(condition_label, exper_type, merge_dist, file_name, file)

                }
                //.groupTuple(by:1, sort:true)
                .set{group_concat_idr_peaks_ch}

            
        }
        else if (params.masterPeak30kb) {
            concat_master_peak_30kb
                .map { file -> 
                
                file_basename = file.baseName
                file_name = file.name

                tokens = file_basename.tokenize("_")

                condition_label = tokens[2]
                exper_type = tokens[3]
                merge_dist = tokens[4]

                tuple(condition_label, exper_type, merge_dist, file_name, file)

                }
                //.groupTuple(by:1, sort:true)
                .set{group_concat_idr_peaks_ch}

            
        }
        else {
            onError:
            throw new IllegalArgumentException("Error: Please put one of the following parameters in CLI when running this pipeline: masterPeak10kb, masterPeak30kb, masterPeak100kb")
            
        }

    

        // find_diff_peaks_R_process(group_concat_idr_peaks_ch, all_bams_paths)

        if (params.PE) {
            find_diff_peaks_R_process(group_concat_idr_peaks_ch, all_bams_paths)

            diff_peaks_tuple = find_diff_peaks_R_process.out.diff_peaks_ch
            //diff_peaks_tuple.view()

            // try keeping the peaks in separate channels
            // then when I put this into a process or workflow, just check to see if the histones match
            master_peaks_list_ch = find_diff_peaks_R_process.out.master_peak_emit.collect()
            up_peaks_list_ch = find_diff_peaks_R_process.out.up_peaks_emit.collect()
            down_peaks_list_ch = find_diff_peaks_R_process.out.down_peaks_emit.collect()
            // unchanging_peaks_list_ch = find_diff_peaks_R_process.out.unchanging_peaks_emit.collect()
            unchanging_peaks_list_ch = find_diff_peaks_R_process.out.other_peaks_emit.collect()
        }
        else if (params.SE) {
            find_diff_peaks_R_process_SE(group_concat_idr_peaks_ch, all_bams_paths)

            diff_peaks_tuple = find_diff_peaks_R_process_SE.out.diff_peaks_ch
            //diff_peaks_tuple.view()

            // try keeping the peaks in separate channels
            // then when I put this into a process or workflow, just check to see if the histones match
            master_peaks_list_ch = find_diff_peaks_R_process_SE.out.master_peak_emit.collect()
            up_peaks_list_ch = find_diff_peaks_R_process_SE.out.up_peaks_emit.collect()
            down_peaks_list_ch = find_diff_peaks_R_process_SE.out.down_peaks_emit.collect()
            // unchanging_peaks_list_ch = find_diff_peaks_R_process_SE.out.unchanging_peaks_emit.collect()
            unchanging_peaks_list_ch = find_diff_peaks_R_process_SE.out.other_peaks_emit.collect()
        }
        else {

            onError:
            throw new IllegalArgumentException("User Parameter Error: Please use the parameter --PE or --SE, if you have pair end reads or single end reads respectively!")
        }

        // using the 10kb merged idr peaks
        //find_diff_peaks_R_process(group_10kb_concat_idr_peaks_ch, all_bams_paths)

        // diff_peaks_tuple = find_diff_peaks_R_process.out.diff_peaks_ch
        // //diff_peaks_tuple.view()

        // // try keeping the peaks in separate channels
        // // then when I put this into a process or workflow, just check to see if the histones match
        // master_peaks_list_ch = find_diff_peaks_R_process.out.master_peak_emit.collect()
        // up_peaks_list_ch = find_diff_peaks_R_process.out.up_peaks_emit.collect()
        // down_peaks_list_ch = find_diff_peaks_R_process.out.down_peaks_emit.collect()
        // unchanging_peaks_list_ch = find_diff_peaks_R_process.out.unchanging_peaks_emit.collect()

        //up_peaks_list_ch
            //.collect()
            //.view()

    }

    if (params.make_html_report = true) {
        // running multiqc on the duplicate log files

        // group by histone 
        dups_log_ch
            .flatten()
            .map { file ->
            
            file_basename = file.baseName
            tokens = file_basename.tokenize("_")
            condition = tokens[0]
            histone = tokens[1]
            
            //grouping_key  = "${condition}_${histone}"
            tuple( condition, histone, file)
            }
            .groupTuple(by:1)
            .set{dups_log_meta_ch}
        multiqc_process(dups_log_meta_ch)
    }


    // I should try seacr now //////////////
    //control_bam_meta_ch.view()
    //wt_bam_meta_ch.view()
    
    // possibly concate the two above so i can run process for getting bedgraph files in parallel
    control_bams_index_tuple_ch
        .concat(wt_bams_index_tuple_ch)
        //.view()
        .set{combined_bam_index_tuple_ch}

    //combined_bam_index_tuple_ch.view()
    // it takes the reference genome also
    mk_bedgraph_process(combined_bam_index_tuple_ch, ref_genome_size_ch)

    bedgraphs_for_seacr_ch = mk_bedgraph_process.out.bedgraph_for_seacr
    
    // then the seacr process 
    //seacr_peakcalls_process(bedgraphs_for_seacr_ch)

    // now make a idr process for seacr peaks
    // well seacr does not give back peaks that when merged with idr will total more that 20
    // and idr needs to have a post merge of 20 peaks.
    
    /* seacr_peaks = seacr_peakcalls_process.out.seacr_peaks

    seacr_peaks
        .map{ peakpath -> 
        
        basename = peakpath.baseName
        file_name = peakpath.name

        

        tokens = basename.tokenize("_")

        condition = tokens[0]
        histone = tokens[1]
        replicate = tokens[2]
        bio_rep = tokens[3]

        grouping_key = "${condition}_${histone}"

        tuple(grouping_key, condition, histone, replicate, bio_rep, file_name, basename, peakpath)
        
        }
        .groupTuple(by:0, sort:true)
        //.view()
        .set{seacrpeaks_gtuple_meta_ch}



    seacr_idr_process(seacrpeaks_gtuple_meta_ch)
    */


    // now making a process to get sicer2 peaks
    // I can use the bam files because I added bedtools in the sicer2 conda environment I made

    //combined_bam_index_tuple_ch.view()
    //sicer2_peakcall_process(combined_bam_index_tuple_ch)

    // ill just make a process to create bed files from the bams and then input it into sicer2 env that is by itself
    // mk_bed_for_sicer2_process(combined_bam_index_tuple_ch)

    // bed_for_sicer2_ch = mk_bed_for_sicer2_process.out.bed_for_sicer2

    //bed_for_sicer2_ch
        //.filter ( ~/.*H1low.*/)
        //.view()
        //.set{hlow_bed_for_sicer2_ch}
    
    //sicer2_peakcall_process(hlow_bed_for_sicer2_ch)

    emit:
    control_bigwig //control_meta_bw_ch
    wt_bigwig //wt_meta_bw_ch
    
    
    macspeaks_ch
    // control_meta_cpm_bw_ch
    // wt_meta_cpm_bw_ch
    group_concat_idr_peaks_ch
    master_peaks_list_ch
    up_peaks_list_ch
    down_peaks_list_ch
    unchanging_peaks_list_ch
    //diff_peaks_tuple
    //concat_idr_peaks

}



workflow mk_bw_call_peaks_workflow_sicer2 {


    take:
    control_bams_index_tuple_ch
    wt_bams_index_tuple_ch
    ref_genome_ch
    ref_genome_size_ch
    dups_log_ch
    control_igg_bam_index_tuple_ch
    wt_igg_bam_index_tuple_ch

    


    main:

    control_bams_index_tuple_ch
            .concat(wt_bams_index_tuple_ch)
            .map { key, tuple ->
        
            
                bam = tuple[0]
                bai = tuple[1]
                basename = bam.baseName
                file_name = bam.name

                tokens = basename.tokenize("_")

                condition = tokens[0]
                histone = tokens[1]

                bam
                //tuple(histone, condition, file_name, bam)
            
            }
            .collect()
            .view{it -> "this is the all_bams_paths channel but are the igg bams here also? they should not be: $it"}
            .set{all_bams_paths}

    // i should just put the igg both in the sicer2 process and manually find which is AO or HC
    // first I need to separate the control and wt by their histone marks

    control_bams_index_tuple_ch
        .map { unique_name, bam_index_tuple -> 
        
        bam_path = bam_index_tuple[0]
        bai_path = bam_index_tuple[1]

        basename = bam_path.baseName
        tokens = basename.tokenize("_")
        
        condition_label = tokens[0]

        histone_label = tokens[1]

        replicate_label = tokens[2]

        bio_label = tokens[3]

        bam_file_name = bam_path.name

        meta_name = "${tokens[0]}_${tokens[1]}_${tokens[2]}_${tokens[3]}"



        tuple(condition_label, histone_label, replicate_label, bio_label, bam_file_name, bam_path, bai_path)//.join(','))


        }
        .groupTuple(by:1, sort: true) // how the grouping would look  control h3k27me3 [[H1low, H1low, H1low], H3k27me3, [r2, r3, r1], [H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam, H1low_H3k27me3_r3_S27_001.trim.st.all.blft.qft.rmdup.sorted.bam, H1low_H3k27me3_r1_S25_001.trim.st.all.blft.qft.rmdup.sorted.bam],
        // .multiMap{row ->

        // histone_key = row[1] // hoping the second element in the grouping list is the histone key. look above at a snippet of how grouping looked

        // result = [:]
        // if (histone_key.contains('H3k27me3')) { result.h3k27me3 = row} 

        // if (histone_key.contains('H3k9me3')) { result.h3k9me3 = row}

        // return result

        // // h3k27me3: histone_key.contains('H3k27me3') ? row : null
        // // h3k9me3: histone_key.contains('H3k9me3') ? row : null

        // }
        .set{control_bam_meta_ch} // how one of the channels look [H1low,H3k27me3,r2,H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam,/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/bin/H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam,/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/bin/H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam.bai]
        //.view()
    //control_bam_meta_ch.view { v -> "control $v"} // not using multimap. here is how it looks.  [[H1low, H1low, H1low], H3k27me3, [r2, r3, r1], [H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam, H1low_H3k27me3_r3_S27_001.trim.st.all.blft.qft.rmdup.sorted.bam, H1low_H3k27me3_r1_S25_001.trim.st.all.blft.qft.rmdup.sorted.bam],

    //control_bam_meta_ch.view()

    wt_bams_index_tuple_ch
        .map { unique_name, bam_index_tuple -> 
        
        bam_path = bam_index_tuple[0]
        bai_path = bam_index_tuple[1]

        basename = bam_path.baseName
        tokens = basename.tokenize("_")
        
        condition_label = tokens[0]

        histone_label = tokens[1]

        replicate_label = tokens[2]

        bio_label = tokens[3]

        bam_file_name = bam_path.name

        meta_name = "${tokens[0]}_${tokens[1]}_${tokens[2]}_${tokens[3]}"

        

        tuple(condition_label, histone_label, replicate_label, bio_label, bam_file_name, bam_path, bai_path)//.join(', '))


        }
        .groupTuple(by:1, sort: true) // how the grouping would look [Scrm, Scrm, Scrm], H3k27me3, [r3, r2, r1], [Scrm_H3k27me3_r3_S3_001.trim.st.all.blft.qft.rmdup.sorted.bam, Scrm_H3k27me3_r2_S2_001.trim.st.all.blft.qft.rmdup.sorted.bam, Scrm_H3k27me3_r1_S1_001.trim.st.all.blft.qft.rmdup.sorted.bam],
        //.view()
        // .multiMap{row ->

        // histone_key = row[1] // hoping the second element in the grouping list is the histone key. look above at a snippet of how grouping looked

        // result = [:]
        // if (histone_key.contains('H3k27me3')) { result.h3k27me3 = row} 

        // if (histone_key.contains('H3k9me3')) { result.h3k9me3 = row}

        // return result

        // //h3k27me3: histone_key.contains('H3k27me3') ? row : null
        // //h3k9me3: histone_key.contains('H3k9me3') ? row : null

        // }
        .set{wt_bam_meta_ch} // how one of the wt channels look with out grouping, but i ended up grouping so look above [Scrm,H3k27me3,r3,Scrm_H3k27me3_r3_S3_001.trim.st.all.blft.qft.rmdup.sorted.bam,/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/bin/Scrm_H3k27me3_r3_S3_001.trim.st.all.blft.qft.rmdup.sorted.bam,/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/bin/Scrm_H3k27me3_r3_S3_001.trim.st.all.blft.qft.rmdup.sorted.bam.bai]
    // wt_bam_meta_ch.h3k27me3.view { v -> "h3k27me3 $v"}
    // wt_bam_meta_ch.h3k9me3.view { v -> "h3k9me3 $v"}


    // first I need to make a process to make bigwig files out of each bam in the histone marks groups.
    // these bigwig files will be made from the bam files 
    // all_control_bams = control_bam_meta_ch.h3k27me3.concat(control_bam_meta_ch.h3k9me3)
    // all_control_bams.view()
    
    //control_bam_meta_ch.view()

    make_alignment_bw_process_control(control_bam_meta_ch)

    make_alignment_bw_process_wt(wt_bam_meta_ch)


    // control_meta_bw_ch = make_alignment_bw_process_control.out.bigwig_meta_ch

    // wt_meta_bw_ch = make_alignment_bw_process_wt.out.bigwig_meta_ch

    // control_meta_cpm_bw_ch = make_alignment_bw_process_control.out.cpm_bigwig_meta_ch
    // wt_meta_cpm_bw_ch = make_alignment_bw_process_wt.out.cpm_bigwig_meta_ch

    if (params.raw_bigwig) {

        control_bigwig = make_alignment_bw_process_control.out.bigwig_meta_ch
        wt_bigwig = make_alignment_bw_process_wt.out.bigwig_meta_ch

    }
    else if (params.cpm_bigwig) {

        control_bigwig = make_alignment_bw_process_control.out.cpm_bigwig_meta_ch
        wt_bigwig = make_alignment_bw_process_wt.out.cpm_bigwig_meta_ch

    }
    else if (params.rpgc_bigwig) {

        control_bigwig = make_alignment_bw_process_control.out.rpgc_bigwig_meta_ch
        wt_bigwig = make_alignment_bw_process_wt.out.rpgc_bigwig_meta_ch

    }
    else {

        control_bigwig = make_alignment_bw_process_control.out.bigwig_meta_ch
        wt_bigwig = make_alignment_bw_process_wt.out.bigwig_meta_ch

    }

    // I want to view how the meta bigwig output channel looks.
    //control_meta_bw_ch.view()

    // now that i have the meta channels where I grouped them by the histone marks, I can put it into a process to call peaks on all the files for that histone mark and another process will spawn calling peaks for the other histone marks in parallel

    // i think if i concat the control_bam_meta_ch and the wt_bam_meta_ch I can parallelize the process
    wt_bam_meta_ch
        .concat(control_bam_meta_ch)
        .transpose()
        //.view()
        .set{concat_wt_control_bam_meta_ch}
    //concat_wt_control_bam_meta_ch.view()

    ////////////////// GIVE THE USER THE OPTION TO USE MACS2 OR SICER2 ////////////////////////////////////

    // how does this look
    control_igg_bam_index_tuple_ch
        // .view{it -> "this is the control igg bam index tuple $it"}
        .map{ key, bam_index_tuple ->

        bam_path = bam_index_tuple[0]
        bai_path = bam_index_tuple[1]
        
        basename = bam_path.baseName
        bam_file_name = bam_path.name

        tokens = basename.tokenize("_")

        condition = tokens[0]
        igg_label = tokens[1]
        tech_rep = tokens[2]
        bio_label = tokens[3]

        tuple(condition, igg_label, tech_rep, bio_label, bam_file_name, bam_path, bai_path)

        }
        .set{control_igg_meta_ch}

    wt_igg_bam_index_tuple_ch
        // .view{it -> "this is the wild type igg bam index tuple $it"}
        .map{ key, bam_index_tuple ->

        bam_path = bam_index_tuple[0]
        bai_path = bam_index_tuple[1]
        
        basename = bam_path.baseName
        bam_file_name = bam_path.name

        tokens = basename.tokenize("_")

        condition = tokens[0]
        igg_label = tokens[1]
        tech_rep = tokens[2]
        bio_label = tokens[3]

        tuple(condition, igg_label, tech_rep, bio_label, bam_file_name, bam_path, bai_path)

        }
        .set{wt_igg_meta_ch}

    
    concat_wt_control_bam_meta_ch
        .concat(control_igg_meta_ch, wt_igg_meta_ch)
        .groupTuple(by:[0,3])
        .view{it -> "this is attempting to get the igg in the correct place $it"} // this is how one looks  [Scr, [H3K27me3, H3K27me3], [mergedReps1, r1], AO.arnold, [Scr_H3K27me3_mergedReps1_AO.arnold.bam, Scr_H3K27me3_r1_AO.arnold_IgG_S19.filt_r1_r2_filt_coor_sorted_BL_filt_sort2.bam], [/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/merged_bams/Scr_H3K27me3_mergedReps1_AO.arnold.bam, /lustre/fs4/risc_lab/scratch/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/merged_bams/igg_bams/Scr_H3K27me3_r1_AO.arnold_IgG_S19.filt_r1_r2_filt_coor_sorted_BL_filt_sort2.bam], [/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/merged_bams/Scr_H3K27me3_mergedReps1_AO.arnold.bam.bai, /lustre/fs4/risc_lab/scratch/rjohnson/pipelines/merged_hera_arnold_analysis/peak_calling_workflow/merged_bams/igg_bams/Scr_H3K27me3_r1_AO.arnold_IgG_S19.filt_r1_r2_filt_coor_sorted_BL_filt_sort2.bam.bai]]
        .set{final_norm_vs_igg_bams_meta_ch}

    // I need the norm wt and norm control bam channels separate, then each replicate combined with the wt or control igg
    wt_bams_index_tuple_ch
        .combine(wt_igg_bam_index_tuple_ch)
        .map { norm_key, norm_bam_index, igg_key, igg_bam_index ->

        norm_bam = norm_bam_index[0]
        norm_index = norm_bam_index[1]

        igg_bam = igg_bam_index[0]
        igg_index = igg_bam_index[1]

        // put all bams together
        // wt_igg_bams = [norm_bam, igg_bam]

        // get the basename of the normal and igg?
        norm_basename = norm_bam.baseName
        norm_tokens = norm_basename.tokenize("_")
        norm_bam_name = norm_bam.name

        norm_condition = norm_tokens[0]
        norm_exper = norm_tokens[1]
        norm_reps = norm_tokens[2]
        norm_grouping_label = norm_tokens[3]




        norm_tuple = tuple(norm_condition, norm_exper, norm_reps, norm_grouping_label, norm_bam_name, norm_bam, norm_index)

        // now i want to make a tuple for the igg then group the two tuples somehow

        igg_basename = igg_bam.baseName
        igg_tokens = igg_basename.tokenize("_")
        igg_bam_name = igg_bam.name

        igg_condition = igg_tokens[0]
        igg_exper = igg_tokens[1]
        igg_reps = igg_tokens[2]
        igg_grouping_label = igg_tokens[3]

        igg_tuple = tuple(igg_condition, igg_exper, igg_reps, igg_grouping_label, igg_bam_name, igg_bam, igg_index)

        [norm_tuple, igg_tuple].transpose()



        }
        // .groupTuple(by:[0,3])
        // .flatten()
        // .transpose()
        // .groupTuple(by:[0,3])
        // .toList()
        // .groupTuple()
        // .flatten() // dont use flatten
        // .groupTuple(by:[0,3])
        .view{it -> "lets see how the norm wt bam index with igg wt bam index looks $it"}
        .set{norm_igg_wt_meta_tuple_ch}

    control_bams_index_tuple_ch
        .combine(control_igg_bam_index_tuple_ch)
        .map { norm_key, norm_bam_index, igg_key, igg_bam_index ->

        norm_bam = norm_bam_index[0]
        norm_index = norm_bam_index[1]

        igg_bam = igg_bam_index[0]
        igg_index = igg_bam_index[1]

        // put all bams together
        // wt_igg_bams = [norm_bam, igg_bam]

        // get the basename of the normal and igg?
        norm_basename = norm_bam.baseName
        norm_tokens = norm_basename.tokenize("_")
        norm_bam_name = norm_bam.name

        norm_condition = norm_tokens[0]
        norm_exper = norm_tokens[1]
        norm_reps = norm_tokens[2]
        norm_grouping_label = norm_tokens[3]




        norm_tuple = tuple(norm_condition, norm_exper, norm_reps, norm_grouping_label, norm_bam_name, norm_bam, norm_index)

        // now i want to make a tuple for the igg then group the two tuples somehow

        igg_basename = igg_bam.baseName
        igg_tokens = igg_basename.tokenize("_")
        igg_bam_name = igg_bam.name

        igg_condition = igg_tokens[0]
        igg_exper = igg_tokens[1]
        igg_reps = igg_tokens[2]
        igg_grouping_label = igg_tokens[3]

        igg_tuple = tuple(igg_condition, igg_exper, igg_reps, igg_grouping_label, igg_bam_name, igg_bam, igg_index)

        [norm_tuple, igg_tuple].transpose()



        }
        // .groupTuple(by:[0,3])
        // .flatten()
        // .transpose()
        // .groupTuple(by:[0,3])
        // .toList()
        // .groupTuple()
        // .flatten() // dont use flatten
        // .groupTuple(by:[0,3])
        .view{it -> "lets see how the norm control bam index with igg control bam index looks $it"}
        .set{norm_igg_control_meta_tuple_ch}
    
    norm_igg_wt_meta_tuple_ch
        .concat(norm_igg_control_meta_tuple_ch)
        .view{it -> "this is the final norm_igg wt and control meta channels concatenated: $it"}
        .set{final_norm_igg_good_meta_ch}


    if (params.sicer2) {

        // sicer2_peakcall_process(final_norm_vs_igg_bams_meta_ch, ref_genome_ch)
        sicer2_peakcall_process(final_norm_igg_good_meta_ch, ref_genome_ch)

        // seeing if sicer2 can output all the bed files without using igg
        // sicer2_peakcall_process_noigg(concat_wt_control_bam_meta_ch, ref_genome_ch)

        sicer2_peaks = sicer2_peakcall_process.out.sicer2_peak_file

        sicer2_peaks
            .map {file -> 
            
            basename = file.baseName
            file_name = file.name
            
            tokens = basename.tokenize("_")

            condition = tokens[0]
            exper_type = tokens[1]
            tech_reps = tokens[2]

            tuple(condition, exper_type, tech_reps, file_name, file )

            }
            .groupTuple(by:1)
            .view{it -> "how the sicer2 peaks meta ch looks: $it"}
            .set{sicer2_peaks_meta_ch}

        concat_sicer2_peaks_process(sicer2_peaks_meta_ch)

        concat_sicer_peaks = concat_sicer2_peaks_process.out.sicer2_master_peak


        concat_sicer_peaks
                //.flatten()
                .map { file -> 
                
                file_basename = file.baseName
                file_name = file.name

                tokens = file_basename.tokenize("_")

                // concat_label = tokens[0]
                // master_label = tokens[1]
                // condition_label = tokens[2]
                // exper_type = tokens[3]
                // //merge_dist = tokens[4]
                // // rep_label = tokens[4]

                // tuple(condition_label, exper_type, rep_label, file_name, file)

                condition_label = tokens[2]
                exper_type = tokens[3]
                // merge_dist = tokens[4]
                merge_dist = "sicer_no_merge"

                tuple(condition_label, exper_type, merge_dist, file_name, file)

                }
                //.groupTuple(by:1, sort:true)
                .view{it -> "these are the sicer2 concat files $it" }
                .set{sicer2_concat_meta_peaks_ch}
        
        if (params.PE) {
            find_diff_peaks_R_process(sicer2_concat_meta_peaks_ch, all_bams_paths)

            diff_peaks_tuple = find_diff_peaks_R_process.out.diff_peaks_ch
            //diff_peaks_tuple.view()

            // try keeping the peaks in separate channels
            // then when I put this into a process or workflow, just check to see if the histones match
            master_peaks_list_ch = find_diff_peaks_R_process.out.master_peak_emit.collect()
            up_peaks_list_ch = find_diff_peaks_R_process.out.up_peaks_emit.collect()
            down_peaks_list_ch = find_diff_peaks_R_process.out.down_peaks_emit.collect()
            // unchanging_peaks_list_ch = find_diff_peaks_R_process.out.unchanging_peaks_emit.collect()
            unchanging_peaks_list_ch = find_diff_peaks_R_process.out.other_peaks_emit.collect()

            // other_peaks_emit
            
        }
        else if (params.SE) {
            find_diff_peaks_R_process_SE(sicer2_concat_meta_peaks_ch, all_bams_paths)

            diff_peaks_tuple = find_diff_peaks_R_process_SE.out.diff_peaks_ch
            //diff_peaks_tuple.view()

            // try keeping the peaks in separate channels
            // then when I put this into a process or workflow, just check to see if the histones match
            master_peaks_list_ch = find_diff_peaks_R_process_SE.out.master_peak_emit.collect()
            up_peaks_list_ch = find_diff_peaks_R_process_SE.out.up_peaks_emit.collect()
            down_peaks_list_ch = find_diff_peaks_R_process_SE.out.down_peaks_emit.collect()
            // unchanging_peaks_list_ch = find_diff_peaks_R_process_SE.out.unchanging_peaks_emit.collect()
            unchanging_peaks_list_ch = find_diff_peaks_R_process_SE.out.other_peaks_emit.collect()
        }
        else {

            onError:
            throw new IllegalArgumentException("User Parameter Error: Please use the parameter --PE or --SE, if you have pair end reads or single end reads respectively!")
        }
        // sicer2_peaks.view{it -> "these are the sicer2 peaks $it"}
    }
    




  

    if (params.make_html_report = true) {
        // running multiqc on the duplicate log files

        // group by histone 
        dups_log_ch
            .flatten()
            .map { file ->
            
            file_basename = file.baseName
            tokens = file_basename.tokenize("_")
            condition = tokens[0]
            histone = tokens[1]
            
            //grouping_key  = "${condition}_${histone}"
            tuple( condition, histone, file)
            }
            .groupTuple(by:1)
            .set{dups_log_meta_ch}
        multiqc_process(dups_log_meta_ch)
    }


    // I should try seacr now //////////////
    //control_bam_meta_ch.view()
    //wt_bam_meta_ch.view()
    
    // possibly concate the two above so i can run process for getting bedgraph files in parallel
    control_bams_index_tuple_ch
        .concat(wt_bams_index_tuple_ch)
        //.view()
        .set{combined_bam_index_tuple_ch}

    //combined_bam_index_tuple_ch.view()
    // it takes the reference genome also
    mk_bedgraph_process(combined_bam_index_tuple_ch, ref_genome_size_ch)

    bedgraphs_for_seacr_ch = mk_bedgraph_process.out.bedgraph_for_seacr
    
    // then the seacr process 
    //seacr_peakcalls_process(bedgraphs_for_seacr_ch)

    // now make a idr process for seacr peaks
    // well seacr does not give back peaks that when merged with idr will total more that 20
    // and idr needs to have a post merge of 20 peaks.
    
    /* seacr_peaks = seacr_peakcalls_process.out.seacr_peaks

    seacr_peaks
        .map{ peakpath -> 
        
        basename = peakpath.baseName
        file_name = peakpath.name

        

        tokens = basename.tokenize("_")

        condition = tokens[0]
        histone = tokens[1]
        replicate = tokens[2]
        bio_rep = tokens[3]

        grouping_key = "${condition}_${histone}"

        tuple(grouping_key, condition, histone, replicate, bio_rep, file_name, basename, peakpath)
        
        }
        .groupTuple(by:0, sort:true)
        //.view()
        .set{seacrpeaks_gtuple_meta_ch}



    seacr_idr_process(seacrpeaks_gtuple_meta_ch)
    */


    // now making a process to get sicer2 peaks
    // I can use the bam files because I added bedtools in the sicer2 conda environment I made

    //combined_bam_index_tuple_ch.view()
    //sicer2_peakcall_process(combined_bam_index_tuple_ch)

    // ill just make a process to create bed files from the bams and then input it into sicer2 env that is by itself
    // mk_bed_for_sicer2_process(combined_bam_index_tuple_ch)

    // bed_for_sicer2_ch = mk_bed_for_sicer2_process.out.bed_for_sicer2

    //bed_for_sicer2_ch
        //.filter ( ~/.*H1low.*/)
        //.view()
        //.set{hlow_bed_for_sicer2_ch}
    
    //sicer2_peakcall_process(hlow_bed_for_sicer2_ch)

    // emit:
    // control_meta_bw_ch
    // wt_meta_bw_ch
    
    
    // macspeaks_ch
    // control_meta_cpm_bw_ch
    // wt_meta_cpm_bw_ch
    // group_concat_idr_peaks_ch
    // master_peaks_list_ch
    // up_peaks_list_ch
    // down_peaks_list_ch
    // unchanging_peaks_list_ch
    //diff_peaks_tuple
    //concat_idr_peaks

    emit:
    
    control_bigwig //control_meta_bw_ch
    wt_bigwig //wt_meta_bw_ch
    
    
    // macspeaks_ch
    // control_meta_cpm_bw_ch
    // wt_meta_cpm_bw_ch
    sicer2_concat_meta_peaks_ch
    master_peaks_list_ch
    up_peaks_list_ch
    down_peaks_list_ch
    unchanging_peaks_list_ch
}


/*
workflow get_diff_peaks_workflow {


    take:

    all_broadpeaks_ch



    main:

    //all_broadpeaks_ch.view{file -> "this is the broadpeak: $file"}




}*/

workflow plot_histone_data_workflow {


    take:
    control_bw_meta_ch
    wt_bw_meta_ch
    wtvslowup_genebody_ch
    wtvslowdown_nochange_ch




    main:

    control_bw_meta_ch
        .map {condition_label, histone_label, replicate_label, bigwig_path ->

        bigwig_name = bigwig_path.name  // this is vectorized. it got the file name for all three in the list removing the paths
        //bigwig_name2 = bigwig_path[[1]].name 
        //bigwig_name3 = bigwig_path[[2]].name 

        file_basename = bigwig_path.baseName

        // tokens = file_basename.tokenize("_")

        // condition = tokens[0]
        // histone = tokens[1]
        // replicate = tokens[2]
        // hopefully i can end with this keeping k9 and k27 separated still

        //tuple(condition_label, histone_label, replicate_label, bigwig_name, bigwig_path)
        //tuple(bigwig_name,bigwig_path)
        bigwig_path

        }
        //.sort{ a, b -> a[2] <=> b[2]}
        //.transpose()
        //.groupTuple(by:1, sort:true)
        //.view()
        .set{control_bw_meta2_ch}

    wt_bw_meta_ch
        .map {condition_label, histone_label, replicate_label, bigwig_path ->

        bigwig_name = bigwig_path.name  // this is vectorized. it got the file name for all three in the list removing the paths
        //bigwig_name2 = bigwig_path[[1]].name 
        //bigwig_name3 = bigwig_path[[2]].name 

        file_basename = bigwig_path.baseName

        // tokens = file_basename.tokenize("_")

        // condition = tokens[0]
        // histone = tokens[1]
        // replicate = tokens[2]
        // hopefully i can end with this keeping k9 and k27 separated still

        //tuple(condition_label, histone_label, replicate_label, bigwig_name, bigwig_path)
        //tuple(bigwig_name,bigwig_path)
        bigwig_path

        }
        //.sort{ a, b -> a[2] <=> b[2]}
        //.transpose()
        //.groupTuple(by:1, sort:true)
        //.view()
        .set{wt_bw_meta2_ch}
        
    wt_bw_meta2_ch
        .concat(control_bw_meta2_ch)
        .flatten()
        .map{ file ->

        file_basename = file.baseName
        file_name = file.name
        tokens = file_basename.tokenize("_")

        condition = tokens[0]
        histone = tokens[1]
        replicate = tokens[2]

        tuple(condition, histone, replicate, file_name, file)


        }
        .groupTuple(by:1)
        //.view()
        .set {experiment_group_meta_ch}
    // making a process that takes both control and wt and makes a plot. I need to give it a gene list also

    plot_histone_at_genes_process(experiment_group_meta_ch, wtvslowup_genebody_ch, wtvslowdown_nochange_ch)





}


workflow plot_histone_calledpeak_workflow {



    take:
    control_bw_meta_ch2
    wt_bw_meta_ch2
    all_broadpeaks_ch2






    main:

    // i need to make a meta channel for the broadPeaks where I have in a list the broadPeaks for each condition and histone
    
    control_bw_meta_ch2
        .concat(wt_bw_meta_ch2)
        .map { condition, histone, replicate, bigwig_paths ->

        
        bigwig_paths
        
        }
        .flatten()
        .concat(all_broadpeaks_ch2)
        .map { paths ->

        file_basenames = paths.baseName

        file_name = paths.name
        
        tokens = file_basenames.tokenize("_")

        condition = tokens[0]

        histone = tokens[1]

        replicate = tokens[2]

        // i need to get the bio rep

        bio_rep = tokens[3]

        grouping_key = "${condition}_${histone}_${replicate}_${bio_rep}"

        tuple(grouping_key, condition, histone, replicate, file_name, paths)

        }
        .groupTuple(by:0)
        //.view()
        .set{meta_bw_peak_ch}
    //meta_bw_peak_ch.view()
    plot_histones_at_peaks_process(meta_bw_peak_ch)
    
    
    // all_broadpeaks_ch2
    //     .map {peak_paths ->

    //     peak_basename = peak_paths.baseName

    //     peak_filename = peak_paths.name

    //     tokens = peak_basename.tokenize("_")

    //     condition_label = tokens[0]
    //     histone_label = tokens[1]
    //     replicate_label = tokens[2]

    //     grouping_key = "${condition_label}_${histone_label}"

    //     tuple( condition_label, histone_label, replicate_label, peak_filename, peak_paths)


    //     }
    //     .groupTuple(by:1)
    //     // .multiMap {tuple ->

    //     // histone = tuple[1]

    //     // k27me3: histone=="H3k27me3"
    //     // k9me3: histone=="H3k9me3"

    //     // }
    //     .filter { tuple ->
        
    //     //tuple(condition, histone, replicate, peak_filename, peak_basename, peak_paths)
    //     tuple[1]=="H3k27me3"
        

    //     }
    //     //.transpose()
    //     //.view()
    //     .set{broadpeak_k27_ch}
    
    // // now getting the k9 channel

    // all_broadpeaks_ch2
    //     .map {peak_paths ->

    //     peak_basename = peak_paths.baseName

    //     peak_filename = peak_paths.name

    //     tokens = peak_basename.tokenize("_")

    //     condition_label = tokens[0]
    //     histone_label = tokens[1]
    //     replicate_label = tokens[2]

    //     grouping_key = "${condition_label}_${histone_label}"

    //     tuple( condition_label, histone_label, replicate_label, peak_filename, peak_paths)


    //     }
    //     .groupTuple(by:1)
    //     // .multiMap {tuple ->

    //     // histone = tuple[1]

    //     // k27me3: histone=="H3k27me3"
    //     // k9me3: histone=="H3k9me3"

    //     // }
    //     .filter { tuple ->
        
    //     //tuple(condition, histone, replicate, peak_filename, peak_basename, peak_paths)
    //     tuple[1]=="H3k9me3"
        

    //     }
    //     //.transpose()
    //     //.view()
    //     .set{broadpeak_k9_ch} // how it looks [[Scrm, Scrm, Scrm, H1low, H1low, H1low], H3k9me3, [r1, r3, r2, r1, r3, r2], 


    //     //.set{broadpeak_meta_ch}

    

    // control_bw_meta_ch2
    //     .map {condition_label, histone_label, replicate_label, bigwig_path ->

    //     bigwig_name = bigwig_path.name  // this is vectorized. it got the file name for all three in the list removing the paths
        

    //     file_basename = bigwig_path.baseName

        
    //     bigwig_path

    //     }
        
    //     .set{control_bw_meta3_ch}

    // wt_bw_meta_ch2
    //     .map {condition_label, histone_label, replicate_label, bigwig_path ->

    //     bigwig_name = bigwig_path.name  // this is vectorized. it got the file name for all three in the list removing the paths
         

    //     file_basename = bigwig_path.baseName

    //     bigwig_path

    //     }
        
    //     .set{wt_bw_meta3_ch}
        
    // wt_bw_meta3_ch
    //     .concat(control_bw_meta3_ch)
    //     //.concat(broadpeak_meta_ch)
    //     //.view()
    //     .flatten()
    //     .map{ file ->

    //     file_basename = file.baseName
    //     file_name = file.name
    //     tokens = file_basename.tokenize("_")

    //     condition = tokens[0]
    //     histone = tokens[1]
    //     replicate = tokens[2]

    //     grouping_key = "${condition}_${histone}_${replicate}"

    //     tuple(grouping_key, condition, histone, replicate, file_name, file_basename, file)


    //     }
    //     .groupTuple(by:2)
    //     .map{grouping_key, condition, histone, replicate, file_name, file_basename, file ->

    //     tuple(condition, histone, replicate, file_name, file)

    //     }
    //     //.view()
    //     // .groupTuple(by:0)
    //     // //.transpose() 
    //     // .map { condition, histone, replicate, file_name, file_basename, file ->

    //     // tuple(file_name, file)

    //     // }
    //     // .view()
    //     .set {experiment_group_bigwigs_meta_ch}

    // // i need to separate k27 and k9 bigwigs

    // experiment_group_bigwigs_meta_ch
    //     .filter{ tuple ->

    //     histone = tuple[1]

    //     histone == "H3k27me3"
    //     }
    //     .transpose()
    //     //.view()
    //     .set{bigwig_k27_meta_ch}

    // // now for the bigwig k9

    // experiment_group_bigwigs_meta_ch
    //     .filter{ tuple ->

    //     histone = tuple[1]

    //     histone == "H3k9me3"
    //     }
    //     .transpose()
    //     //.view()
    //     .set{bigwig_k9_meta_ch}


    // // if i now concat the two channels i can hopefully be sure the first will be k27 and the second will be k9 if i put it that way

    // bigwig_k27_meta_ch
    //     .concat(bigwig_k9_meta_ch)
    //     .set{concat_bw_k27_k9_ch}

    // broadpeak_k27_ch.transpose()
    //     .concat(broadpeak_k9_ch.transpose())
    //     .set{concat_peak_k27_k9_ch}

    //concat_bw_k27_k9_ch.view()
    //concat_peak_k27_k9_ch.view()


    // see how this looks and how i can get the broad peaks to be in here

    //experiment_group_meta_ch.view()

    //plot_histones_at_peaks_process(concat_bw_k27_k9_ch, concat_peak_k27_k9_ch)


}

// this will be very similar to the plotting histone signal at genes workflow
workflow plot_signal_up_down_peaks_workflow {




    take:

    control_meta_cpm_bw
    wt_meta_cpm_bw
    up_peaks_ch
    down_peaks_ch
    bisulfate_bigwig_ch
    master_peaks_ch
    unchanging_peaks_list_true
    cpg_island_unmasked_ch



    main:

    // need to combine the two bigwig channels

    control_meta_cpm_bw
        .map {condition_label, histone_label, replicate_label, bigwig_path ->

        bigwig_name = bigwig_path.name  // this is vectorized. it got the file name for all three in the list removing the paths
        //bigwig_name2 = bigwig_path[[1]].name 
        //bigwig_name3 = bigwig_path[[2]].name 

        file_basename = bigwig_path.baseName

        // tokens = file_basename.tokenize("_")

        // condition = tokens[0]
        // histone = tokens[1]
        // replicate = tokens[2]
        // hopefully i can end with this keeping k9 and k27 separated still

        //tuple(condition_label, histone_label, replicate_label, bigwig_name, bigwig_path)
        //tuple(bigwig_name,bigwig_path)
        bigwig_path

        }
        //.sort{ a, b -> a[2] <=> b[2]}
        //.transpose()
        //.groupTuple(by:1, sort:true)
        //.view()
        .set{control_bw_cpm_meta2_ch}

    wt_meta_cpm_bw
        .map {condition_label, histone_label, replicate_label, bigwig_path ->

        bigwig_name = bigwig_path.name  // this is vectorized. it got the file name for all three in the list removing the paths
        //bigwig_name2 = bigwig_path[[1]].name 
        //bigwig_name3 = bigwig_path[[2]].name 

        file_basename = bigwig_path.baseName

        // tokens = file_basename.tokenize("_")

        // condition = tokens[0]
        // histone = tokens[1]
        // replicate = tokens[2]
        // hopefully i can end with this keeping k9 and k27 separated still

        //tuple(condition_label, histone_label, replicate_label, bigwig_name, bigwig_path)
        //tuple(bigwig_name,bigwig_path)
        bigwig_path

        }
        //.sort{ a, b -> a[2] <=> b[2]}
        //.transpose()
        //.groupTuple(by:1, sort:true)
        //.view()
        .set{wt_bw_cpm_meta2_ch}
        
    wt_bw_cpm_meta2_ch
        .concat(control_bw_cpm_meta2_ch)
        .flatten()
        .map{ file ->

        file_basename = file.baseName
        file_name = file.name
        tokens = file_basename.tokenize("_")

        condition = tokens[0]
        histone = tokens[1]
        replicate = tokens[2]

        tuple(condition, histone, replicate, file_name, file)


        }
        .groupTuple(by: [1,2])  // here i can try groupting tuple by 2 fields or columns. By the histone and then by the replicate !!!
        //.view()
        .set {experiment_group_meta_cpm_ch}

    // new group cpm ch for joining
    wt_bw_cpm_meta2_ch
        .concat(control_bw_cpm_meta2_ch)
        .flatten()
        .map{ file ->

        file_basename = file.baseName
        file_name = file.name
        tokens = file_basename.tokenize("_")

        condition = tokens[0]
        histone = tokens[1]
        replicate = tokens[2]

        tuple(condition, histone, replicate, file_name, file)


        }
        .groupTuple(by: [1,2])  // here i can try groupting tuple by 2 fields or columns. By the histone and then by the replicate !!!
        //.view()
        .set {experiment_group_meta_to_join_ch}

    // try to make a meta map of up down masterpeaks
    up_peaks_ch
        .concat(down_peaks_ch)
        .concat(unchanging_peaks_list_true)
        .concat(master_peaks_ch)
        .flatten()
        .map {peaks -> 
        
        peak_basename = peaks.baseName

        tokens = peak_basename.tokenize("_")

        diff_peak = tokens[0] // this will be up, down or master
        exper_type = tokens[1] // this will be k27me3 or k9me3 or anytype of experiment

        tuple(diff_peak, exper_type, peaks)


        }
        .groupTuple(by:1)
        //.view()
        .set{peaks_group_meta_ch}

    experiment_group_meta_to_join_ch
        //.map { histone, condition, replicate, file_name, file ->
        //tuple(histone, tuple(histone, condition, replicate, file_name, file))

        //}
        .combine(peaks_group_meta_ch, by:1)

        //.map{histone, meta_tuple, peaks ->

        //def (expr_type, conditions, replicate, file_name, files) = meta_tuple

        //tuple(expr_type, conditions, replicate, file_name, files, peaks)

        //}

        .view()
        .set{exper_rep_bigwig_peak_group}
    // now using .join to group the two channels

    // now to plot the signal at the peaks
    
    // this below is the version that used the geo data with the h1low to make masterpeaks
    //plot_at_up_down_peaks_process(experiment_group_meta_cpm_ch, up_peaks_ch, down_peaks_ch, bisulfate_bigwig_ch, master_peaks_ch, cpg_island_unmasked_ch)

    ////////// uncomment when using only master peaks generated from data only in this pipeline //////////
    plot_at_up_down_peaks_process(exper_rep_bigwig_peak_group, bisulfate_bigwig_ch, cpg_island_unmasked_ch)
    // here i can try groupting tuple by 2 fields or columns. By the histone and then by the replicate !!!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////



    emit:

    experiment_group_meta_cpm_ch

    // this is the new combined bigwig meta ch that is grouped by histone and replicate but also has the peaks files there also
    exper_rep_bigwig_peak_group

    // i need this channel
    experiment_group_meta_to_join_ch

}

workflow plot_diff_peaks_over_diff_genes_workflow {



    take:
    up_peaks_ch
    down_peaks_ch
    unchanging_peaks_ch
    master_peaks_ch
    up_genes_ch
    down_genes_ch
    nochanging_genes_ch
    gtf_ch
    ref_genome_size_ch
    knownGene_ch
    proseq_up_gene_ch
    proseq_down_gene_ch
    proseq_unchanging_gene_ch
    combined_bigwig_meta_2grouped_ch // this will have the peak type and the peaks paths if using the version that only uses data generated in the pipeline
    


    main:

    // first lets get the the strandedness of the up and down genes then use bedtools slop to get 5kb based on strand

    //bedtools_stranded_create_process(up_genes_ch, down_genes_ch, nochanging_genes_ch, gtf_ch, ref_genome_size_ch) 

    // using the proseq genes instead
    bedtools_stranded_create_process(proseq_up_gene_ch, proseq_down_gene_ch, proseq_unchanging_gene_ch, gtf_ch, ref_genome_size_ch) 

    // changed from 5kb to 20kb
    up_genes_with_20kb = bedtools_stranded_create_process.out.up_genes_20kb_stranded
    down_genes_with_20kb = bedtools_stranded_create_process.out.down_genes_20kb_stranded
    nochange_genes_with_20kb = bedtools_stranded_create_process.out.nochange_genes_20kb_stranded


    // I might just concat the up and down peak files to get a large diff peak file then plot over up and down genes
    // the combined_bigwig_meta_2grouped_ch is grouped by histone and replicate and has this format tuple(condition, histone, replicate, file_name, file)
    // this will have independent channels where in each histone, you have a single replicate with its treatment bigwig and its control bigwig. this compares the control and treatment in the same replicate and same histone
    
    //signal_over_gene_tss_process(up_peaks_ch, down_peaks_ch, up_genes_with_20kb, down_genes_with_20kb, nochange_genes_with_20kb, up_genes_ch, down_genes_ch, combined_bigwig_meta_2grouped_ch) //version for outside data
    //signal_over_gene_tss_process( up_genes_with_20kb, down_genes_with_20kb, nochange_genes_with_20kb, up_genes_ch, down_genes_ch, combined_bigwig_meta_2grouped_ch)

    // i wasnt using the 20kb files from above.
    signal_over_gene_tss_process( proseq_up_gene_ch, proseq_down_gene_ch, proseq_unchanging_gene_ch, up_genes_ch, down_genes_ch, nochanging_genes_ch, combined_bigwig_meta_2grouped_ch)

    // now make a process to find which peaks intersect the up and down genes
    
    //not using for now
    //diff_peaks_intersect_diff_genes_process(up_peaks_ch, down_peaks_ch, unchanging_peaks_ch, master_peaks_ch, up_genes_with_20kb, down_genes_with_20kb, nochange_genes_with_20kb, knownGene_ch, proseq_up_gene_ch, proseq_down_gene_ch, proseq_unchanging_gene_ch)



}


workflow plot_atac_signal_over_diff_peaks_workflow {


    take:

    control_atac_bigwig

    treatment_atac_bigwig

    up_peaks

    down_peaks

    unchanging_peaks

    // new additions
    down_atac_peaks_ch
    up_atac_peaks_ch
    nochange_atac_peaks_ch
    combined_bigwig_meta_2grouped_ch
    cpg_island_unmasked_ch

    main:

    // now make a process that will plot this information

    // with the histone meta bigwig channels, just merge the three replicates into one for each condition in the process using cat

    // version using peaks from outside pipeline
    // atac_signal_over_peaks_process(control_atac_bigwig, treatment_atac_bigwig, up_peaks, down_peaks, unchanging_peaks,  down_atac_peaks_ch, up_atac_peaks_ch, combined_bigwig_meta_2grouped_ch, cpg_island_unmasked_ch)

    // version using peaks from the pipeline
    atac_signal_over_peaks_process(control_atac_bigwig, treatment_atac_bigwig, down_atac_peaks_ch, up_atac_peaks_ch, nochange_atac_peaks_ch, combined_bigwig_meta_2grouped_ch, cpg_island_unmasked_ch)



    //emit:
}

workflow find_then_plot_cpgIslands_in_peaks_workflow {


    take:

    up_peaks_ch
    down_peaks_ch
    unchanging_peaks_ch
    master_peak_list_true
    cpg_island_unmasked_ch
    combined_bigwig_meta_2grouped_ch
    combined_bigwig_peak_2grouped_ch

    // testing this
    bigwig_meta_ch_to_join




    main:

    // now first lets find which CpG islands overlap with the different peaks

    // have to add the grouped peak bigwig channel so i can get the peaks for up down and unchanging
    //get_CpG_islands_in_peaks_process(up_peaks_ch, down_peaks_ch, unchanging_peaks_ch, cpg_island_unmasked_ch, combined_bigwig_meta_2grouped_ch)

    //cpg_island_unmasked_ch.view(file -> "this is the cpg_island unmasked file: $file")
    
    // i just need all the true peaks in one channel. do not need the combined peak meta ch here

    up_peaks_ch
        .concat(down_peaks_ch, unchanging_peaks_ch, master_peak_list_true)
        .flatten()
        .map{ peaks -> 
        
        basename = peaks.baseName

        tokens = basename.tokenize("_")

        peak_type = tokens[0]
        exper_type = tokens[1]

        tuple(peak_type, exper_type, peaks)

        }
        .groupTuple(by:1)
        .set{all_true_peaks_ch}
    

    get_CpG_islands_in_peaks_process( cpg_island_unmasked_ch, all_true_peaks_ch)

    cpg_up_ch = get_CpG_islands_in_peaks_process.out.cpg_up_regions
    cpg_down_ch = get_CpG_islands_in_peaks_process.out.cpg_down_regions
    cpg_unchanging_ch = get_CpG_islands_in_peaks_process.out.cpg_unchanging_regions
    cpg_master_ch = get_CpG_islands_in_peaks_process.out.cpg_masterpeak_regions

    // now to plot the signal over these up, down, unchanging cpg regions
    // i will have to make a new version where i merge the cpg_types with the combined bigiwig meta ch like i did with the normal peak types
    //plot_over_diff_cpg_regions_process(cpg_up_ch, cpg_down_ch, cpg_unchanging_ch, combined_bigwig_meta_2grouped_ch)

    //cpg_up_ch.view(file -> "this should be a single up cpg file for each exper: $file")

    cpg_up_ch
        .concat(cpg_down_ch, cpg_unchanging_ch, cpg_master_ch)
        //.view()
        .map {cpg_type ->

        basename = cpg_type.baseName

        tokens = basename.tokenize("_")

        type = tokens[0]
        experiment = tokens[1]
        

        tuple(cpg_type, experiment)



        }
        .groupTuple(by:[1])
        //.view(tuple -> "this is the channel so far for cpg type: $tuple")
        .set{changing_cpg_grouped_meta_ch}
        //.combine(combined_bigwig_meta_2grouped_ch, by:0)
        //.view(meta_ch -> "this contains 3 replicates of experiments grouped with the cpgs that are changing: $meta_ch")

    bigwig_meta_ch_to_join
        .combine(changing_cpg_grouped_meta_ch, by:[1])
        //.groupTuple(by:2)
        //.view(new_tuple -> "this should be the bigwig channel grouped with the changing cpg peaks $new_tuple")
        .set{bigwig_cpg_changing_peaks_mets_ch}

    plot_over_diff_cpg_regions_process( bigwig_cpg_changing_peaks_mets_ch)


}

workflow get_proximal_distal_atac_peaks_workflow {


    take:
    up_peaks_list_ch
    down_peaks_list_ch
    unchanging_peaks_list_ch
    master_peak_list_true
    proseq_up_gene_ch
    proseq_down_gene_ch
    proseq_unchanging_gene_ch
    ref_genome_size_ch

    // I want to get the histone (H3k27me3) bigwigs in here for the calculations

    control_histone_bams
    treatment_histone_bams



    main:

    up_peaks_list_ch
        .concat(down_peaks_list_ch, unchanging_peaks_list_ch, master_peak_list_true)
        .flatten()
        .map{ peaks -> 
        
        basename = peaks.baseName

        tokens = basename.tokenize("_")

        peak_type = tokens[0]
        exper_type = tokens[1]

        tuple(peak_type, exper_type, peaks)

        }
        //.groupTuple(by:1)
        .view{it -> "this is the atac peaks that should be flattened without grouping $it"}
        .set{all_true_atac_peaks_ch}
    
    // what i need to do is get the peaks and use bedtools intersect to subset for the atac-seq peaks that are in these diff genes
    // get_atacPeaks_in_genetss_process(all_true_atac_peaks_ch, proseq_up_gene_ch, proseq_down_gene_ch, proseq_unchanging_gene_ch, ref_genome_size_ch)

    // proximal_up_atac = get_atacPeaks_in_genetss_process.out.diff_peaks_up_genes
    // proximal_down_atac = get_atacPeaks_in_genetss_process.out.diff_peaks_down_genes
    // proximal_unchanging_atac = get_atacPeaks_in_genetss_process.out.diff_peaks_unchanging_genes

    // I think what I should do now is get the differential proximal ATAC-peaks and see the fold change of H3k27me3 scrambled and H1low signal in them
    if (params.atac_analysis ) {

        get_atacPeaks_in_genetss_process(all_true_atac_peaks_ch, proseq_up_gene_ch, proseq_down_gene_ch, proseq_unchanging_gene_ch, ref_genome_size_ch)

        proximal_up_atac = get_atacPeaks_in_genetss_process.out.diff_peaks_up_genes
        proximal_up_atac.view{it -> "proximal up atac peaks in genes: $it"}

        proximal_down_atac = get_atacPeaks_in_genetss_process.out.diff_peaks_down_genes
        proximal_unchanging_atac = get_atacPeaks_in_genetss_process.out.diff_peaks_unchanging_genes

        // now using the up atac peaks anymore
        // proximal_up_atac
        //     .concat(proximal_down_atac, proximal_unchanging_atac)
        proximal_down_atac
            .concat( proximal_unchanging_atac)
            .toList()
            .map {file ->
            
            basename = file.baseName
            name = file.name

            tuple(name, basename, file)


            }
            .view{it -> "these are the atac analysis pipeline peaks $it"}
            .set{proximal_atac_peaks_meta_ch}

        // now I will have to merge the h1low bams and scrm bams to get an all merge bam
        // making a samtools process to merge the histone experiment bams to get one merged bam per condition

        // make a meta channel so I can grab which histone or experiment type I am working with
        control_histone_bams
            .map { tuple_key, bam_index_tuple -> 

            key_tokens = tuple_key.tokenize("_")
            exper_type = key_tokens[1]

            bam = bam_index_tuple[0]
            index = bam_index_tuple[1]
            
            basename = bam.baseName

            tokens = basename.tokenize("_")

            condition = tokens[0]
            //exper_type = tokens[1]

            tuple(exper_type, condition, bam)

            }
            .groupTuple(by:0)
            .set{control_histone_meta_bams}
        
        treatment_histone_bams
            //.view{it -> "this is the treatment histone bam: $it"}
            .map { tuple_key, bam_index_tuple ->

            key_tokens = tuple_key.tokenize("_")
            exper_type = key_tokens[1]

            bam = bam_index_tuple[0]
            index = bam_index_tuple[1]

            basename = bam.baseName

            tokens = basename.tokenize("_")

            condition = tokens[0]
            //exper_type = tokens[1]

            tuple(exper_type, condition, bam)

            }
            .groupTuple(by:0)
            .view{it -> "this is the treatment histone tuple after grouping: $it"}
            .set{treatment_histone_meta_bams}

        merge_bams_on_condition_process(control_histone_meta_bams, treatment_histone_meta_bams)
        
        allmerge_control_tuple = merge_bams_on_condition_process.out.merged_control_bam_index_tuple
        allmerge_treatment_tuple = merge_bams_on_condition_process.out.merged_treatment_bam_index_tuple

        both_all_merge_tuple = allmerge_control_tuple.concat(allmerge_treatment_tuple)
        // now to get the bigwig of each of these bams
        get_merged_bigwig_process(both_all_merge_tuple)

        allmerged_bigwigs = get_merged_bigwig_process.out.allmerged_bigwig.collect()

        // now to do the merging of bams channel and get the histone enrichment log2FC

        allmerged_bigwigs.view{it -> "these are the all merged bigwigs $it"}
        
        // not doing this right now

        atac_enrich_counts_2nd_version_process(proximal_atac_peaks_meta_ch, both_all_merge_tuple)

        enrich_counts_tab_ch = atac_enrich_counts_2nd_version_process.out.raw_enrichment_counts.collect()

        enrich_counts_tab_ch
            .flatten()
            .map{ file ->
            // i need to get the tokens and find which used broad and narrow then group them by that, so i have scrm with h1low broad, and scrm with h1low narrow
            basename = file.baseName

            filename = file.name

            tokens = basename.tokenize("_")

            condition = tokens[0]
            peak_type = tokens[1]

            tuple(peak_type, basename, filename, file)

            }
            .groupTuple(by:0)
            .view(tuple -> "This is the meta channel for enrichment counts tab file grouped by peak type: $tuple")
            .set{enrich_counts_tab_meta_ch}

        r_atac_enrich_plot_2nd_version_process(enrich_counts_tab_meta_ch)

    }
    else {

    }

    //emit:
}


workflow get_roadmap_histone_enrichment_workflow {


    take:
    roadmap_broad_histones
    roadmap_narrow_histones
    idr_merged_peaks
    control_atac_bigwig_ch // these are actually bam files not bigwig files 
    treatment_atac_bigwig_ch
    control_proseq_bam
    treatment_proseq_bam
    up_peaks_list_true
    down_peaks_list_true



    main:

    // checking how the idr_merged_peaks look
    idr_merged_peaks
        //.view{peaks -> "these should be the idr merged peaks for each condition: $peaks"}
        .map {condition, exper_type, merge_dist, file_name, path ->
        
        path
        }
        .concat( down_peaks_list_true) // not using the up peaks here???
        // .concat(up_peaks_list_true, down_peaks_list_true) // new part to test
        .flatten()
        .toList()
        //.view{path -> "how does flatten then toList look? $path"}
        .map{merged_peak_path ->

        basename = merged_peak_path.baseName
        name = merged_peak_path.name

        tuple(name, basename, merged_peak_path)

        }
        .view{peak_tuple -> "this is the new peak tuple with combined up and down peaks for the in pipeline peaks: $peak_tuple"}
        .set{pipeline_merged_idr_peak_meta_ch}

    // up_peaks_list_true  // .view{it -> "see true up peaks: $it"}
    //     .combine(down_peaks_list_true)
    //     .flatten()
    //     .map { file ->

    //     basename = file.baseName

    //     name = file.name

    //     tokens = basename.tokenize("_")

    //     peak_type = tokens[0]

    //     exper_type = tokens[1]

    //     tuple(name, basename, file)

    //     }
    //     .view{it -> "true changing peaks meta channel: $it"}
    //     .set{changing_peaks_meta_ch}
        

    // now I should combine the master and changing peaks into one channel
    // pipeline_merged_idr_peak_meta_ch
    //     .concat(changing_peaks_meta_ch)
    //     .view{it -> "combining two ch the changing and master peaks: $it"}
    //     .set{final_pipeline_peaks_meta_ch}

    roadmap_broad_histones
        .map {file ->
        
        basename = file.baseName
        name = file.name

        tuple(name, basename, file)


        }
        .view{it -> "this is how the broad roadmap histone channel looks: $it"}
        .set{roadmap_broad_histone_meta_ch}

    roadmap_narrow_histones
        .map {file ->
        
        basename = file.baseName
        name = file.name

        tuple(name, basename, file)


        }
        //.view()
        .set{roadmap_narrowhistone_meta_ch}
    // now i should combine the two atac bigwig files then flatten them so they are run in parallel

    // let's just concate all the meta peaks into one channel and use that
    pipeline_merged_idr_peak_meta_ch
        .concat(roadmap_broad_histone_meta_ch, roadmap_narrowhistone_meta_ch)
        .set{ all_peaks_roadmap_pipeline_meta_ch}

    control_atac_bigwig_ch
        .concat(treatment_atac_bigwig_ch)
        //.flatten()
        // .map{ file ->
        
        // basename = file.baseName

        // tokens = basename.tokenize("_")

        // condition_label = tokens[0]
        // experiment_label = tokens[1]
        // replicate_label = tokens[2]
        // tuple(condition_label, experiment_label, replicate_label)
        // }
        .map { key, bam_bai_list ->
            
            bam_bai_list

        }
        .view{it -> "this should be the atac bigwig/ bams : $it"}
        .set{atac_bigwig_cat}

    control_proseq_bam
        .concat(treatment_proseq_bam)

        .view{it -> "Proseq bams concatenated: $it"}
        .set{proseq_bams_cat}

    proseq_bams_cat
        .concat(atac_bigwig_cat) // this is actually the atac bam files
        .set{bam_files_to_plot}

    bam_files_to_plot.view{it -> "Combined proseq and atac bams: $it"}
    // now i want to make a channel that will have the two or more different experiments bams but they will go in parrallel for each experiment 
    
    // now make a deeptools process that uses multibigwigsummary to get the counts

    
    // for this process to work, the atac-seq bam files must only have one field in its name. scrm-allmerge.bam as an example
    // while the other bam files must be in the normal format.
    atac_enrich_counts_process(roadmap_broad_histone_meta_ch, roadmap_narrowhistone_meta_ch, pipeline_merged_idr_peak_meta_ch, atac_bigwig_cat)

    proseq_enrich_counts_process(roadmap_broad_histone_meta_ch, roadmap_narrowhistone_meta_ch, pipeline_merged_idr_peak_meta_ch, proseq_bams_cat)

    enrich_counts_tab_ch = atac_enrich_counts_process.out.raw_enrichment_counts.collect()

    enrich_proseq_counts_tab_ch = proseq_enrich_counts_process.out.raw_enrichment_counts.collect()

    enrich_counts_tab_ch
        .flatten()
        .map{ file ->
        // i need to get the tokens and find which used broad and narrow then group them by that, so i have scrm with h1low broad, and scrm with h1low narrow
        basename = file.baseName

        filename = file.name

        tokens = basename.tokenize("_")

        condition = tokens[0]
        peak_type = tokens[3]

        tuple(peak_type, basename, filename, file)

        }
        .groupTuple(by:0)
        .view(tuple -> "This is the meta channel for enrichment counts tab file grouped by peak type: $tuple")
        .set{enrich_counts_tab_meta_ch}

    enrich_proseq_counts_tab_ch
        .flatten()
        .map{ file ->
        // i need to get the tokens and find which used broad and narrow then group them by that, so i have scrm with h1low broad, and scrm with h1low narrow
        basename = file.baseName

        filename = file.name

        tokens = basename.tokenize("_")

        condition = tokens[0]
        peak_type = tokens[4] // this will be broad or narrow peak epigenetic states

        tuple(peak_type, basename, filename, file)

        }
        .groupTuple(by:0)
        .view(tuple -> "This is the meta channel for proseq enrichment counts tab file grouped by peak type: $tuple")
        .set{enrich_proseq_counts_tab_meta_ch}

    // then maybe I can concat both channels now
    enrich_counts_tab_meta_ch // this is the atac-seq channel
        .concat(enrich_proseq_counts_tab_meta_ch) // this is the proseq channel
        .set{final_enrichment_counts_meta_ch}

    // now I want to make a process that will take the output tab files from above and make the enrichemnt in R

    // r_atac_enrich_plot_process(enrich_counts_tab_meta_ch)
    r_atac_enrich_plot_process(final_enrichment_counts_meta_ch)

    //emit:
}


workflow find_then_plot_atacseqPeaks_in_experiment_peaks_workflow {


    take:

    up_peaks_list_true
    down_peaks_list_true
    unchanging_peaks_list_true
    master_peak_list_true
    nochange_atac_peaks_ch
    up_atac_peaks_ch
    roadmap_broad_histones
    roadmap_narrow_histones
    idr_merged_peaks
    control_atac_bam_ch
    treatment_atac_bam_ch



    main:

    roadmap_broad_histones
        .map {file ->
        
        basename = file.baseName
        name = file.name


        tuple(name, basename, file)


        }
        //.view()
        .set{roadmap_broad_histone_meta_ch}

    roadmap_narrow_histones
        .map {file ->
        
        basename = file.baseName
        name = file.name

        tuple(name, basename, file)


        }
        //.view()
        .set{roadmap_narrowhistone_meta_ch}

    // first make a process that finds which atac peaks are in the (up down unchanging peaks maybe) and in the roadmap peaks 
    // I will start with the roadmap and then hopefull go on to the idr_merged_peaks instead

    get_atacPeaks_in_roadmapPeaks_process(nochange_atac_peaks_ch, up_atac_peaks_ch, roadmap_broad_histones.flatten(), roadmap_narrow_histones.flatten() )

    // so i am actually doing percent overlap, which means, I now take the number of  up atac peaks in a histone (ex H3k27me3), and divide that by the total number of up atac peaks in general
    // do that for all up and unchanging atac peaks in histones
    // then divide the percentages up/unchanging and take the log2 of that number to see which is increasing or decreasing

    // i will make an r process below to do that after I figure it out in rstudio


    //emit:
}