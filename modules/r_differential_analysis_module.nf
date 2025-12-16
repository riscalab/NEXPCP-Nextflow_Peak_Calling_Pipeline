
process find_diff_peaks_R_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/r_language'


    publishDir "./nextflow_R_script_outputs/${idr_histone}/", mode: 'copy', pattern: '*'
    // if (params.narrowPeak_data) {
    //     publishDir "./nextflow_R_script_outputs/narrow_Peaks_called/${idr_histone}/", mode: 'copy', pattern: '*'
    // }
    // else {
    //     publishDir "./nextflow_R_script_outputs/broad_Peaks_called/${idr_histone}/", mode: 'copy', pattern: '*'
    // }
    
    label 'normal_big_resources'

    //debug true

    input:

    tuple val(idr_condition), val(idr_histone), val(merge_dist), val(idr_peak_name), path(idr_peak_path)

    // these are all experiment bams.
    // i need to use the histone label, and condition label to get the correct bams in this process
    path(bam_path_list)
    


    
    output:

    path("*.{png,pdf,bed}"), emit: all_r_plots

    path("${master_peak_export_out}"), emit: master_peak_emit

    tuple val("${full_condition}"), val("${idr_histone}"), path("${up_peaks_out}"), path("${down_peaks_out}"), path("${unchanging_peaks_out}"), emit: diff_peaks_ch 

    path("${up_peaks_out}"), emit: up_peaks_emit
    path("${down_peaks_out}"), emit: down_peaks_emit
    path("${unchanging_peaks_out}"), emit: unchanging_peaks_emit
    path("${other_peaks_out}"), emit: other_peaks_emit



    script:

    //full_condition = "${idr_condition[0]}vs${idr_condition[1]}"
    full_condition = "${idr_condition}"

    //first_idr_peak = idr_peak_name[0]
    //second_idr_peak = idr_peak_name[1]

    bam_name_list = bam_path_list.collect { "\"${it.getName()}\"" }.join(',')

    //peaks_list = [first_idr_peak, second_idr_peak]

    // need the output file name for masterpeak export

    if (params.narrowPeak_data ) {

        merge_dist = 'ATAC_or_narrow_peak'
    }

    master_peak_export_out = "masterpeak_${idr_histone}_${merge_dist}_maxgap_${idr_condition}.bed"

    up_peaks_out = "up_${idr_histone}_${full_condition}_regulated_peaks.bed"
    down_peaks_out = "down_${idr_histone}_${full_condition}_regulated_peaks.bed"
    unchanging_peaks_out = "unchanging_${idr_histone}_${full_condition}_regulated_peaks.bed"
    other_peaks_out = "others_${idr_histone}_${full_condition}_regulated_peaks.bed"

    // masterPeak100kb = false
    // masterPeak10kb = false
    // masterPeak30kb = false

    resize_num = params.masterPeak100kb ? 100000 : (params.masterPeak10kb ? 10000: (params.masterPeak30kb ? 30000 : 1))

    // if (params.masterPeak100kb) {
    //     resize_num = 100000
    // }
    // if (params.masterPeak10kb) {
    //     resize_num = 10000
    // }
    // if (params.masterPeak30kb) {
    //     resize_num = 30000
    // }
    // if (!params.masterPeak100kb & !params.masterPeak10kb & !masterPeak30kb) {
    //     resize_num = 1
    // }

    if (params.narrowPeak_data || params.sicer2) {



        """
        #!/usr/bin/env Rscript

        bam_list = c(${bam_name_list})

        print(bam_list)

        #peaklist = list("\${first_idr_peak}", "\${second_idr_peak}")
        #peaklist = \${peaks_list}

        #print(peaklist)

        #best_hlow_idr
        #best_scrm_idr

        library(DESeq2)
        library(GenomicRanges)
        library(chromVAR)
        library(tidyr)
        library(EnhancedVolcano)
        library(readr)

        print(dir())
        # making the master peak genomic ranges object
        
        #mPeak = GRanges()

        #for (peakfile in c("./\${first_idr_peak}", "./\${second_idr_peak}")) {
        
            #peaktable = read.table(peakfile, header = FALSE, sep = "\t")
        
            #gr_object = GRanges(seqnames = peaktable\$V1, IRanges(start = peaktable\$V2, end = peaktable\$V3), strand = "*" )
        
            #gr_object = rtracklayer::import(peakfile, format = "BED")

            #mPeak = append(mPeak, gr_object)
        #}

        #peaktable = read.table("\${idr_peak_name}", header = FALSE)
        #mPeak = GRanges(seqnames = peaktable\$V1, IRanges(start = peaktable\$V2, end = peaktable\$V3), strand = "*" )

        # might have to let it auto detect the file type when using narrow peaks that arent merged into a strict bed file
        mPeak = rtracklayer::import("${idr_peak_name}", format = "BED")

        # making sure there are no redundant peaks
        # this will allow me to resizd all of the regions by 50kb on each side, so any regions that are less than 100kb away from eachother will be merged by reduce because they overlap
        
        
        #masterPeak2 = reduce(mPeak)

        # using resize to extnd regions
        extended_master_peak = resize(mPeak, width = width(mPeak)+1, fix = "center")

        # before trimming i need to add the sequence lengths of the correct genome
        # will have to find a way to automate this step for other genomes

        library(BSgenome.Hsapiens.UCSC.hg38)
        #seq_lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
        #seqlengths(extended_master_peak) <- seq_lengths[seqlevels(extended_master_peak)]
        #trimmed_master_peak = trim(extended_master_peak)
        #masterPeak_beta = reduce(trimmed_master_peak)

        masterPeak_beta = reduce(extended_master_peak)
        # now to resize again to remove the artificial 50kb from the start and end of the newly reduced peak
        
        new_peak_size = width(masterPeak_beta)-1
        masterPeak = resize(masterPeak_beta,  width = new_peak_size , fix = "center")
        
        #print(masterPeak)

        # now I want to keep the standard chromosomes
        #seqnames_to_keep =masterPeak@seqnames@values[1:23]

        #masterPeak = masterPeak[seqnames(masterPeak) %in% seqnames_to_keep]

        # hoping to export my GRanges object master peaks to a bed file
        #write_tsv("\${master_peak_export_out}")
        rtracklayer::export.bed(masterPeak, con = "${master_peak_export_out}")

        

        list_of_mBams = c()
        
        for (bam in bam_list) {
    
            tokens = strsplit(bam, split = "_")[[1]]
            #print(tokens)
            
            histone = tokens[2]
            #print(histone)

            #list_of_mBams = list()

            if ("${idr_histone}" == histone) {
                #print(bam)
                list_of_mBams = c(list_of_mBams, bam)

            }
        }

        print(list_of_mBams)
    
        
        # now doing the next section to make the count matrix and the condition and experiment design

        # making the matrix 

        countsMatrix = matrix(data = NA, length(masterPeak), length(list_of_mBams) )


        # I want to get the chr start end of the peaks and put them in the matrix as column names so I know which peaks I am looking at from the deseq2 results in the differential peaks.
        seqnames(masterPeak)

        # i should put the unique ids from the master peak object in the matrix as row names.
        unique_masterPeak_ids = paste0(as.character(seqnames(masterPeak)), ":" ,start(masterPeak), "-", end(masterPeak))

        rownames(countsMatrix) = unique_masterPeak_ids

        ####################### new version ##########################
        library(Rsubread)
        library(dplyr)

        # featureCounts in rsubread needs a dataframe and in the correct format, with a GeneID column first then chr1 column
        df_masterPeak = as.data.frame(masterPeak)

        # now I need to change the second column name to chr
        df_masterPeak = df_masterPeak %>% rename(chr = seqnames)

        # next i need to add the unique ids as the columns in this dataframe
        df_masterPeak\$GeneID = unique_masterPeak_ids
        ##############################################################

        #getting the list of bam base names to add to the matrix column names
        list_bam_basenames = list()


        # for deseq2 i need the condition design. so hlow and scrm
        # then i will find a way to tally how many of each are there so i can automate the rep count
        condition_design = list()

        type_design = list()

        for (x in c(1:length(list_of_mBams))) {
        
        
        path_bam = list_of_mBams[[x]]
        print(path_bam)
        bam_basename = basename(as.character(path_bam))
        bam_tokens = strsplit(basename(as.character(path_bam)), split = "_")[[1]]
        
        # labeling the important tokens so it is easier to keep track of
        condition = bam_tokens[1]
        histone = bam_tokens[2]
        replicate = bam_tokens[3]
        
        
        # for later parts I need the condition and to know how many times to repeat them
        condition_design = append(condition_design, paste(condition,histone,sep="_"))
        
        
        # also get the replicate too for type design
        type_design = append(type_design, replicate)
        #type_design = append(type_design, paste(histone,replicate, sep="_"))
        
        
        
        # using chromVAR getCounts function
        # not using chromVAR anymore, for now
        #fragment_counts = getCounts(path_bam, masterPeak, paired= TRUE, by_rg = FALSE, format = "bam")

        fragment_counts <- featureCounts(
            files = path_bam,
            annot.ext = df_masterPeak,
            isPairedEnd = TRUE,      # TRUE if your data were paired-end
            strandSpecific = 0,       # 0 = unstranded, 1/2 if strand-specific
            useMetaFeatures = TRUE
        )
        
        
        # putting the fragment counts in the column labeled by the bam name
        countsMatrix[,x] = fragment_counts\$counts
        
        list_bam_basenames = append(list_bam_basenames, bam_basename)
        
        }

        colnames(countsMatrix) = list_bam_basenames

        ## this part i would have to tell the use to name their conditions as control and treatment so it can make the design treatment vs control

        #first removing the low count fragments. only keeping the rows that are above 5 count in total

        keep_Rows = which(rowSums(countsMatrix) > 5)

        filt_countmatrix = countsMatrix[keep_Rows,]

        # now to get the condition_design
        condition_counts = table(unlist(condition_design))

        # this gives back the names and the counts give back the counts for each of the names
        condition_names = names(condition_counts)
        condition_num = as.numeric(condition_counts)




        # now i can put both lists in to get back the experiment design
        #condition_factor = factor(rep(condition_names, times=condition_num))

        # Build condition factor in column order
        condition_factor <- factor(unlist(condition_design))

        # I want it so the treatment is second in the list here so it is first in the experimental design in deseq2. This way the experimental design will have the diff results going up or down in the treatment vs control

        # for treatment I need to make nextflow take the users input so they can relevel this section
        #treatment = "H1low_H3k27me3"
        treatment = "${params.treatment_name}_${idr_histone}"

        # this doesnt work because you need to put ref = levels(condition_factor)[2] to assign the actual relevel
        #if (levels(condition_factor)[1] == treatment ) {
        #condition_factor = relevel(condition_factor, levels(condition_factor)[2])
        #}else {
        #condition_factor
        #}

        # have to do this instead to make the treatment always the baseline
        #condition_factor <- relevel(condition_factor, ref = treatment)

        levels_now <- levels(condition_factor)
        control <- levels_now[ levels_now != treatment ]  # the other level

        # Set control as reference → LFC = treatment / control
        condition_factor <- relevel(condition_factor, ref = control)


        print(condition_factor)

        # repeating the above to have another column with type (replicates)
        #type_counts = table(unlist(type_design))
        #type_names = names(type_counts)
        #type_num = as.numeric(type_counts)

        #type_factor = factor(rep(type_names, times=type_counts))
        #type_factor = factor(rep(type_names, times=type_num[1]))

        # Build type factor in column order
        type_factor <- factor(unlist(type_design))

        # I want to get the idr threshold and use that as input for the file names and other things

        #peak_file = basename(peaklist[[1]])
        peak_file = basename("./${idr_peak_name}")

        idr_used = strsplit(peak_file, split = "_")[[1]][10]
        idr_used

        # now for deseq2 workflow


        # now to do the normal deseq2 workflow

        dds = DESeqDataSetFromMatrix(countData = filt_countmatrix,
                                    colData = DataFrame(condition_factor, type_factor),
                                    design = ~ condition_factor)


        # using the function on our data
        DDS = DESeq(dds)

        norm_DDS = counts(DDS, normalized = TRUE) # normalization with respect to the sequencing depth

        # adding _norm onto the column names in the normalized matrix
        colnames(norm_DDS) = paste0(colnames(norm_DDS), "_norm")


        # provides independent filtering using the mean of normalized counts
        res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")


        # this is looking at the differences between the 3 deseq analyzed options
        countMatDiff = cbind(filt_countmatrix, norm_DDS, res)

        head(countMatDiff)




        # getting the results name and addding to the coef we want to shrink
        experiment_design_name = resultsNames(DDS)[2]

        # useful for visualization and ranking of genes or in this case peaks
        resLFC = lfcShrink(DDS, coef= resultsNames(DDS)[2], type = "apeglm")


        # finding the up and down regulated counts that pass the threshold
        

        

        up_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange >= 1) ,]
        
        
        down_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange <= -1) ,]
        

        #unchanging_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange <1 & resLFC\$log2FoldChange > -1) ,]
        unchanging_reg = resLFC[which(resLFC\$log2FoldChange <1 & resLFC\$log2FoldChange > -1) ,]

        #others_reg = resLFC[which(resLFC\$padj > 0.05), ] 
        others_reg = resLFC[which(resLFC\$padj > 0.05 & resLFC\$log2FoldChange <1 & resLFC\$log2FoldChange > -1), ]        

        # testing the chat gpt code to make the plot look publication ready

        # Relabel categories (for peaks, not genes)
        resLFC\$Label <- ifelse(resLFC\$padj > 0.05, "Not significant",
                            ifelse(resLFC\$log2FoldChange >= 1 & resLFC\$padj < 0.05, "Upregulated peaks",
                            ifelse(resLFC\$log2FoldChange <= -1 & resLFC\$padj < 0.05, "Downregulated peaks",
                            ifelse(resLFC\$log2FoldChange >= 1 & resLFC\$padj > 0.05, "Not significant",
                            ifelse(resLFC\$log2FoldChange <= -1 & resLFC\$padj > 0.05, "Not significant",
                                    "padj < 0.05 only")))))

        # Counts & percentages
        total <- nrow(resLFC)
        #up_count <- sum(resLFCLabel == "Upregulated peaks")
        #down_count <- sum(resLFCLabel == "Downregulated peaks")
        #ns_count <- sum(resLFCLabel == "Not significant")
        #padj_only_count <- sum(resLFCLabel == "padj < 0.05 only")

        up_pct <- round(100 * nrow(up_reg) / total, 1)
        down_pct <- round(100 * nrow(down_reg) / total, 1)
        ns_pct <- round(100 * nrow(others_reg) / total, 1)
        padj_only_pct <- round(100 * nrow(unchanging_reg) / total, 1)

        # Annotation text
        annotation_text <- paste0(
        "Upregulated peaks: ", nrow(up_reg), " (", up_pct, "%)\n",
        "Downregulated peaks: ", nrow(down_reg), " (", down_pct, "%)\n",
        "padj < 0.05 only: ", nrow(unchanging_reg), " (", padj_only_pct, "%)\n",
        "Not significant: ", nrow(others_reg), " (", ns_pct, "%)"
        )

        # Plot
        ma_plot_labeled <- ggplot(resLFC, aes(x = baseMean, y = log2FoldChange, color = Label)) +
        geom_point(alpha = 0.7, size = 0.5) +
        geom_hline(yintercept = 0, linetype = "solid", color = "black") +
        geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black") +
        scale_x_continuous(trans = "log10", 
                            breaks = scales::trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        scale_colour_manual(values = c("Upregulated peaks" = "red",
                                        "Downregulated peaks" = "red",
                                        "padj < 0.05 only" = "lightpink",
                                        "Not significant" = "grey70"),
                            breaks = c("Upregulated peaks", "Downregulated peaks", "padj < 0.05 only", "Not significant")) +
        labs(title = experiment_design_name,
            subtitle = paste("IDR threshold for masterPeaks:", "${params.return_idr}"),
            x = "Mean peak signal (log10 scale)",
            y = expression(Log[2]~Fold~Change),
            color = "Peak status") +
        theme_classic(base_size = 4) +
        theme(legend.position = "top",
                legend.title = element_text(size = 4),
                legend.text = element_text(size = 5),
                plot.title = element_text(size = 6, face = "bold"),
                plot.subtitle = element_text(size = 6),
                axis.title = element_text(size = 6),
                axis.text = element_text(size = 6)) +
        annotate("text", x = min(resLFC\$baseMean, na.rm = TRUE), 
                y = max(resLFC\$log2FoldChange, na.rm = TRUE), 
                hjust = -3, vjust = 1, 
                label = annotation_text, 
                size = 1.5)



        print(ma_plot_labeled)


        pdf(file = paste(experiment_design_name,"IDR", idr_used,"MA_plot_our_data_counts.pdf", sep = "_"), width = 4, height = 4)

        print(ma_plot_labeled)
        dev.off()

        png(filename = paste(experiment_design_name,"IDR", idr_used,"MA_plot_our_data_counts.png", sep = "_"), width = 4.3, height = 4, antialias = "subpixel", units = "in", res = 2000)

        print(ma_plot_labeled)
        dev.off()


        volcano_plot_removed_reps = EnhancedVolcano(resLFC,
                        lab = rownames(resLFC),
                        title = experiment_design_name, subtitle = paste("This means the first condition has these points more than in the second condition. The idr threshold for masterPeaks is",idr_used, sep = " "),
                        x = 'log2FoldChange', FCcutoff = 1,
                        y = 'padj', pCutoff = 0.05, labSize = 0 # making this 0 so it doesn't show the region names in the volcano plot 
                        
                        )

        print(volcano_plot_removed_reps)

        png(filename = paste(experiment_design_name,"IDR",idr_used,"volcano_plot_our_data_counts.png", sep = "_"), width = 1200, height = 900, antialias = "subpixel")
        print(volcano_plot_removed_reps)
        dev.off()

        pdf(file = paste(experiment_design_name,"IDR",idr_used,"volcano_plot_our_data_counts.pdf", sep = "_"), width = 16, height = 12)
        print(volcano_plot_removed_reps)
        dev.off()
        



        # using rlog over vst for transformation

        rld = rlog(DDS, blind=FALSE)

        #it is clear that the rlog transformation inherently accounts for differences in sequencing depth *
        head(assay(rld), 5)


        library(ggplot2)

        pcaData = plotPCA(rld, intgroup=c("condition_factor","type_factor"), returnData = TRUE )
        pcaData

        percentVar <- round(100 * attr(pcaData, "percentVar"))


        pca_plot_rlog = ggplot(pcaData, aes(PC1, PC2, color=condition_factor, shape=type_factor)) +
        ggtitle(paste("PCA plot using Rlog transform IDR", idr_used, sep = " ")) +
        theme(plot.caption = element_text(size = 30, face = "bold"), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15))+
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed()
        #pca_plot
        name_for_rlog_png_file =paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used, "rlog.png", sep = "_")

        png(filename = name_for_rlog_png_file, width = 900, height = 900, antialias = "subpixel")
        print(pca_plot_rlog)
        dev.off()

        name_for_rlog_pdf_file =paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used, "rlog.pdf", sep = "_")

        pdf(file = name_for_rlog_pdf_file, width = 10, height = 10)
        print(pca_plot_rlog)
        dev.off()


        # testing with vst
        #vsd_t = vst(DDS, blind = FALSE)

        # for histone mark k36me2 vst fails so i have to use the direct function
        vsd_t <- tryCatch({
            vst(DDS, blind = FALSE)
        }, error = function(e) {
            message("vst() failed, using varianceStabilizingTransformation() instead.")
            varianceStabilizingTransformation(DDS, blind = FALSE)
        })

        head(assay(vsd_t), 5)

        pcaData2 = plotPCA(vsd_t, intgroup=c("condition_factor","type_factor"), returnData = TRUE )
        pcaData2

        percentVar <- round(100 * attr(pcaData2, "percentVar"))
        pca_plot_vst = ggplot(pcaData2, aes(PC1, PC2, color=condition_factor, shape=type_factor)) +
        ggtitle(paste("PCA plot using VST transform IDR", idr_used, sep = " ")) +
        theme(plot.caption = element_text(size = 30, face = "bold"), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15))+
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed()


        #name_for_vst_png_file =paste(experiment_design_name,"PCA_plot_IDR", idr_used, "vst.png", sep = "_")

        pdf(file = paste(experiment_design_name,"PCA_plot_IDR_our_data_counts",idr_used,"vst_all_reps.pdf", sep = "_"), width = 10, height = 10)
        print(pca_plot_vst)
        dev.off()

        png(file = paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used,"vst_all_reps.png", sep = "_"), width = 900, height = 900)
        print(pca_plot_vst)
        dev.off()

        pca_plot_rlog

        pca_plot_vst


        rtracklayer::export.bed(row.names(up_reg), con = "${up_peaks_out}")

        rtracklayer::export.bed(row.names(down_reg), con = "${down_peaks_out}")

        # now exporting the unchanging peaks

        rtracklayer::export.bed(row.names(unchanging_reg), con = "${unchanging_peaks_out}")

        # export the other peaks that are greater than 0.05 and less than |1|
        rtracklayer::export.bed(row.names(others_reg), con = "${other_peaks_out}")

        # now hoping to get the peak lengths histone


        ###### if any errors happen here then dont do anything ######
        tryCatch({
        up_peaks = read.table(file = './${up_peaks_out}', header = FALSE, sep = "\t")

        up_peak_lengths = up_peaks\$V3 - up_peaks\$V2
        print(max(up_peak_lengths))

        png(filename = "hist_up_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
        hist.default(up_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        pdf(file = "hist_up_regulated_${idr_histone}_peak_lengths.pdf", height = 8, width = 8)
        hist.default(up_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        # just to view it here, not needed in nextflow here.
        #hist.default(up_peak_lengths, xaxt = "n", breaks = 1e2)
        #axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))

        down_peaks = read.table(file = './${down_peaks_out}', header = FALSE, sep = "\t")

        down_peak_lengths = down_peaks\$V3 - down_peaks\$V2
        print(max(down_peak_lengths))

        png(filename = "hist_down_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
        hist.default(down_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        pdf(file = "hist_down_regulated_${idr_histone}_peak_lengths.pdf", height = 8, width = 8)
        hist.default(down_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()


        unchanging_peaks = read.table(file = './${unchanging_peaks_out}', header = FALSE, sep = "\t")

        unchanging_peak_lengths = unchanging_peaks\$V3 - unchanging_peaks\$V2
        print(max(unchanging_peak_lengths))

        png(filename = "hist_unchanging_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
        hist.default(unchanging_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        pdf(file = "hist_unchanging_regulated_${idr_histone}_peak_lengths.pdf", height = 8, width = 8)
        hist.default(unchanging_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        ########## now violin plots for peak lengths ################

        up_peak_df = DataFrame(peak_lengths_${idr_histone} = up_peak_lengths, category = "up_peaks_lengths")
        down_peak_df = DataFrame(peak_lengths_${idr_histone} = down_peak_lengths, category = "down_peaks_lengths")
        unchanging_peak_df = DataFrame(peak_lengths_${idr_histone} = unchanging_peak_lengths, category = "unchanging_peaks_lengths")


        df_all_peak_lengths = rbind(up_peak_df, down_peak_df, unchanging_peak_df)

        #df_all_peak_lengths

        df_all_peak_lengths_gg =  ggplot(df_all_peak_lengths, aes( category, peak_lengths_${idr_histone}))
        
        all_peak_lengths_violin = df_all_peak_lengths_gg+
            geom_violin()+
            scale_y_continuous(labels = scales::label_number())
        ggsave("${idr_histone}_up_down_unchanging_peak_length_violin_plot.pdf", plot = all_peak_lengths_violin)
        
        
        }, error = function(x) {
        
        message("some of the peak files had no lenght so plotting is pointless")
        })


        tryCatch({
        # now for the annotated peaks, if an error occurs, don't do anything
        # lets get the annotated peaks

        library(ChIPseeker)

        library(TxDb.Hsapiens.UCSC.hg38.knownGene)

        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

        # now to read in the peak files
        peak_files = c(up_peaks = './${up_peaks_out}', down_peaks = './${down_peaks_out}', unchanging_peaks = './${unchanging_peaks_out}')
        up_peak_file = readPeakFile(peakfile = './${up_peaks_out}')


        unchanging_peak_file = readPeakFile(peakfile = './${unchanging_peaks_out}')
        unchanging_annotated_peaks = annotatePeak(peak= unchanging_peak_file, tssRegion = c(-5000, 5000), TxDb = txdb)

        up_annotated_peaks = annotatePeak(peak= up_peak_file, tssRegion = c(-5000, 5000), TxDb = txdb)
        # I would have to plot this pie chart as many times as there are peak files, but not yet.
        #plotAnnoPie(up_annotated_peaks)

        down_peak_file = readPeakFile(peakfile = './${down_peaks_out}')
        down_annotated_peaks = annotatePeak(peak= down_peak_file, tssRegion = c(-5000, 5000), TxDb = txdb)

        # not plotting this right now
        #plotAnnoPie(down_annotated_peaks)


        # not plotting this right now also
        #plotDistToTSS(down_annotated_peaks,
        #            title="Distribution of H3K27me3 peaks relative to TSS")


        annotated_peaks_list = lapply(peak_files, annotatePeak, tssRegion = c(-5000, 5000), TxDb = txdb) 


        # now get the number of total annotated peaks
        down_annotated_df <- as.data.frame(down_annotated_peaks)
        #num_down_annotated_peaks = count(down_annotated_df)[[1]]
        num_down_annotated_peaks = length(down_annotated_df[,1])

        up_annotated_df <- as.data.frame(up_annotated_peaks)
        #num_up_annotated_peaks = count(up_annotated_df)[[1]]
        num_up_annotated_peaks = length(up_annotated_df[,1])

        unchanging_annotated_df <- as.data.frame(unchanging_annotated_peaks)
        #num_unchanging_annotated_peaks = count(unchanging_annotated_df)[[1]]
        num_unchanging_annotated_peaks = length(unchanging_annotated_df[,1])

        # then plot with bar because it uses a ggplot object
        plotAnnoBar(annotated_peaks_list) + ggtitle("Percentage of ${idr_histone} peaks in Genomic Regions", subtitle = paste("Number of down ${idr_histone} peaks annotated ", num_down_annotated_peaks, "\nNumber of up ${idr_histone} peaks annotated ",num_up_annotated_peaks, "\nNumber of unchanging ${idr_histone} peaks annotated ", num_unchanging_annotated_peaks ) )

        ggsave("annotation_bar_graph_${idr_condition}_${idr_histone}_peaks.pdf", plot = last_plot())

        }, error = function(x) {
        
        message("for making annotated peaks, some of the files might have no differential peaks")
        })

        """


    }
    else {

    
        """
        #!/usr/bin/env Rscript

        bam_list = c(${bam_name_list})

        print(bam_list)

        #peaklist = list("\${first_idr_peak}", "\${second_idr_peak}")
        #peaklist = \${peaks_list}

        #print(peaklist)

        #best_hlow_idr
        #best_scrm_idr

        library(DESeq2)
        library(GenomicRanges)
        library(chromVAR)
        library(tidyr)
        library(EnhancedVolcano)
        library(readr)

        print(dir())
        # making the master peak genomic ranges object
        
        #mPeak = GRanges()

        #for (peakfile in c("./\${first_idr_peak}", "./\${second_idr_peak}")) {
        
            #peaktable = read.table(peakfile, header = FALSE, sep = "\t")
        
            #gr_object = GRanges(seqnames = peaktable\$V1, IRanges(start = peaktable\$V2, end = peaktable\$V3), strand = "*" )
        
            #gr_object = rtracklayer::import(peakfile, format = "BED")

            #mPeak = append(mPeak, gr_object)
        #}

        #peaktable = read.table("\${idr_peak_name}", header = FALSE)
        #mPeak = GRanges(seqnames = peaktable\$V1, IRanges(start = peaktable\$V2, end = peaktable\$V3), strand = "*" )

        mPeak = rtracklayer::import("${idr_peak_name}", format = "BED")

        # making sure there are no redundant peaks
        # this will allow me to resizd all of the regions by 50kb on each side, so any regions that are less than 100kb away from eachother will be merged by reduce because they overlap
        
        
        #masterPeak2 = reduce(mPeak)

        # using resize to extnd regions
        extended_master_peak = resize(mPeak, width = width(mPeak)+${resize_num}, fix = "center")

        # before trimming i need to add the sequence lengths of the correct genome
        # will have to find a way to automate this step for other genomes

        library(BSgenome.Hsapiens.UCSC.hg38)
        #seq_lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
        #seqlengths(extended_master_peak) <- seq_lengths[seqlevels(extended_master_peak)]
        #trimmed_master_peak = trim(extended_master_peak)
        #masterPeak_beta = reduce(trimmed_master_peak)

        masterPeak_beta = reduce(extended_master_peak)
        # now to resize again to remove the artificial 50kb from the start and end of the newly reduced peak
        
        new_peak_size = width(masterPeak_beta)-${resize_num}
        masterPeak = resize(masterPeak_beta,  width = new_peak_size , fix = "center")
        
        #print(masterPeak)

        # now I want to keep the standard chromosomes
        #seqnames_to_keep =masterPeak@seqnames@values[1:23]

        #masterPeak = masterPeak[seqnames(masterPeak) %in% seqnames_to_keep]

        # hoping to export my GRanges object master peaks to a bed file
        #write_tsv("\${master_peak_export_out}")
        rtracklayer::export.bed(masterPeak, con = "${master_peak_export_out}")

        

        list_of_mBams = c()
        
        for (bam in bam_list) {
    
            tokens = strsplit(bam, split = "_")[[1]]
            #print(tokens)
            
            histone = tokens[2]
            #print(histone)

            #list_of_mBams = list()

            if ("${idr_histone}" == histone) {
                #print(bam)
                list_of_mBams = c(list_of_mBams, bam)

            }
        }

        print(list_of_mBams)
    
        
        # now doing the next section to make the count matrix and the condition and experiment design

        # making the matrix 

        countsMatrix = matrix(data = NA, length(masterPeak), length(list_of_mBams) )


        # I want to get the chr start end of the peaks and put them in the matrix as column names so I know which peaks I am looking at from the deseq2 results in the differential peaks.
        seqnames(masterPeak)

        # i should put the unique ids from the master peak object in the matrix as row names.
        unique_masterPeak_ids = paste0(as.character(seqnames(masterPeak)), ":" ,start(masterPeak), "-", end(masterPeak))

        rownames(countsMatrix) = unique_masterPeak_ids

        ####################### new version ##########################
        library(Rsubread)
        library(dplyr)

        # featureCounts in rsubread needs a dataframe and in the correct format, with a GeneID column first then chr1 column
        df_masterPeak = as.data.frame(masterPeak)

        # now I need to change the second column name to chr
        df_masterPeak = df_masterPeak %>% rename(chr = seqnames)

        # next i need to add the unique ids as the columns in this dataframe
        df_masterPeak\$GeneID = unique_masterPeak_ids
        ##############################################################

        #getting the list of bam base names to add to the matrix column names
        list_bam_basenames = list()


        # for deseq2 i need the condition design. so hlow and scrm
        # then i will find a way to tally how many of each are there so i can automate the rep count
        condition_design = list()

        type_design = list()

        for (x in c(1:length(list_of_mBams))) {
        
        path_bam = list_of_mBams[[x]]
        print(path_bam)
        bam_basename = basename(as.character(path_bam))
        bam_tokens = strsplit(basename(as.character(path_bam)), split = "_")[[1]]
        
        # labeling the important tokens so it is easier to keep track of
        condition = bam_tokens[1]
        histone = bam_tokens[2]
        replicate = bam_tokens[3]
        
        
        # for later parts I need the condition and to know how many times to repeat them
        condition_design = append(condition_design, paste(condition,histone,sep="_"))
        
        
        # also get the replicate too for type design
        type_design = append(type_design, replicate)
        #type_design = append(type_design, paste(histone,replicate, sep="_"))
        
        
        
        # using chromVAR getCounts function
        # not using chromVAR anymore, for now
        #fragment_counts = getCounts(path_bam, masterPeak, paired= TRUE, by_rg = FALSE, format = "bam")

        fragment_counts <- featureCounts(
            files = path_bam,
            annot.ext = df_masterPeak,
            isPairedEnd = TRUE,      # TRUE if your data were paired-end
            strandSpecific = 0,       # 0 = unstranded, 1/2 if strand-specific
            useMetaFeatures = TRUE
        )
        
        # putting the fragment counts in the column labeled by the bam name
        countsMatrix[,x] = fragment_counts\$counts
        
        list_bam_basenames = append(list_bam_basenames, bam_basename)
        
        }

        colnames(countsMatrix) = list_bam_basenames

        ## this part i would have to tell the use to name their conditions as control and treatment so it can make the design treatment vs control

        #first removing the low count fragments. only keeping the rows that are above 5 count in total

        keep_Rows = which(rowSums(countsMatrix) > 5)

        filt_countmatrix = countsMatrix[keep_Rows,]

        # now to get the condition_design
        condition_counts = table(unlist(condition_design))

        # this gives back the names and the counts give back the counts for each of the names
        condition_names = names(condition_counts)
        condition_num = as.numeric(condition_counts)




        # now i can put both lists in to get back the experiment design
        #condition_factor = factor(rep(condition_names, times=condition_num))

        # Build condition factor in column order
        condition_factor <- factor(unlist(condition_design))

        # I want it so the treatment is second in the list here so it is first in the experimental design in deseq2. This way the experimental design will have the diff results going up or down in the treatment vs control

        # for treatment I need to make nextflow take the users input so they can relevel this section
        #treatment = "H1low_H3k27me3"
        treatment = "${params.treatment_name}_${idr_histone}"

        # this doesnt work because you need to put ref = levels(condition_factor)[2] to assign the actual relevel
        #if (levels(condition_factor)[1] == treatment ) {
        #condition_factor = relevel(condition_factor, levels(condition_factor)[2])
        #}else {
        #condition_factor
        #}

        # have to do this instead to make the treatment always the baseline
        #condition_factor <- relevel(condition_factor, ref = treatment)

        levels_now <- levels(condition_factor)
        control <- levels_now[ levels_now != treatment ]  # the other level

        # Set control as reference → LFC = treatment / control
        condition_factor <- relevel(condition_factor, ref = control)


        print(condition_factor)

        # repeating the above to have another column with type (replicates)
        #type_counts = table(unlist(type_design))
        #type_names = names(type_counts)
        #type_num = as.numeric(type_counts)

        #type_factor = factor(rep(type_names, times=type_counts))
        #type_factor = factor(rep(type_names, times=type_num[1]))

        # Build type factor in column order
        type_factor <- factor(unlist(type_design))

        # I want to get the idr threshold and use that as input for the file names and other things

        #peak_file = basename(peaklist[[1]])
        peak_file = basename("./${idr_peak_name}")

        idr_used = strsplit(peak_file, split = "_")[[1]][10]
        idr_used

        # now for deseq2 workflow


        # now to do the normal deseq2 workflow

        dds = DESeqDataSetFromMatrix(countData = filt_countmatrix,
                                    colData = DataFrame(condition_factor, type_factor),
                                    design = ~ condition_factor)


        # using the function on our data
        DDS = DESeq(dds)

        norm_DDS = counts(DDS, normalized = TRUE) # normalization with respect to the sequencing depth

        # adding _norm onto the column names in the normalized matrix
        colnames(norm_DDS) = paste0(colnames(norm_DDS), "_norm")


        # provides independent filtering using the mean of normalized counts
        res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")


        # this is looking at the differences between the 3 deseq analyzed options
        countMatDiff = cbind(filt_countmatrix, norm_DDS, res)

        head(countMatDiff)




        # getting the results name and addding to the coef we want to shrink
        experiment_design_name = resultsNames(DDS)[2]

        # useful for visualization and ranking of genes or in this case peaks
        resLFC = lfcShrink(DDS, coef= resultsNames(DDS)[2], type = "apeglm")


        # finding the up and down regulated counts that pass the threshold

        

        up_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange >= 1) ,]
        
        
        down_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange <= -1) ,]
        

        #unchanging_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange <1 & resLFC\$log2FoldChange > -1) ,]
        unchanging_reg = resLFC[which(resLFC\$log2FoldChange <1 & resLFC\$log2FoldChange > -1) ,]

        #others_reg = resLFC[which(resLFC\$padj > 0.05), ] 
        others_reg = resLFC[which(resLFC\$padj > 0.05 & resLFC\$log2FoldChange <1 & resLFC\$log2FoldChange > -1), ]        

        # testing the chat gpt code to make the plot look publication ready

        # Relabel categories (for peaks, not genes)
        resLFC\$Label <- ifelse(resLFC\$padj > 0.05, "Not significant",
                            ifelse(resLFC\$log2FoldChange >= 1 & resLFC\$padj < 0.05, "Upregulated peaks",
                            ifelse(resLFC\$log2FoldChange <= -1 & resLFC\$padj < 0.05, "Downregulated peaks",
                            ifelse(resLFC\$log2FoldChange >= 1 & resLFC\$padj > 0.05, "Not significant",
                            ifelse(resLFC\$log2FoldChange <= -1 & resLFC\$padj > 0.05, "Not significant",
                                    "padj < 0.05 only")))))

        # Counts & percentages
        total <- nrow(resLFC)
        #up_count <- sum(resLFCLabel == "Upregulated peaks")
        #down_count <- sum(resLFCLabel == "Downregulated peaks")
        #ns_count <- sum(resLFCLabel == "Not significant")
        #padj_only_count <- sum(resLFCLabel == "padj < 0.05 only")

        up_pct <- round(100 * nrow(up_reg) / total, 1)
        down_pct <- round(100 * nrow(down_reg) / total, 1)
        ns_pct <- round(100 * nrow(others_reg) / total, 1)
        padj_only_pct <- round(100 * nrow(unchanging_reg) / total, 1)

        # Annotation text
        annotation_text <- paste0(
        "Upregulated peaks: ", nrow(up_reg), " (", up_pct, "%)\n",
        "Downregulated peaks: ", nrow(down_reg), " (", down_pct, "%)\n",
        "padj < 0.05 only: ", nrow(unchanging_reg), " (", padj_only_pct, "%)\n",
        "Not significant: ", nrow(others_reg), " (", ns_pct, "%)"
        )

        # Plot
        ma_plot_labeled <- ggplot(resLFC, aes(x = baseMean, y = log2FoldChange, color = Label)) +
        geom_point(alpha = 0.7, size = 0.5) +
        geom_hline(yintercept = 0, linetype = "solid", color = "black") +
        geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black") +
        scale_x_continuous(trans = "log10", 
                            breaks = scales::trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        scale_colour_manual(values = c("Upregulated peaks" = "red",
                                        "Downregulated peaks" = "red",
                                        "padj < 0.05 only" = "lightpink",
                                        "Not significant" = "grey70"),
                            breaks = c("Upregulated peaks", "Downregulated peaks", "padj < 0.05 only", "Not significant")) +
        labs(title = experiment_design_name,
            subtitle = paste("IDR threshold for masterPeaks:", "${params.return_idr}"),
            x = "Mean peak signal (log10 scale)",
            y = expression(Log[2]~Fold~Change),
            color = "Peak status") +
        theme_classic(base_size = 4) +
        theme(legend.position = "top",
                legend.title = element_text(size = 4),
                legend.text = element_text(size = 5),
                plot.title = element_text(size = 6, face = "bold"),
                plot.subtitle = element_text(size = 6),
                axis.title = element_text(size = 6),
                axis.text = element_text(size = 6)) +
        annotate("text", x = min(resLFC\$baseMean, na.rm = TRUE), 
                y = max(resLFC\$log2FoldChange, na.rm = TRUE), 
                hjust = -3, vjust = 1, 
                label = annotation_text, 
                size = 1.5)



        print(ma_plot_labeled)


        pdf(file = paste(experiment_design_name,"mergedBy", "${merge_dist}","MA_plot_our_data_counts.pdf", sep = "_"), width = 4, height = 4)

        print(ma_plot_labeled)
        dev.off()

        png(filename = paste(experiment_design_name,"mergedBy", "${merge_dist}","MA_plot_our_data_counts.png", sep = "_"), width = 4.3, height = 4, antialias = "subpixel", units = "in", res = 2000)

        print(ma_plot_labeled)
        dev.off()


        volcano_plot_removed_reps = EnhancedVolcano(resLFC,
                        lab = rownames(resLFC),
                        title = experiment_design_name, subtitle = paste("This means the first condition has these points more than in the second condition. The idr threshold for masterPeaks is",idr_used, sep = " "),
                        x = 'log2FoldChange', FCcutoff = 1,
                        y = 'padj', pCutoff = 0.05, labSize = 0 # making this 0 so it doesn't show the region names in the volcano plot 
                        
                        )

        print(volcano_plot_removed_reps)

        png(filename = paste(experiment_design_name,"mergedBy","${merge_dist}","volcano_plot_our_data_counts.png", sep = "_"), width = 1200, height = 900, antialias = "subpixel")
        print(volcano_plot_removed_reps)
        dev.off()

        pdf(file = paste(experiment_design_name,"mergedBy","${merge_dist}","volcano_plot_our_data_counts.pdf", sep = "_"), width = 16, height = 12)
        print(volcano_plot_removed_reps)
        dev.off()
        



        # using rlog over vst for transformation

        rld = rlog(DDS, blind=FALSE)

        #it is clear that the rlog transformation inherently accounts for differences in sequencing depth *
        head(assay(rld), 5)


        library(ggplot2)

        pcaData = plotPCA(rld, intgroup=c("condition_factor","type_factor"), returnData = TRUE )
        pcaData

        percentVar <- round(100 * attr(pcaData, "percentVar"))


        pca_plot_rlog = ggplot(pcaData, aes(PC1, PC2, color=condition_factor, shape=type_factor)) +
        ggtitle(paste("PCA plot using Rlog transform IDR", idr_used, sep = " ")) +
        theme(plot.caption = element_text(size = 30, face = "bold"), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15))+
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed()
        #pca_plot
        name_for_rlog_png_file =paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used, "rlog.png", sep = "_")

        png(filename = name_for_rlog_png_file, width = 900, height = 900, antialias = "subpixel")
        print(pca_plot_rlog)
        dev.off()

        name_for_rlog_pdf_file =paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used, "rlog.pdf", sep = "_")

        pdf(file = name_for_rlog_pdf_file, width = 10, height = 10)
        print(pca_plot_rlog)
        dev.off()


        # testing with vst
        #vsd_t = vst(DDS, blind = FALSE)

        # for histone mark k36me2 vst fails so i have to use the direct function
        vsd_t <- tryCatch({
            vst(DDS, blind = FALSE)
        }, error = function(e) {
            message("vst() failed, using varianceStabilizingTransformation() instead.")
            varianceStabilizingTransformation(DDS, blind = FALSE)
        })

        head(assay(vsd_t), 5)

        pcaData2 = plotPCA(vsd_t, intgroup=c("condition_factor","type_factor"), returnData = TRUE )
        pcaData2

        percentVar <- round(100 * attr(pcaData2, "percentVar"))
        pca_plot_vst = ggplot(pcaData2, aes(PC1, PC2, color=condition_factor, shape=type_factor)) +
        ggtitle(paste("PCA plot using VST transform IDR", idr_used, sep = " ")) +
        theme(plot.caption = element_text(size = 30, face = "bold"), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15))+
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed()


        #name_for_vst_png_file =paste(experiment_design_name,"PCA_plot_IDR", idr_used, "vst.png", sep = "_")

        pdf(file = paste(experiment_design_name,"PCA_plot_IDR_our_data_counts",idr_used,"vst_all_reps.pdf", sep = "_"), width = 10, height = 10)
        print(pca_plot_vst)
        dev.off()

        png(file = paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used,"vst_all_reps.png", sep = "_"), width = 900, height = 900)
        print(pca_plot_vst)
        dev.off()

        pca_plot_rlog

        pca_plot_vst


        rtracklayer::export.bed(row.names(up_reg), con = "${up_peaks_out}")

        rtracklayer::export.bed(row.names(down_reg), con = "${down_peaks_out}")

        # now exporting the unchanging peaks

        rtracklayer::export.bed(row.names(unchanging_reg), con = "${unchanging_peaks_out}")

        # export the other peaks that are greater than 0.05 and less than |1|
        rtracklayer::export.bed(row.names(others_reg), con = "${other_peaks_out}")

        # now hoping to get the peak lengths histone


        ###### if any errors happen here then dont do anything ######
        tryCatch({
        up_peaks = read.table(file = './${up_peaks_out}', header = FALSE, sep = "\t")

        up_peak_lengths = up_peaks\$V3 - up_peaks\$V2
        print(max(up_peak_lengths))

        png(filename = "hist_up_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
        hist.default(up_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        pdf(file = "hist_up_regulated_${idr_histone}_peak_lengths.pdf", height = 8, width = 8)
        hist.default(up_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        # just to view it here, not needed in nextflow here.
        #hist.default(up_peak_lengths, xaxt = "n", breaks = 1e2)
        #axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))

        down_peaks = read.table(file = './${down_peaks_out}', header = FALSE, sep = "\t")

        down_peak_lengths = down_peaks\$V3 - down_peaks\$V2
        print(max(down_peak_lengths))

        png(filename = "hist_down_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
        hist.default(down_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        pdf(file = "hist_down_regulated_${idr_histone}_peak_lengths.pdf", height = 8, width = 8)
        hist.default(down_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()


        unchanging_peaks = read.table(file = './${unchanging_peaks_out}', header = FALSE, sep = "\t")

        unchanging_peak_lengths = unchanging_peaks\$V3 - unchanging_peaks\$V2
        print(max(unchanging_peak_lengths))

        png(filename = "hist_unchanging_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
        hist.default(unchanging_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        pdf(file = "hist_unchanging_regulated_${idr_histone}_peak_lengths.pdf", height = 8, width = 8)
        hist.default(unchanging_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        ########## now violin plots for peak lengths ################

        up_peak_df = DataFrame(peak_lengths_${idr_histone} = up_peak_lengths, category = "up_peaks_lengths")
        down_peak_df = DataFrame(peak_lengths_${idr_histone} = down_peak_lengths, category = "down_peaks_lengths")
        unchanging_peak_df = DataFrame(peak_lengths_${idr_histone} = unchanging_peak_lengths, category = "unchanging_peaks_lengths")


        df_all_peak_lengths = rbind(up_peak_df, down_peak_df, unchanging_peak_df)

        #df_all_peak_lengths

        df_all_peak_lengths_gg =  ggplot(df_all_peak_lengths, aes( category, peak_lengths_${idr_histone}))
        
        all_peak_lengths_violin = df_all_peak_lengths_gg+
            geom_violin()+
            scale_y_continuous(labels = scales::label_number())
        ggsave("${idr_histone}_up_down_unchanging_peak_length_violin_plot.pdf", plot = all_peak_lengths_violin)
        
        
        }, error = function(x) {
        
        message("some of the peak files had no lenght so plotting is pointless")
        })


        tryCatch({
        # now for the annotated peaks, if an error occurs, don't do anything

        # lets get the annotated peaks

        library(ChIPseeker)

        library(TxDb.Hsapiens.UCSC.hg38.knownGene)

        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

        # now to read in the peak files
        peak_files = c(up_peaks = './${up_peaks_out}', down_peaks = './${down_peaks_out}', unchanging_peaks = './${unchanging_peaks_out}')
        up_peak_file = readPeakFile(peakfile = './${up_peaks_out}')


        unchanging_peak_file = readPeakFile(peakfile = './${unchanging_peaks_out}')
        unchanging_annotated_peaks = annotatePeak(peak= unchanging_peak_file, tssRegion = c(-5000, 5000), TxDb = txdb)

        up_annotated_peaks = annotatePeak(peak= up_peak_file, tssRegion = c(-5000, 5000), TxDb = txdb)
        # I would have to plot this pie chart as many times as there are peak files, but not yet.
        #plotAnnoPie(up_annotated_peaks)

        down_peak_file = readPeakFile(peakfile = './${down_peaks_out}')
        down_annotated_peaks = annotatePeak(peak= down_peak_file, tssRegion = c(-5000, 5000), TxDb = txdb)

        # not plotting this right now
        #plotAnnoPie(down_annotated_peaks)


        # not plotting this right now also
        #plotDistToTSS(down_annotated_peaks,
        #            title="Distribution of H3K27me3 peaks relative to TSS")


        annotated_peaks_list = lapply(peak_files, annotatePeak, tssRegion = c(-5000, 5000), TxDb = txdb) 


        # now get the number of total annotated peaks
        down_annotated_df <- as.data.frame(down_annotated_peaks)
        #num_down_annotated_peaks = count(down_annotated_df)[[1]]
        num_down_annotated_peaks = length(down_annotated_df[,1])

        up_annotated_df <- as.data.frame(up_annotated_peaks)
        #num_up_annotated_peaks = count(up_annotated_df)[[1]]
        num_up_annotated_peaks = length(up_annotated_df[,1])

        unchanging_annotated_df <- as.data.frame(unchanging_annotated_peaks)
        #num_unchanging_annotated_peaks = count(unchanging_annotated_df)[[1]]
        num_unchanging_annotated_peaks = length(unchanging_annotated_df[,1])

        # then plot with bar because it uses a ggplot object
        plotAnnoBar(annotated_peaks_list) + ggtitle("Percentage of ${idr_histone} peaks in Genomic Regions", subtitle = paste("Number of down ${idr_histone} peaks annotated ", num_down_annotated_peaks, "\nNumber of up ${idr_histone} peaks annotated ",num_up_annotated_peaks, "\nNumber of unchanging ${idr_histone} peaks annotated ", num_unchanging_annotated_peaks ) )

        ggsave("annotation_bar_graph_${idr_condition}_${idr_histone}_peaks.pdf", plot = last_plot())

        }, error = function(x) {
        
        message("for making annotated peaks, some of the files might have no differential peaks")
        })

        """
    }
    
    
}

process find_diff_peaks_R_process_SE {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/r_language'


    publishDir "./nextflow_R_script_outputs/${idr_histone}/", mode: 'copy', pattern: '*'
    // if (params.narrowPeak_data) {
    //     publishDir "./nextflow_R_script_outputs/narrow_Peaks_called/${idr_histone}/", mode: 'copy', pattern: '*'
    // }
    // else {
    //     publishDir "./nextflow_R_script_outputs/broad_Peaks_called/${idr_histone}/", mode: 'copy', pattern: '*'
    // }
    
    label 'normal_big_resources'

    //debug true

    input:

    tuple val(idr_condition), val(idr_histone), val(merge_dist), val(idr_peak_name), path(idr_peak_path)

    // these are all experiment bams.
    // i need to use the histone label, and condition label to get the correct bams in this process
    path(bam_path_list)
    


    
    output:

    path("*.{png,pdf,bed}"), emit: all_r_plots

    path("${master_peak_export_out}"), emit: master_peak_emit

    tuple val("${full_condition}"), val("${idr_histone}"), path("${up_peaks_out}"), path("${down_peaks_out}"), path("${unchanging_peaks_out}"), emit: diff_peaks_ch 

    path("${up_peaks_out}"), emit: up_peaks_emit
    path("${down_peaks_out}"), emit: down_peaks_emit
    path("${unchanging_peaks_out}"), emit: unchanging_peaks_emit
    path("${other_peaks_out}"), emit: other_peaks_emit



    script:

    //full_condition = "${idr_condition[0]}vs${idr_condition[1]}"
    full_condition = "${idr_condition}"

    //first_idr_peak = idr_peak_name[0]
    //second_idr_peak = idr_peak_name[1]

    bam_name_list = bam_path_list.collect { "\"${it.getName()}\"" }.join(',')

    //peaks_list = [first_idr_peak, second_idr_peak]

    // need the output file name for masterpeak export

    master_peak_export_out = "masterpeak_${idr_histone}_${merge_dist}_maxgap_${idr_condition}.bed"

    up_peaks_out = "up_${idr_histone}_${full_condition}_regulated_peaks.bed"
    down_peaks_out = "down_${idr_histone}_${full_condition}_regulated_peaks.bed"
    unchanging_peaks_out = "unchanging_${idr_histone}_${full_condition}_regulated_peaks.bed"
    other_peaks_out = "others_${idr_histone}_${full_condition}_regulated_peaks.bed"

    resize_num = params.masterPeak100kb ? 100000 : (params.masterPeak10kb ? 10000 : (params.masterPeak30kb ? 30000 : 1))

    // if (params.masterPeak100kb) {
    //     resize_num = 100000
    // }
    // if (params.masterPeak10kb) {
    //     resize_num = 10000
    // }
    // if (params.masterPeak30kb) {
    //     resize_num = 30000
    // }
    // if (!params.masterPeak100kb & !params.masterPeak10kb & !masterPeak30kb) {
    //     resize_num = 1
    // }

    if (params.narrowPeak_data || params.sicer2) {

        """
        #!/usr/bin/env Rscript

        bam_list = c(${bam_name_list})

        print(bam_list)

        #peaklist = list("\${first_idr_peak}", "\${second_idr_peak}")
        #peaklist = \${peaks_list}

        #print(peaklist)

        #best_hlow_idr
        #best_scrm_idr

        library(DESeq2)
        library(GenomicRanges)
        library(chromVAR)
        library(tidyr)
        library(EnhancedVolcano)
        library(readr)

        print(dir())
        # making the master peak genomic ranges object
        
        #mPeak = GRanges()

        #for (peakfile in c("./\${first_idr_peak}", "./\${second_idr_peak}")) {
        
            #peaktable = read.table(peakfile, header = FALSE, sep = "\t")
        
            #gr_object = GRanges(seqnames = peaktable\$V1, IRanges(start = peaktable\$V2, end = peaktable\$V3), strand = "*" )
        
            #gr_object = rtracklayer::import(peakfile, format = "BED")

            #mPeak = append(mPeak, gr_object)
        #}

        #peaktable = read.table("\${idr_peak_name}", header = FALSE)
        #mPeak = GRanges(seqnames = peaktable\$V1, IRanges(start = peaktable\$V2, end = peaktable\$V3), strand = "*" )

        # might have to let it auto detect the file type when using narrow peaks that arent merged into a strict bed file
        mPeak = rtracklayer::import("${idr_peak_name}", format = "BED")

        # making sure there are no redundant peaks
        # this will allow me to resizd all of the regions by 50kb on each side, so any regions that are less than 100kb away from eachother will be merged by reduce because they overlap
        
        
        #masterPeak2 = reduce(mPeak)

        # using resize to extnd regions
        extended_master_peak = resize(mPeak, width = width(mPeak)+1, fix = "center")

        # before trimming i need to add the sequence lengths of the correct genome
        # will have to find a way to automate this step for other genomes

        library(BSgenome.Hsapiens.UCSC.hg38)
        #seq_lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
        #seqlengths(extended_master_peak) <- seq_lengths[seqlevels(extended_master_peak)]
        #trimmed_master_peak = trim(extended_master_peak)
        #masterPeak_beta = reduce(trimmed_master_peak)

        masterPeak_beta = reduce(extended_master_peak)
        # now to resize again to remove the artificial 50kb from the start and end of the newly reduced peak
        
        new_peak_size = width(masterPeak_beta)-1
        masterPeak = resize(masterPeak_beta,  width = new_peak_size , fix = "center")
        
        #print(masterPeak)

        # now I want to keep the standard chromosomes
        #seqnames_to_keep =masterPeak@seqnames@values[1:23]

        #masterPeak = masterPeak[seqnames(masterPeak) %in% seqnames_to_keep]

        # hoping to export my GRanges object master peaks to a bed file
        #write_tsv("\${master_peak_export_out}")
        rtracklayer::export.bed(masterPeak, con = "${master_peak_export_out}")

        

        list_of_mBams = c()
        
        for (bam in bam_list) {
    
            tokens = strsplit(bam, split = "_")[[1]]
            #print(tokens)
            
            histone = tokens[2]
            #print(histone)

            #list_of_mBams = list()

            if ("${idr_histone}" == histone) {
                #print(bam)
                list_of_mBams = c(list_of_mBams, bam)

            }
        }

        print(list_of_mBams)
    
        
        # now doing the next section to make the count matrix and the condition and experiment design

        # making the matrix 

        countsMatrix = matrix(data = NA, length(masterPeak), length(list_of_mBams) )


        # I want to get the chr start end of the peaks and put them in the matrix as column names so I know which peaks I am looking at from the deseq2 results in the differential peaks.
        seqnames(masterPeak)

        # i should put the unique ids from the master peak object in the matrix as row names.
        unique_masterPeak_ids = paste0(as.character(seqnames(masterPeak)), ":" ,start(masterPeak), "-", end(masterPeak))

        rownames(countsMatrix) = unique_masterPeak_ids

        ####################### new version ##########################
        library(Rsubread)
        library(dplyr)

        # featureCounts in rsubread needs a dataframe and in the correct format, with a GeneID column first then chr1 column
        df_masterPeak = as.data.frame(masterPeak)

        # now I need to change the second column name to chr
        df_masterPeak = df_masterPeak %>% rename(chr = seqnames)

        # next i need to add the unique ids as the columns in this dataframe
        df_masterPeak\$GeneID = unique_masterPeak_ids
        ##############################################################

        #getting the list of bam base names to add to the matrix column names
        list_bam_basenames = list()


        # for deseq2 i need the condition design. so hlow and scrm
        # then i will find a way to tally how many of each are there so i can automate the rep count
        condition_design = list()

        type_design = list()

        for (x in c(1:length(list_of_mBams))) {
        
        
        path_bam = list_of_mBams[[x]]
        print(path_bam)
        bam_basename = basename(as.character(path_bam))
        bam_tokens = strsplit(basename(as.character(path_bam)), split = "_")[[1]]
        
        # labeling the important tokens so it is easier to keep track of
        condition = bam_tokens[1]
        histone = bam_tokens[2]
        replicate = bam_tokens[3]
        
        
        # for later parts I need the condition and to know how many times to repeat them
        condition_design = append(condition_design, paste(condition,histone,sep="_"))
        
        
        # also get the replicate too for type design
        type_design = append(type_design, replicate)
        #type_design = append(type_design, paste(histone,replicate, sep="_"))
        
        
        
        # using chromVAR getCounts function
        # not using chromVAR anymore, for now
        #fragment_counts = getCounts(path_bam, masterPeak, paired= TRUE, by_rg = FALSE, format = "bam")

        fragment_counts <- featureCounts(
            files = path_bam,
            annot.ext = df_masterPeak,
            isPairedEnd = FALSE,      # TRUE if your data were paired-end
            strandSpecific = 0,       # 0 = unstranded, 1/2 if strand-specific
            useMetaFeatures = TRUE
        )
        
        
        # putting the fragment counts in the column labeled by the bam name
        countsMatrix[,x] = fragment_counts\$counts
        
        list_bam_basenames = append(list_bam_basenames, bam_basename)
        
        }

        colnames(countsMatrix) = list_bam_basenames

        ## this part i would have to tell the use to name their conditions as control and treatment so it can make the design treatment vs control

        #first removing the low count fragments. only keeping the rows that are above 5 count in total

        keep_Rows = which(rowSums(countsMatrix) > 5)

        filt_countmatrix = countsMatrix[keep_Rows,]

        # now to get the condition_design
        condition_counts = table(unlist(condition_design))

        # this gives back the names and the counts give back the counts for each of the names
        condition_names = names(condition_counts)
        condition_num = as.numeric(condition_counts)




        # now i can put both lists in to get back the experiment design
        #condition_factor = factor(rep(condition_names, times=condition_num))

        # Build condition factor in column order
        condition_factor <- factor(unlist(condition_design))

        # I want it so the treatment is second in the list here so it is first in the experimental design in deseq2. This way the experimental design will have the diff results going up or down in the treatment vs control

        # for treatment I need to make nextflow take the users input so they can relevel this section
        #treatment = "H1low_H3k27me3"
        treatment = "${params.treatment_name}_${idr_histone}"

        # this doesnt work because you need to put ref = levels(condition_factor)[2] to assign the actual relevel
        #if (levels(condition_factor)[1] == treatment ) {
        #condition_factor = relevel(condition_factor, levels(condition_factor)[2])
        #}else {
        #condition_factor
        #}

        # have to do this instead to make the treatment always the baseline
        #condition_factor <- relevel(condition_factor, ref = treatment)

        levels_now <- levels(condition_factor)
        control <- levels_now[ levels_now != treatment ]  # the other level

        # Set control as reference → LFC = treatment / control
        condition_factor <- relevel(condition_factor, ref = control)


        print(condition_factor)

        # repeating the above to have another column with type (replicates)
        #type_counts = table(unlist(type_design))
        #type_names = names(type_counts)
        #type_num = as.numeric(type_counts)

        #type_factor = factor(rep(type_names, times=type_counts))
        #type_factor = factor(rep(type_names, times=type_num[1]))

        # Build type factor in column order
        type_factor <- factor(unlist(type_design))

        # I want to get the idr threshold and use that as input for the file names and other things

        #peak_file = basename(peaklist[[1]])
        peak_file = basename("./${idr_peak_name}")

        idr_used = strsplit(peak_file, split = "_")[[1]][10]
        idr_used

        # now for deseq2 workflow


        # now to do the normal deseq2 workflow

        dds = DESeqDataSetFromMatrix(countData = filt_countmatrix,
                                    colData = DataFrame(condition_factor, type_factor),
                                    design = ~ condition_factor)


        # using the function on our data
        DDS = DESeq(dds)

        norm_DDS = counts(DDS, normalized = TRUE) # normalization with respect to the sequencing depth

        # adding _norm onto the column names in the normalized matrix
        colnames(norm_DDS) = paste0(colnames(norm_DDS), "_norm")


        # provides independent filtering using the mean of normalized counts
        res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")


        # this is looking at the differences between the 3 deseq analyzed options
        countMatDiff = cbind(filt_countmatrix, norm_DDS, res)

        head(countMatDiff)




        # getting the results name and addding to the coef we want to shrink
        experiment_design_name = resultsNames(DDS)[2]

        # useful for visualization and ranking of genes or in this case peaks
        resLFC = lfcShrink(DDS, coef= resultsNames(DDS)[2], type = "apeglm")


        # finding the up and down regulated counts that pass the threshold
        

        

        up_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange >= 1) ,]
        
        
        down_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange <= -1) ,]
        

        #unchanging_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange <1 & resLFC\$log2FoldChange > -1) ,]
        unchanging_reg = resLFC[which(resLFC\$log2FoldChange <1 & resLFC\$log2FoldChange > -1) ,]

        #others_reg = resLFC[which(resLFC\$padj > 0.05), ] 
        others_reg = resLFC[which(resLFC\$padj > 0.05 & resLFC\$log2FoldChange <1 & resLFC\$log2FoldChange > -1), ]         

        # testing the chat gpt code to make the plot look publication ready

        # Relabel categories (for peaks, not genes)
        resLFC\$Label <- ifelse(resLFC\$padj > 0.05, "Not significant",
                            ifelse(resLFC\$log2FoldChange >= 1 & resLFC\$padj < 0.05, "Upregulated peaks",
                            ifelse(resLFC\$log2FoldChange <= -1 & resLFC\$padj < 0.05, "Downregulated peaks",
                            ifelse(resLFC\$log2FoldChange >= 1 & resLFC\$padj > 0.05, "Not significant",
                            ifelse(resLFC\$log2FoldChange <= -1 & resLFC\$padj > 0.05, "Not significant",
                                    "padj < 0.05 only")))))

        # Counts & percentages
        total <- nrow(resLFC)
        #up_count <- sum(resLFCLabel == "Upregulated peaks")
        #down_count <- sum(resLFCLabel == "Downregulated peaks")
        #ns_count <- sum(resLFCLabel == "Not significant")
        #padj_only_count <- sum(resLFCLabel == "padj < 0.05 only")

        up_pct <- round(100 * nrow(up_reg) / total, 1)
        down_pct <- round(100 * nrow(down_reg) / total, 1)
        ns_pct <- round(100 * nrow(others_reg) / total, 1)
        padj_only_pct <- round(100 * nrow(unchanging_reg) / total, 1)

        # Annotation text
        annotation_text <- paste0(
        "Upregulated peaks: ", nrow(up_reg), " (", up_pct, "%)\n",
        "Downregulated peaks: ", nrow(down_reg), " (", down_pct, "%)\n",
        "padj < 0.05 only: ", nrow(unchanging_reg), " (", padj_only_pct, "%)\n",
        "Not significant: ", nrow(others_reg), " (", ns_pct, "%)"
        )

        # Plot
        ma_plot_labeled <- ggplot(resLFC, aes(x = baseMean, y = log2FoldChange, color = Label)) +
        geom_point(alpha = 0.7, size = 0.5) +
        geom_hline(yintercept = 0, linetype = "solid", color = "black") +
        geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black") +
        scale_x_continuous(trans = "log10", 
                            breaks = scales::trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        scale_colour_manual(values = c("Upregulated peaks" = "red",
                                        "Downregulated peaks" = "red",
                                        "padj < 0.05 only" = "lightpink",
                                        "Not significant" = "grey70"),
                            breaks = c("Upregulated peaks", "Downregulated peaks", "padj < 0.05 only", "Not significant")) +
        labs(title = experiment_design_name,
            subtitle = paste("IDR threshold for masterPeaks:", "${params.return_idr}"),
            x = "Mean peak signal (log10 scale)",
            y = expression(Log[2]~Fold~Change),
            color = "Peak status") +
        theme_classic(base_size = 4) +
        theme(legend.position = "top",
                legend.title = element_text(size = 4),
                legend.text = element_text(size = 5),
                plot.title = element_text(size = 6, face = "bold"),
                plot.subtitle = element_text(size = 6),
                axis.title = element_text(size = 6),
                axis.text = element_text(size = 6)) +
        annotate("text", x = min(resLFC\$baseMean, na.rm = TRUE), 
                y = max(resLFC\$log2FoldChange, na.rm = TRUE), 
                hjust = -3, vjust = 1, 
                label = annotation_text, 
                size = 1.5)



        print(ma_plot_labeled)


        pdf(file = paste(experiment_design_name,"IDR", idr_used,"MA_plot_our_data_counts.pdf", sep = "_"), width = 4, height = 4)

        print(ma_plot_labeled)
        dev.off()

        png(filename = paste(experiment_design_name,"IDR", idr_used,"MA_plot_our_data_counts.png", sep = "_"), width = 4.3, height = 4, antialias = "subpixel", units = "in", res = 2000)

        print(ma_plot_labeled)
        dev.off()


        volcano_plot_removed_reps = EnhancedVolcano(resLFC,
                        lab = rownames(resLFC),
                        title = experiment_design_name, subtitle = paste("This means the first condition has these points more than in the second condition. The idr threshold for masterPeaks is",idr_used, sep = " "),
                        x = 'log2FoldChange', FCcutoff = 1,
                        y = 'padj', pCutoff = 0.05, labSize = 0 # making this 0 so it doesn't show the region names in the volcano plot 
                        
                        )

        print(volcano_plot_removed_reps)

        png(filename = paste(experiment_design_name,"IDR",idr_used,"volcano_plot_our_data_counts.png", sep = "_"), width = 1200, height = 900, antialias = "subpixel")
        print(volcano_plot_removed_reps)
        dev.off()

        pdf(file = paste(experiment_design_name,"IDR",idr_used,"volcano_plot_our_data_counts.pdf", sep = "_"), width = 16, height = 12)
        print(volcano_plot_removed_reps)
        dev.off()
        



        # using rlog over vst for transformation

        rld = rlog(DDS, blind=FALSE)

        #it is clear that the rlog transformation inherently accounts for differences in sequencing depth *
        head(assay(rld), 5)


        library(ggplot2)

        pcaData = plotPCA(rld, intgroup=c("condition_factor","type_factor"), returnData = TRUE )
        pcaData

        percentVar <- round(100 * attr(pcaData, "percentVar"))


        pca_plot_rlog = ggplot(pcaData, aes(PC1, PC2, color=condition_factor, shape=type_factor)) +
        ggtitle(paste("PCA plot using Rlog transform IDR", idr_used, sep = " ")) +
        theme(plot.caption = element_text(size = 30, face = "bold"), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15))+
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed()
        #pca_plot
        name_for_rlog_png_file =paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used, "rlog.png", sep = "_")

        png(filename = name_for_rlog_png_file, width = 900, height = 900, antialias = "subpixel")
        print(pca_plot_rlog)
        dev.off()

        name_for_rlog_pdf_file =paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used, "rlog.pdf", sep = "_")

        pdf(file = name_for_rlog_pdf_file, width = 10, height = 10)
        print(pca_plot_rlog)
        dev.off()


        # testing with vst
        #vsd_t = vst(DDS, blind = FALSE)

        # for histone mark k36me2 vst fails so i have to use the direct function
        vsd_t <- tryCatch({
            vst(DDS, blind = FALSE)
        }, error = function(e) {
            message("vst() failed, using varianceStabilizingTransformation() instead.")
            varianceStabilizingTransformation(DDS, blind = FALSE)
        })

        head(assay(vsd_t), 5)

        pcaData2 = plotPCA(vsd_t, intgroup=c("condition_factor","type_factor"), returnData = TRUE )
        pcaData2

        percentVar <- round(100 * attr(pcaData2, "percentVar"))
        pca_plot_vst = ggplot(pcaData2, aes(PC1, PC2, color=condition_factor, shape=type_factor)) +
        ggtitle(paste("PCA plot using VST transform IDR", idr_used, sep = " ")) +
        theme(plot.caption = element_text(size = 30, face = "bold"), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15))+
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed()


        #name_for_vst_png_file =paste(experiment_design_name,"PCA_plot_IDR", idr_used, "vst.png", sep = "_")

        pdf(file = paste(experiment_design_name,"PCA_plot_IDR_our_data_counts",idr_used,"vst_all_reps.pdf", sep = "_"), width = 10, height = 10)
        print(pca_plot_vst)
        dev.off()

        png(file = paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used,"vst_all_reps.png", sep = "_"), width = 900, height = 900)
        print(pca_plot_vst)
        dev.off()

        pca_plot_rlog

        pca_plot_vst


        rtracklayer::export.bed(row.names(up_reg), con = "${up_peaks_out}")

        rtracklayer::export.bed(row.names(down_reg), con = "${down_peaks_out}")

        # now exporting the unchanging peaks

        rtracklayer::export.bed(row.names(unchanging_reg), con = "${unchanging_peaks_out}")

        # export the other peaks that are greater than 0.05 and less than |1|
        rtracklayer::export.bed(row.names(others_reg), con = "${other_peaks_out}")

        # now hoping to get the peak lengths histone


        ###### if any errors happen here then dont do anything ######
        tryCatch({
        up_peaks = read.table(file = './${up_peaks_out}', header = FALSE, sep = "\t")

        up_peak_lengths = up_peaks\$V3 - up_peaks\$V2
        print(max(up_peak_lengths))

        png(filename = "hist_up_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
        hist.default(up_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        pdf(file = "hist_up_regulated_${idr_histone}_peak_lengths.pdf", height = 8, width = 8)
        hist.default(up_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        # just to view it here, not needed in nextflow here.
        #hist.default(up_peak_lengths, xaxt = "n", breaks = 1e2)
        #axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))

        down_peaks = read.table(file = './${down_peaks_out}', header = FALSE, sep = "\t")

        down_peak_lengths = down_peaks\$V3 - down_peaks\$V2
        print(max(down_peak_lengths))

        png(filename = "hist_down_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
        hist.default(down_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        pdf(file = "hist_down_regulated_${idr_histone}_peak_lengths.pdf", height = 8, width = 8)
        hist.default(down_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()


        unchanging_peaks = read.table(file = './${unchanging_peaks_out}', header = FALSE, sep = "\t")

        unchanging_peak_lengths = unchanging_peaks\$V3 - unchanging_peaks\$V2
        print(max(unchanging_peak_lengths))

        png(filename = "hist_unchanging_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
        hist.default(unchanging_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        pdf(file = "hist_unchanging_regulated_${idr_histone}_peak_lengths.pdf", height = 8, width = 8)
        hist.default(unchanging_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        ########## now violin plots for peak lengths ################

        up_peak_df = DataFrame(peak_lengths_${idr_histone} = up_peak_lengths, category = "up_peaks_lengths")
        down_peak_df = DataFrame(peak_lengths_${idr_histone} = down_peak_lengths, category = "down_peaks_lengths")
        unchanging_peak_df = DataFrame(peak_lengths_${idr_histone} = unchanging_peak_lengths, category = "unchanging_peaks_lengths")


        df_all_peak_lengths = rbind(up_peak_df, down_peak_df, unchanging_peak_df)

        #df_all_peak_lengths

        df_all_peak_lengths_gg =  ggplot(df_all_peak_lengths, aes( category, peak_lengths_${idr_histone}))
        
        all_peak_lengths_violin = df_all_peak_lengths_gg+
            geom_violin()+
            scale_y_continuous(labels = scales::label_number())
        ggsave("${idr_histone}_up_down_unchanging_peak_length_violin_plot.pdf", plot = all_peak_lengths_violin)
        
        
        }, error = function(x) {
        
        message("some of the peak files had no lenght so plotting is pointless")
        })


        tryCatch({
        # now for the annotated peaks, if an error occurs, don't do anything
        # lets get the annotated peaks

        library(ChIPseeker)

        library(TxDb.Hsapiens.UCSC.hg38.knownGene)

        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

        # now to read in the peak files
        peak_files = c(up_peaks = './${up_peaks_out}', down_peaks = './${down_peaks_out}', unchanging_peaks = './${unchanging_peaks_out}')
        up_peak_file = readPeakFile(peakfile = './${up_peaks_out}')


        unchanging_peak_file = readPeakFile(peakfile = './${unchanging_peaks_out}')
        unchanging_annotated_peaks = annotatePeak(peak= unchanging_peak_file, tssRegion = c(-5000, 5000), TxDb = txdb)

        up_annotated_peaks = annotatePeak(peak= up_peak_file, tssRegion = c(-5000, 5000), TxDb = txdb)
        # I would have to plot this pie chart as many times as there are peak files, but not yet.
        #plotAnnoPie(up_annotated_peaks)

        down_peak_file = readPeakFile(peakfile = './${down_peaks_out}')
        down_annotated_peaks = annotatePeak(peak= down_peak_file, tssRegion = c(-5000, 5000), TxDb = txdb)

        # not plotting this right now
        #plotAnnoPie(down_annotated_peaks)


        # not plotting this right now also
        #plotDistToTSS(down_annotated_peaks,
        #            title="Distribution of H3K27me3 peaks relative to TSS")


        annotated_peaks_list = lapply(peak_files, annotatePeak, tssRegion = c(-5000, 5000), TxDb = txdb) 


        # now get the number of total annotated peaks
        down_annotated_df <- as.data.frame(down_annotated_peaks)
        #num_down_annotated_peaks = count(down_annotated_df)[[1]]
        num_down_annotated_peaks = length(down_annotated_df[,1])

        up_annotated_df <- as.data.frame(up_annotated_peaks)
        #num_up_annotated_peaks = count(up_annotated_df)[[1]]
        num_up_annotated_peaks = length(up_annotated_df[,1])

        unchanging_annotated_df <- as.data.frame(unchanging_annotated_peaks)
        #num_unchanging_annotated_peaks = count(unchanging_annotated_df)[[1]]
        num_unchanging_annotated_peaks = length(unchanging_annotated_df[,1])

        # then plot with bar because it uses a ggplot object
        plotAnnoBar(annotated_peaks_list) + ggtitle("Percentage of ${idr_histone} peaks in Genomic Regions", subtitle = paste("Number of down ${idr_histone} peaks annotated ", num_down_annotated_peaks, "\nNumber of up ${idr_histone} peaks annotated ",num_up_annotated_peaks, "\nNumber of unchanging ${idr_histone} peaks annotated ", num_unchanging_annotated_peaks ) )

        ggsave("annotation_bar_graph_${idr_condition}_${idr_histone}_peaks.pdf", plot = last_plot())

        }, error = function(x) {
        
        message("for making annotated peaks, some of the files might have no differential peaks")
        })

        """


    }
    else {

    
        """
        #!/usr/bin/env Rscript

        bam_list = c(${bam_name_list})

        print(bam_list)

        #peaklist = list("\${first_idr_peak}", "\${second_idr_peak}")
        #peaklist = \${peaks_list}

        #print(peaklist)

        #best_hlow_idr
        #best_scrm_idr

        library(DESeq2)
        library(GenomicRanges)
        library(chromVAR)
        library(tidyr)
        library(EnhancedVolcano)
        library(readr)

        print(dir())
        # making the master peak genomic ranges object
        
        #mPeak = GRanges()

        #for (peakfile in c("./\${first_idr_peak}", "./\${second_idr_peak}")) {
        
            #peaktable = read.table(peakfile, header = FALSE, sep = "\t")
        
            #gr_object = GRanges(seqnames = peaktable\$V1, IRanges(start = peaktable\$V2, end = peaktable\$V3), strand = "*" )
        
            #gr_object = rtracklayer::import(peakfile, format = "BED")

            #mPeak = append(mPeak, gr_object)
        #}

        #peaktable = read.table("\${idr_peak_name}", header = FALSE)
        #mPeak = GRanges(seqnames = peaktable\$V1, IRanges(start = peaktable\$V2, end = peaktable\$V3), strand = "*" )

        mPeak = rtracklayer::import("${idr_peak_name}", format = "BED")

        # making sure there are no redundant peaks
        # this will allow me to resizd all of the regions by 50kb on each side, so any regions that are less than 100kb away from eachother will be merged by reduce because they overlap
        
        
        #masterPeak2 = reduce(mPeak)

        # using resize to extnd regions
        extended_master_peak = resize(mPeak, width = width(mPeak)+${resize_num}, fix = "center")

        # before trimming i need to add the sequence lengths of the correct genome
        # will have to find a way to automate this step for other genomes

        library(BSgenome.Hsapiens.UCSC.hg38)
        #seq_lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
        #seqlengths(extended_master_peak) <- seq_lengths[seqlevels(extended_master_peak)]
        #trimmed_master_peak = trim(extended_master_peak)
        #masterPeak_beta = reduce(trimmed_master_peak)

        masterPeak_beta = reduce(extended_master_peak)
        # now to resize again to remove the artificial 50kb from the start and end of the newly reduced peak
        
        new_peak_size = width(masterPeak_beta)-${resize_num}
        masterPeak = resize(masterPeak_beta,  width = new_peak_size , fix = "center")
        
        #print(masterPeak)

        # now I want to keep the standard chromosomes
        #seqnames_to_keep =masterPeak@seqnames@values[1:23]

        #masterPeak = masterPeak[seqnames(masterPeak) %in% seqnames_to_keep]

        # hoping to export my GRanges object master peaks to a bed file
        #write_tsv("\${master_peak_export_out}")
        rtracklayer::export.bed(masterPeak, con = "${master_peak_export_out}")

        

        list_of_mBams = c()
        
        for (bam in bam_list) {
    
            tokens = strsplit(bam, split = "_")[[1]]
            #print(tokens)
            
            histone = tokens[2]
            #print(histone)

            #list_of_mBams = list()

            if ("${idr_histone}" == histone) {
                #print(bam)
                list_of_mBams = c(list_of_mBams, bam)

            }
        }

        print(list_of_mBams)
    
        
        # now doing the next section to make the count matrix and the condition and experiment design

        # making the matrix 

        countsMatrix = matrix(data = NA, length(masterPeak), length(list_of_mBams) )


        # I want to get the chr start end of the peaks and put them in the matrix as column names so I know which peaks I am looking at from the deseq2 results in the differential peaks.
        seqnames(masterPeak)

        # i should put the unique ids from the master peak object in the matrix as row names.
        unique_masterPeak_ids = paste0(as.character(seqnames(masterPeak)), ":" ,start(masterPeak), "-", end(masterPeak))

        rownames(countsMatrix) = unique_masterPeak_ids

        ####################### new version ##########################
        library(Rsubread)
        library(dplyr)

        # featureCounts in rsubread needs a dataframe and in the correct format, with a GeneID column first then chr1 column
        df_masterPeak = as.data.frame(masterPeak)

        # now I need to change the second column name to chr
        df_masterPeak = df_masterPeak %>% rename(chr = seqnames)

        # next i need to add the unique ids as the columns in this dataframe
        df_masterPeak\$GeneID = unique_masterPeak_ids
        ##############################################################

        #getting the list of bam base names to add to the matrix column names
        list_bam_basenames = list()


        # for deseq2 i need the condition design. so hlow and scrm
        # then i will find a way to tally how many of each are there so i can automate the rep count
        condition_design = list()

        type_design = list()

        for (x in c(1:length(list_of_mBams))) {
        
        path_bam = list_of_mBams[[x]]
        print(path_bam)
        bam_basename = basename(as.character(path_bam))
        bam_tokens = strsplit(basename(as.character(path_bam)), split = "_")[[1]]
        
        # labeling the important tokens so it is easier to keep track of
        condition = bam_tokens[1]
        histone = bam_tokens[2]
        replicate = bam_tokens[3]
        
        
        # for later parts I need the condition and to know how many times to repeat them
        condition_design = append(condition_design, paste(condition,histone,sep="_"))
        
        
        # also get the replicate too for type design
        type_design = append(type_design, replicate)
        #type_design = append(type_design, paste(histone,replicate, sep="_"))
        
        
        
        # using chromVAR getCounts function
        # not using chromVAR anymore, for now
        #fragment_counts = getCounts(path_bam, masterPeak, paired= TRUE, by_rg = FALSE, format = "bam")

        fragment_counts <- featureCounts(
            files = path_bam,
            annot.ext = df_masterPeak,
            isPairedEnd = FALSE,      # TRUE if your data were paired-end
            strandSpecific = 0,       # 0 = unstranded, 1/2 if strand-specific
            useMetaFeatures = TRUE
        )
        
        # putting the fragment counts in the column labeled by the bam name
        countsMatrix[,x] = fragment_counts\$counts
        
        list_bam_basenames = append(list_bam_basenames, bam_basename)
        
        }

        colnames(countsMatrix) = list_bam_basenames

        ## this part i would have to tell the use to name their conditions as control and treatment so it can make the design treatment vs control

        #first removing the low count fragments. only keeping the rows that are above 5 count in total

        keep_Rows = which(rowSums(countsMatrix) > 5)

        filt_countmatrix = countsMatrix[keep_Rows,]

        # now to get the condition_design
        condition_counts = table(unlist(condition_design))

        # this gives back the names and the counts give back the counts for each of the names
        condition_names = names(condition_counts)
        condition_num = as.numeric(condition_counts)




        # now i can put both lists in to get back the experiment design
        #condition_factor = factor(rep(condition_names, times=condition_num))

        # Build condition factor in column order
        condition_factor <- factor(unlist(condition_design))

        # I want it so the treatment is second in the list here so it is first in the experimental design in deseq2. This way the experimental design will have the diff results going up or down in the treatment vs control

        # for treatment I need to make nextflow take the users input so they can relevel this section
        #treatment = "H1low_H3k27me3"
        treatment = "${params.treatment_name}_${idr_histone}"

        # this doesnt work because you need to put ref = levels(condition_factor)[2] to assign the actual relevel
        #if (levels(condition_factor)[1] == treatment ) {
        #condition_factor = relevel(condition_factor, levels(condition_factor)[2])
        #}else {
        #condition_factor
        #}

        # have to do this instead to make the treatment always the baseline
        #condition_factor <- relevel(condition_factor, ref = treatment)

        levels_now <- levels(condition_factor)
        control <- levels_now[ levels_now != treatment ]  # the other level

        # Set control as reference → LFC = treatment / control
        condition_factor <- relevel(condition_factor, ref = control)


        print(condition_factor)

        # repeating the above to have another column with type (replicates)
        #type_counts = table(unlist(type_design))
        #type_names = names(type_counts)
        #type_num = as.numeric(type_counts)

        #type_factor = factor(rep(type_names, times=type_counts))
        #type_factor = factor(rep(type_names, times=type_num[1]))

        # Build type factor in column order
        type_factor <- factor(unlist(type_design))

        # I want to get the idr threshold and use that as input for the file names and other things

        #peak_file = basename(peaklist[[1]])
        peak_file = basename("./${idr_peak_name}")

        idr_used = strsplit(peak_file, split = "_")[[1]][10]
        idr_used

        # now for deseq2 workflow


        # now to do the normal deseq2 workflow

        dds = DESeqDataSetFromMatrix(countData = filt_countmatrix,
                                    colData = DataFrame(condition_factor, type_factor),
                                    design = ~ condition_factor)


        # using the function on our data
        DDS = DESeq(dds)

        norm_DDS = counts(DDS, normalized = TRUE) # normalization with respect to the sequencing depth

        # adding _norm onto the column names in the normalized matrix
        colnames(norm_DDS) = paste0(colnames(norm_DDS), "_norm")


        # provides independent filtering using the mean of normalized counts
        res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")


        # this is looking at the differences between the 3 deseq analyzed options
        countMatDiff = cbind(filt_countmatrix, norm_DDS, res)

        head(countMatDiff)




        # getting the results name and addding to the coef we want to shrink
        experiment_design_name = resultsNames(DDS)[2]

        # useful for visualization and ranking of genes or in this case peaks
        resLFC = lfcShrink(DDS, coef= resultsNames(DDS)[2], type = "apeglm")


        # finding the up and down regulated counts that pass the threshold

        

        up_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange >= 1) ,]
        
        
        down_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange <= -1) ,]
        

        #unchanging_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange <1 & resLFC\$log2FoldChange > -1) ,]
        unchanging_reg = resLFC[which(resLFC\$log2FoldChange <1 & resLFC\$log2FoldChange > -1) ,]

        #others_reg = resLFC[which(resLFC\$padj > 0.05), ] 
        others_reg = resLFC[which(resLFC\$padj > 0.05 & resLFC\$log2FoldChange <1 & resLFC\$log2FoldChange > -1), ]        

        # testing the chat gpt code to make the plot look publication ready

        # Relabel categories (for peaks, not genes)
        resLFC\$Label <- ifelse(resLFC\$padj > 0.05, "Not significant",
                            ifelse(resLFC\$log2FoldChange >= 1 & resLFC\$padj < 0.05, "Upregulated peaks",
                            ifelse(resLFC\$log2FoldChange <= -1 & resLFC\$padj < 0.05, "Downregulated peaks",
                            ifelse(resLFC\$log2FoldChange >= 1 & resLFC\$padj > 0.05, "Not significant",
                            ifelse(resLFC\$log2FoldChange <= -1 & resLFC\$padj > 0.05, "Not significant",
                                    "padj < 0.05 only")))))

        # Counts & percentages
        total <- nrow(resLFC)
        #up_count <- sum(resLFCLabel == "Upregulated peaks")
        #down_count <- sum(resLFCLabel == "Downregulated peaks")
        #ns_count <- sum(resLFCLabel == "Not significant")
        #padj_only_count <- sum(resLFCLabel == "padj < 0.05 only")

        up_pct <- round(100 * nrow(up_reg) / total, 1)
        down_pct <- round(100 * nrow(down_reg) / total, 1)
        ns_pct <- round(100 * nrow(others_reg) / total, 1)
        padj_only_pct <- round(100 * nrow(unchanging_reg) / total, 1)

        # Annotation text
        annotation_text <- paste0(
        "Upregulated peaks: ", nrow(up_reg), " (", up_pct, "%)\n",
        "Downregulated peaks: ", nrow(down_reg), " (", down_pct, "%)\n",
        "padj < 0.05 only: ", nrow(unchanging_reg), " (", padj_only_pct, "%)\n",
        "Not significant: ", nrow(others_reg), " (", ns_pct, "%)"
        )

        # Plot
        ma_plot_labeled <- ggplot(resLFC, aes(x = baseMean, y = log2FoldChange, color = Label)) +
        geom_point(alpha = 0.7, size = 0.5) +
        geom_hline(yintercept = 0, linetype = "solid", color = "black") +
        geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black") +
        scale_x_continuous(trans = "log10", 
                            breaks = scales::trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        scale_colour_manual(values = c("Upregulated peaks" = "red",
                                        "Downregulated peaks" = "red",
                                        "padj < 0.05 only" = "lightpink",
                                        "Not significant" = "grey70"),
                            breaks = c("Upregulated peaks", "Downregulated peaks", "padj < 0.05 only", "Not significant")) +
        labs(title = experiment_design_name,
            subtitle = paste("IDR threshold for masterPeaks:", "${params.return_idr}"),
            x = "Mean peak signal (log10 scale)",
            y = expression(Log[2]~Fold~Change),
            color = "Peak status") +
        theme_classic(base_size = 4) +
        theme(legend.position = "top",
                legend.title = element_text(size = 4),
                legend.text = element_text(size = 5),
                plot.title = element_text(size = 6, face = "bold"),
                plot.subtitle = element_text(size = 6),
                axis.title = element_text(size = 6),
                axis.text = element_text(size = 6)) +
        annotate("text", x = min(resLFC\$baseMean, na.rm = TRUE), 
                y = max(resLFC\$log2FoldChange, na.rm = TRUE), 
                hjust = -3, vjust = 1, 
                label = annotation_text, 
                size = 1.5)



        print(ma_plot_labeled)


        pdf(file = paste(experiment_design_name,"mergedBy", "${merge_dist}","MA_plot_our_data_counts.pdf", sep = "_"), width = 4, height = 4)

        print(ma_plot_labeled)
        dev.off()

        png(filename = paste(experiment_design_name,"mergedBy", "${merge_dist}","MA_plot_our_data_counts.png", sep = "_"), width = 4.3, height = 4, antialias = "subpixel", units = "in", res = 2000)

        print(ma_plot_labeled)
        dev.off()


        volcano_plot_removed_reps = EnhancedVolcano(resLFC,
                        lab = rownames(resLFC),
                        title = experiment_design_name, subtitle = paste("This means the first condition has these points more than in the second condition. The idr threshold for masterPeaks is",idr_used, sep = " "),
                        x = 'log2FoldChange', FCcutoff = 1,
                        y = 'padj', pCutoff = 0.05, labSize = 0 # making this 0 so it doesn't show the region names in the volcano plot 
                        
                        )

        print(volcano_plot_removed_reps)

        png(filename = paste(experiment_design_name,"mergedBy","${merge_dist}","volcano_plot_our_data_counts.png", sep = "_"), width = 1200, height = 900, antialias = "subpixel")
        print(volcano_plot_removed_reps)
        dev.off()

        pdf(file = paste(experiment_design_name,"mergedBy","${merge_dist}","volcano_plot_our_data_counts.pdf", sep = "_"), width = 16, height = 12)
        print(volcano_plot_removed_reps)
        dev.off()
        



        # using rlog over vst for transformation

        rld = rlog(DDS, blind=FALSE)

        #it is clear that the rlog transformation inherently accounts for differences in sequencing depth *
        head(assay(rld), 5)


        library(ggplot2)

        pcaData = plotPCA(rld, intgroup=c("condition_factor","type_factor"), returnData = TRUE )
        pcaData

        percentVar <- round(100 * attr(pcaData, "percentVar"))


        pca_plot_rlog = ggplot(pcaData, aes(PC1, PC2, color=condition_factor, shape=type_factor)) +
        ggtitle(paste("PCA plot using Rlog transform IDR", idr_used, sep = " ")) +
        theme(plot.caption = element_text(size = 30, face = "bold"), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15))+
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed()
        #pca_plot
        name_for_rlog_png_file =paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used, "rlog.png", sep = "_")

        png(filename = name_for_rlog_png_file, width = 900, height = 900, antialias = "subpixel")
        print(pca_plot_rlog)
        dev.off()

        name_for_rlog_pdf_file =paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used, "rlog.pdf", sep = "_")

        pdf(file = name_for_rlog_pdf_file, width = 10, height = 10)
        print(pca_plot_rlog)
        dev.off()


        # testing with vst
        #vsd_t = vst(DDS, blind = FALSE)

        # for histone mark k36me2 vst fails so i have to use the direct function
        vsd_t <- tryCatch({
            vst(DDS, blind = FALSE)
        }, error = function(e) {
            message("vst() failed, using varianceStabilizingTransformation() instead.")
            varianceStabilizingTransformation(DDS, blind = FALSE)
        })

        head(assay(vsd_t), 5)

        pcaData2 = plotPCA(vsd_t, intgroup=c("condition_factor","type_factor"), returnData = TRUE )
        pcaData2

        percentVar <- round(100 * attr(pcaData2, "percentVar"))
        pca_plot_vst = ggplot(pcaData2, aes(PC1, PC2, color=condition_factor, shape=type_factor)) +
        ggtitle(paste("PCA plot using VST transform IDR", idr_used, sep = " ")) +
        theme(plot.caption = element_text(size = 30, face = "bold"), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15))+
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed()


        #name_for_vst_png_file =paste(experiment_design_name,"PCA_plot_IDR", idr_used, "vst.png", sep = "_")

        pdf(file = paste(experiment_design_name,"PCA_plot_IDR_our_data_counts",idr_used,"vst_all_reps.pdf", sep = "_"), width = 10, height = 10)
        print(pca_plot_vst)
        dev.off()

        png(file = paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used,"vst_all_reps.png", sep = "_"), width = 900, height = 900)
        print(pca_plot_vst)
        dev.off()

        pca_plot_rlog

        pca_plot_vst


        rtracklayer::export.bed(row.names(up_reg), con = "${up_peaks_out}")

        rtracklayer::export.bed(row.names(down_reg), con = "${down_peaks_out}")

        # now exporting the unchanging peaks

        rtracklayer::export.bed(row.names(unchanging_reg), con = "${unchanging_peaks_out}")

        # export the other peaks that are greater than 0.05 and less than |1|
        rtracklayer::export.bed(row.names(others_reg), con = "${other_peaks_out}")

        # now hoping to get the peak lengths histone


        ###### if any errors happen here then dont do anything ######
        tryCatch({
        up_peaks = read.table(file = './${up_peaks_out}', header = FALSE, sep = "\t")

        up_peak_lengths = up_peaks\$V3 - up_peaks\$V2
        print(max(up_peak_lengths))

        png(filename = "hist_up_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
        hist.default(up_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        pdf(file = "hist_up_regulated_${idr_histone}_peak_lengths.pdf", height = 8, width = 8)
        hist.default(up_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        # just to view it here, not needed in nextflow here.
        #hist.default(up_peak_lengths, xaxt = "n", breaks = 1e2)
        #axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))

        down_peaks = read.table(file = './${down_peaks_out}', header = FALSE, sep = "\t")

        down_peak_lengths = down_peaks\$V3 - down_peaks\$V2
        print(max(down_peak_lengths))

        png(filename = "hist_down_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
        hist.default(down_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        pdf(file = "hist_down_regulated_${idr_histone}_peak_lengths.pdf", height = 8, width = 8)
        hist.default(down_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()


        unchanging_peaks = read.table(file = './${unchanging_peaks_out}', header = FALSE, sep = "\t")

        unchanging_peak_lengths = unchanging_peaks\$V3 - unchanging_peaks\$V2
        print(max(unchanging_peak_lengths))

        png(filename = "hist_unchanging_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
        hist.default(unchanging_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        pdf(file = "hist_unchanging_regulated_${idr_histone}_peak_lengths.pdf", height = 8, width = 8)
        hist.default(unchanging_peak_lengths, xaxt = "n", breaks = 1e2)
        axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
        dev.off()

        ########## now violin plots for peak lengths ################

        up_peak_df = DataFrame(peak_lengths_${idr_histone} = up_peak_lengths, category = "up_peaks_lengths")
        down_peak_df = DataFrame(peak_lengths_${idr_histone} = down_peak_lengths, category = "down_peaks_lengths")
        unchanging_peak_df = DataFrame(peak_lengths_${idr_histone} = unchanging_peak_lengths, category = "unchanging_peaks_lengths")


        df_all_peak_lengths = rbind(up_peak_df, down_peak_df, unchanging_peak_df)

        #df_all_peak_lengths

        df_all_peak_lengths_gg =  ggplot(df_all_peak_lengths, aes( category, peak_lengths_${idr_histone}))
        
        all_peak_lengths_violin = df_all_peak_lengths_gg+
            geom_violin()+
            scale_y_continuous(labels = scales::label_number())
        ggsave("${idr_histone}_up_down_unchanging_peak_length_violin_plot.pdf", plot = all_peak_lengths_violin)
        
        
        }, error = function(x) {
        
        message("some of the peak files had no lenght so plotting is pointless")
        })


        tryCatch({
        # now for the annotated peaks, if an error occurs, don't do anything

        # lets get the annotated peaks

        library(ChIPseeker)

        library(TxDb.Hsapiens.UCSC.hg38.knownGene)

        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

        # now to read in the peak files
        peak_files = c(up_peaks = './${up_peaks_out}', down_peaks = './${down_peaks_out}', unchanging_peaks = './${unchanging_peaks_out}')
        up_peak_file = readPeakFile(peakfile = './${up_peaks_out}')


        unchanging_peak_file = readPeakFile(peakfile = './${unchanging_peaks_out}')
        unchanging_annotated_peaks = annotatePeak(peak= unchanging_peak_file, tssRegion = c(-5000, 5000), TxDb = txdb)

        up_annotated_peaks = annotatePeak(peak= up_peak_file, tssRegion = c(-5000, 5000), TxDb = txdb)
        # I would have to plot this pie chart as many times as there are peak files, but not yet.
        #plotAnnoPie(up_annotated_peaks)

        down_peak_file = readPeakFile(peakfile = './${down_peaks_out}')
        down_annotated_peaks = annotatePeak(peak= down_peak_file, tssRegion = c(-5000, 5000), TxDb = txdb)

        # not plotting this right now
        #plotAnnoPie(down_annotated_peaks)


        # not plotting this right now also
        #plotDistToTSS(down_annotated_peaks,
        #            title="Distribution of H3K27me3 peaks relative to TSS")


        annotated_peaks_list = lapply(peak_files, annotatePeak, tssRegion = c(-5000, 5000), TxDb = txdb) 


        # now get the number of total annotated peaks
        down_annotated_df <- as.data.frame(down_annotated_peaks)
        #num_down_annotated_peaks = count(down_annotated_df)[[1]]
        num_down_annotated_peaks = length(down_annotated_df[,1])

        up_annotated_df <- as.data.frame(up_annotated_peaks)
        #num_up_annotated_peaks = count(up_annotated_df)[[1]]
        num_up_annotated_peaks = length(up_annotated_df[,1])

        unchanging_annotated_df <- as.data.frame(unchanging_annotated_peaks)
        #num_unchanging_annotated_peaks = count(unchanging_annotated_df)[[1]]
        num_unchanging_annotated_peaks = length(unchanging_annotated_df[,1])

        # then plot with bar because it uses a ggplot object
        plotAnnoBar(annotated_peaks_list) + ggtitle("Percentage of ${idr_histone} peaks in Genomic Regions", subtitle = paste("Number of down ${idr_histone} peaks annotated ", num_down_annotated_peaks, "\nNumber of up ${idr_histone} peaks annotated ",num_up_annotated_peaks, "\nNumber of unchanging ${idr_histone} peaks annotated ", num_unchanging_annotated_peaks ) )

        ggsave("annotation_bar_graph_${idr_condition}_${idr_histone}_peaks.pdf", plot = last_plot())

        }, error = function(x) {
        
        message("for making annotated peaks, some of the files might have no differential peaks")
        })

        """
    }
    
    
}