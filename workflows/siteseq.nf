/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPAREGENOME                          } from '../subworkflows/local/preparegenome'
include { PREPROCESSREADS                        } from '../subworkflows/local/preprocessreads'
include { ALIGNREADS                             } from '../subworkflows/local/alignreads'
include { SITESEQPROCESSING                      } from '../subworkflows/local/siteseqprocessing'
include { SITESEQREPORTING                       } from '../subworkflows/local/siteseqreporting'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                 } from '../subworkflows/local/utils_nfcore_siteseq_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SITESEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Pre-process reads
    //
    PREPROCESSREADS(
        ch_samplesheet,
        params.subsample_fastqs,
        params.subsample_target_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESSREADS.out.fastp_json.collect{it[1]})
    ch_versions = ch_versions.mix(PREPROCESSREADS.out.versions)

    //
    // SUBWORKFLOW: Prepare genome
    //
    PREPAREGENOME(params.fasta, params.bowtie2_index)
    ch_versions = ch_versions.mix(PREPAREGENOME.out.versions)

    //
    // SUBWORKFLOW: Align reads
    //
    ALIGNREADS(
        PREPROCESSREADS.out.reads,
        PREPAREGENOME.out.bowtie2_index,
        PREPAREGENOME.out.fasta,
        params.use_umis
    )
    ch_versions = ch_versions.mix(ALIGNREADS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNREADS.out.bowtie2_log.collect{it[1]})
    ch_multiqc_files.mix(ALIGNREADS.out.umitools_dedup_log.collect{it[1]})

    //
    // SUBWORKFLOW: SITE-Seq® assay processing (calculate stats, call sites)
    //
    SITESEQPROCESSING(
        PREPROCESSREADS.out.reads,
        ALIGNREADS.out.bam,
        PREPAREGENOME.out.fasta,
        PREPAREGENOME.out.fai,
        params.apply_normalization,
        params.min_aln_length,
        params.min_map_qual,
        params.max_clipping,
        params.stat_window,
        params.control_noise_window,
        params.peak_window,
        params.agg_peak_window,
        params.test_pval_reporting_thresh,
        params.test_pval_calling_thresh,
        params.ctrl_noise_pval_calling_thresh,
        params.include_noise_pval_fails,
        params.motif_search_factor,
        params.motif_match_score,
        params.motif_mismatch_pen,
        params.motif_gap_pen
    )

    //
    // SUBWORKFLOW: SITE-Seq® assay reporting (create final outputs)
    //
    SITESEQREPORTING(
        PREPROCESSREADS.out.reads,
        PREPROCESSREADS.out.fastp_json,
        PREPROCESSREADS.out.umi_filter_report,
        ALIGNREADS.out.bam,
        ALIGNREADS.out.bowtie2_log,
        ALIGNREADS.out.umitools_dedup_log,
        SITESEQPROCESSING.out.count_reads_json,
        PREPAREGENOME.out.fasta,
        SITESEQPROCESSING.out.agg_site_tsv,
        SITESEQPROCESSING.out.window_diff_bigwig,
        SITESEQPROCESSING.out.peaks_dist_tsv,
        SITESEQPROCESSING.out.ctrl_peaks_dist_tsv,
        SITESEQPROCESSING.out.ctrl_wide_sums_dist_tsv,
        params.test_pval_reporting_thresh,
        params.test_pval_calling_thresh,
        params.ctrl_noise_pval_calling_thresh,
        params.bam_dir,
        params.bigwig_dir,
        params.fasta_dir,
        params.igv_dir,
        params.file_server_url,
        params.file_server_prefixes?.tokenize(","),
        params.use_relative_igv_paths
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  ''  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
