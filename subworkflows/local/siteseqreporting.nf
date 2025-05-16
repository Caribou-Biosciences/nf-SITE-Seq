
include { COMPILESTATS                    } from '../../modules/local/compilestats'
include { IGV                             } from '../../modules/local/igv'
include { CREATEHTMLREPORT                } from '../../modules/local/createhtmlreport'

include { createReplicateGroups           } from './utils_nfcore_siteseq_pipeline'
include { getReplicateGroupSizes          } from './utils_nfcore_siteseq_pipeline'
include { getControlReplicateGroupLimits  } from './utils_nfcore_siteseq_pipeline'
include { createSampleGroups              } from './utils_nfcore_siteseq_pipeline'
include { getSampleGroupSizes             } from './utils_nfcore_siteseq_pipeline'


//
// Utility function for constructing a path to a file resource for use in an IGV
// session file. Supports absolute paths, relative paths, and file server URLs.
//
def getIgvFilePath(
    file, file_dir_path, igv_dir_path, file_server_url, file_server_prefixes, use_relative_paths
) {
    def file_path = java.nio.file.Paths.get(file_dir_path, file.name)
    def igv_dir = java.nio.file.Paths.get(igv_dir_path)
    def igv_file_path
    if (file_server_url) {
        if (file_server_prefixes.size() > 0) {
            igv_file_path = file_path.toString()
            file_server_prefixes.collect { file_server_prefix ->
                igv_file_path = igv_file_path.replace(file_server_prefix, file_server_url)
            }
        } else {
            igv_file_path = java.nio.file.Paths.get(file_server_url, file_path.toString())
        }
    } else if (use_relative_paths) {
        igv_file_path = igv_dir.relativize(file_path)
    } else {
        igv_file_path = file_path.toString()
    }
    return igv_file_path
}


process CONCATENATE_TSV_WITH_HEADER {
    input:
    path tsv_files

    output:
    path "all_processing_metrics.tsv", emit: tsv

    script:
    def sorted_tsv_files = tsv_files.sort { it.name }
    """
    # Extract header from the first file
    head -n 1 ${sorted_tsv_files[0]} > all_processing_metrics.tsv

    # Concatenate the content of all files, excluding the header
    for file in $sorted_tsv_files
    do
        tail -n +2 \$file >> all_processing_metrics.tsv
    done
    """
}


workflow SITESEQREPORTING {

    take:
    ch_reads                       // channel: [ val(meta), [ fastq  ] ]
    ch_fastp_json                  // channel: [ val(meta), [ json   ] ]
    ch_umi_filter_report           // channel: [ val(meta), [ json   ] ]
    ch_bam                         // channel: [ val(meta), [ bam    ] ]
    ch_bowtie2_log                 // channel: [ val(meta), [ log    ] ]
    ch_umitools_dedup_log          // channel: [ val(meta), [ log    ] ]
    ch_count_reads_json            // channel: [ val(meta), [ json   ] ]
    ch_fasta                       // channel: [ val(meta), [ fasta  ] ]
    ch_agg_site_tsv                // channel: [ val(meta), [ tsv    ] ]
    ch_window_diff_bigwig          // channel: [ val(meta), [ bigwig ] ]
    ch_peaks_dist_tsv              // channel: [ val(meta), [ tsv    ] ]
    ch_ctrl_peaks_dist_tsv         // channel: [ val(meta), [ tsv    ] ]
    ch_ctrl_wide_sums_dist_tsv     // channel: [ val(meta), [ tsv    ] ]

    test_pval_reporting_thresh     // float: p-val threshold below which sites will be reported
    test_pval_calling_thresh       // float: p-val threshold below which sites will be called
    ctrl_noise_pval_calling_thresh // float: control noise p-val threshold below which sites will be filtered

    bam_dir                        // value: path to the parent folder where BAMs will be published
    bigwig_dir                     // value: path to the parent folder where BigWigs will be published
    fasta_dir                      // value: path to the parent folder where genome FASTA will be published
    igv_dir                        // value: path to the parent folder where IGV session files will be published
    file_server_url                // value: base URL of file server
    file_server_prefixes           // [value]: path prefixes which are indexed on file server
    use_relative_igv_paths         // value: whether to use relative file paths in IGV session files

    main:

    ch_versions = Channel.empty()

    // Group peak dists by sample group
    ch_meta = ch_reads.map { meta, _fastq -> meta }
    ch_sample_group_sizes = getSampleGroupSizes(ch_meta)
    ch_grouped_peaks_dists = createSampleGroups(ch_peaks_dist_tsv, ch_sample_group_sizes)

    // Gather various per-sample logs/reports and create a single-row TSV with
    // read/UMI/alignment stats
    ch_stats_inputs = ch_fastp_json.map{ meta, files -> [ meta.id, meta, files] }
        .join(ch_umi_filter_report.map{ meta, files -> [ meta.id, files] }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_bowtie2_log.map{ meta, files -> [ meta.id, files] }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_umitools_dedup_log.map{ meta, files -> [ meta.id, files] }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_count_reads_json.map{ meta, files -> [ meta.id, files] }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { it -> it[1..it.size()-1] }
    COMPILESTATS(ch_stats_inputs)

    // Concatenate per-sample stat TSVs into a single TSV
    CONCATENATE_TSV_WITH_HEADER(COMPILESTATS.out.stats_tsv.map { _meta, tsv -> tsv }.collect())

    // Create HTML report summarizing sample group metadata, site calling summary,
    // threshold summary, and called/reported sites
    ch_report_inputs = ch_agg_site_tsv
        .join(ch_grouped_peaks_dists, by: [0], failOnMismatch: true, failOnDuplicate: true)
        .map { meta, agg_site_tsv, peak_dist_tsvs -> [ meta.control_group, meta, agg_site_tsv, peak_dist_tsvs ] }
        .combine(
            ch_ctrl_peaks_dist_tsv.map { meta, tsv -> [ meta.sample_group, tsv ] },
            by: 0
        )
        .combine(
            ch_ctrl_wide_sums_dist_tsv.map { meta, tsv -> [ meta.sample_group, tsv ] },
            by: 0
        )
        .map { _ctrl_group, meta, agg_site_tsv, peak_dists, ctrl_peaks_dist, wide_sums_dist ->
            [ meta, agg_site_tsv, peak_dists, ctrl_peaks_dist, wide_sums_dist ]
        }
    CREATEHTMLREPORT(
        ch_report_inputs,
        CONCATENATE_TSV_WITH_HEADER.out.tsv,
        test_pval_reporting_thresh,
        test_pval_calling_thresh,
        ctrl_noise_pval_calling_thresh
    )

    // Gather BAMs and BigWigs into replicate groups (same sample group and concentration)
    // and sample groups (same sample group)
    ch_replicate_group_sizes = getReplicateGroupSizes(ch_meta)
    ch_control_group_limits = getControlReplicateGroupLimits(ch_meta)
    ch_rep_group_bams = createReplicateGroups(
        ch_bam, ch_replicate_group_sizes, ch_control_group_limits
    )
    ch_test_rep_group_bams = ch_rep_group_bams.filter { meta, _files -> !meta.is_control }
    ch_control_rep_group_bams = ch_rep_group_bams.filter { meta, _files -> meta.is_control }
    ch_test_sample_group_bams = createSampleGroups(ch_test_rep_group_bams, ch_sample_group_sizes)
    ch_test_sample_group_bigwigs = createSampleGroups(ch_window_diff_bigwig, ch_sample_group_sizes)

    // Convert BAM/BigWig paths to file server URLs (if file_server_url is defined)
    // and group into IGV panels to create an IGV session file
    ch_igv_panels = ch_test_sample_group_bams
        .map { meta, bams -> [ meta.sample_group, meta, bams ] }
        .join(
            ch_test_sample_group_bigwigs.map { meta, bigwigs -> [ meta.sample_group, bigwigs ] },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true
        )
        .map { _sample_group, meta, bams, bigwigs -> [ meta.control_group, meta, bams, bigwigs ] }
        .combine(
            ch_control_rep_group_bams.map { meta, bams -> [ meta.sample_group, bams ]},
            by: 0
        )
        .map { _control_group, meta, test_bams, test_bigwigs, control_bams ->
            def test_bam_urls = test_bams
                .flatten()
                .collect { bam -> [getIgvFilePath(bam, bam_dir, igv_dir, file_server_url, file_server_prefixes, use_relative_igv_paths)] }
            def test_bigwig_urls = test_bigwigs
                .flatten()
                .collect { bigwig -> [getIgvFilePath(bigwig, bigwig_dir, igv_dir, file_server_url, file_server_prefixes, use_relative_igv_paths)] }
            def control_bam_urls = control_bams
                .collect { bam ->
                    [getIgvFilePath(bam, bam_dir, igv_dir, file_server_url, file_server_prefixes, use_relative_igv_paths)]
                }
            return [ meta, [test_bigwig_urls] + test_bam_urls + control_bam_urls ]
        }
    ch_fasta_url = ch_fasta
        .map { meta, fasta ->
            [ meta, getIgvFilePath(fasta, fasta_dir, igv_dir, file_server_url, file_server_prefixes, use_relative_igv_paths) ]
        }
    IGV(
        ch_igv_panels,
        ch_fasta_url
    )

    emit:
    html     = CREATEHTMLREPORT.out.html       // channel: [ val(meta), [ html ] ]
    igv      = IGV.out.xml                     // channel: [ val(meta), [ xml ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

