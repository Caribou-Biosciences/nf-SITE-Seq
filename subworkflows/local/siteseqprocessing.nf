include { COUNTREADSTARTS                                     } from '../../modules/local/countreadstarts'
include { MERGEREADSTARTS                                     } from '../../modules/local/mergereadstarts'
include { CALCULATEWINDOWDIFFS                                } from '../../modules/local/calculatewindowdiffs'
include { CALCULATEWINDOWDIFFS as CALCULATECTRLWINDOWDIFFS    } from '../../modules/local/calculatewindowdiffs'
include { CALCULATEWINDOWSUMS                                 } from '../../modules/local/calculatewindowsums'
include { CALCULATEWINDOWSUMS as CALCULATECTRLWINDOWSUMS      } from '../../modules/local/calculatewindowsums'
include { CALCULATEWINDOWSUMS as CALCULATECTRLWIDEWINDOWSUMS  } from '../../modules/local/calculatewindowsums'
include { EXTRACTWINDOWEDPEAKS                                } from '../../modules/local/extractwindowedpeaks'
include { EXTRACTWINDOWEDPEAKS as EXTRACTCTRLWINDOWEDPEAKS    } from '../../modules/local/extractwindowedpeaks'
include { REPORTANDCALLSITES                                  } from '../../modules/local/reportandcallsites'
include { AGGREGATESITES                                      } from '../../modules/local/aggregatesites'
include { TSVTOBIGWIG as TSVTOBIGWIGDIFFS                     } from '../../modules/local/tsvtobigwig'

include { TABIX_TABIX as TABIX_MERGEDSTARTS                   } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_WINDOWDIFFS                    } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_WINDOWSUMS                     } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_CTRLWINDOWSUMS                 } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_CTRLWIDEWINDOWSUMS             } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_WINDOWEDPEAKS                  } from '../../modules/nf-core/tabix/tabix/main'

include { createReplicateGroups                               } from './utils_nfcore_siteseq_pipeline'
include { getReplicateGroupSizes                              } from './utils_nfcore_siteseq_pipeline'
include { getControlReplicateGroupLimits                      } from './utils_nfcore_siteseq_pipeline'
include { createSampleGroups                                  } from './utils_nfcore_siteseq_pipeline'
include { getSampleGroupSizes                                 } from './utils_nfcore_siteseq_pipeline'
include { createReplicateGroupCombos                          } from './utils_nfcore_siteseq_pipeline'
include { getControlReplicateGroupSizes                       } from './utils_nfcore_siteseq_pipeline'

workflow SITESEQPROCESSING {

    take:
    ch_reads                       // channel: [ val(meta), [ fastq ] ]
    ch_bam                         // channel: [ val(meta), [ bam   ] ]
    ch_fasta                       // channel: [ val(meta), [ fasta ] ]
    ch_fai                         // channel: [ val(meta), [ fai   ] ]
    apply_normalization            // boolean: Normalize read counts by reads per million
    min_aln_length                 // integer: minimum reference length for an alignment to be included
    min_map_qual                   // integer: minimum mapping quality for an alignment to be included
    max_clipping                   // integer: maximum clipping at read start for an alignment to be included
    stat_window                    // integer: window size (in bp) used to calculate test statistic
    control_noise_window           // integer: window size (in bp) used to calculate control noise
    peak_window                    // integer: window size (in bp) used to collapse sites into peaks
    agg_peak_window                // integer: window size (in bp) use to collapse peaks across concentrations

    test_pval_reporting_thresh     // float: p-val threshold below which sites will be reported
    test_pval_calling_thresh       // float: p-val threshold below which sites will be called
    ctrl_noise_pval_calling_thresh // float: control noise p-val threshold below which sites will be filtered
    include_noise_pval_fails       // float: include control noise p-vals below threshold (still uncalled)

    motif_search_factor            // float: search for motifs within motif_search_factor * <length of on-target motif>
    motif_match_score              // integer: score for matching a motif base
    motif_mismatch_pen             // integer: penalty for mismatching a motif base
    motif_gap_pen                  // integer: penalty for a gap in motif

    main:

    ch_versions = Channel.empty()

    // Count the number of reads originating from each position in the genome
    COUNTREADSTARTS(
        ch_bam,
        apply_normalization,
        min_aln_length,
        min_map_qual,
        max_clipping
    )

    // Merge read starts by replicate group (replicates of the same sample group and RNP concentration)
    // Figure out how many replicates there are per group (allows groups to emit as they're ready)
    // Figure out how many control replicates to use per group (same number as referencing test groups)
    ch_meta = ch_reads.map { meta, _fastq -> meta }
    ch_replicate_group_sizes = getReplicateGroupSizes(ch_meta)
    ch_control_group_limits = getControlReplicateGroupLimits(ch_meta)
    ch_replicate_groups = createReplicateGroups(
        COUNTREADSTARTS.out.tsv, ch_replicate_group_sizes, ch_control_group_limits
    )
    MERGEREADSTARTS(ch_replicate_groups)
    ch_merged_counts = MERGEREADSTARTS.out.tsv
        .branch { meta, _tsv ->
            test: !meta.is_control
            control: meta.is_control
        }

    // Calculate sliding window sums of merged read starts for both test and control
    // using stat_window as window size
    CALCULATEWINDOWSUMS(ch_merged_counts.test, ch_fai, stat_window)
    CALCULATECTRLWINDOWSUMS(ch_merged_counts.control, ch_fai, stat_window)

    // Calculate sliding window sums of control merged read starts using control_noise_window
    // as window size (a wider window used to identify high-background regions)
    CALCULATECTRLWIDEWINDOWSUMS(ch_merged_counts.control, ch_fai, control_noise_window)

    // Create combinations of merged read starts for test and control replicate groups
    ch_replicate_group_combos = createReplicateGroupCombos(MERGEREADSTARTS.out.tsv)

    // Calculate window diffs (i.e. signal statistic)
    CALCULATEWINDOWDIFFS(ch_replicate_group_combos, ch_fai, stat_window, false)

    // Create BigWig of window diffs for use in IGV
    TSVTOBIGWIGDIFFS(CALCULATEWINDOWDIFFS.out.tsv, ch_fai)

    // Identify peaks of window diff statistic using sliding window
    EXTRACTWINDOWEDPEAKS(CALCULATEWINDOWDIFFS.out.tsv, peak_window)

    // Create spoofed replicate group combos and calculate control vs control window diffs
    // (we're comparing control to itself)
    ch_ctrl_replicate_group_combos = MERGEREADSTARTS.out.tsv
        .filter { meta, _files -> meta.is_control }
        .map { meta, files -> [ meta, files, null, [] ] }
    CALCULATECTRLWINDOWDIFFS(ch_ctrl_replicate_group_combos, ch_fai, stat_window, true)

    // Identify peaks of control window diff statistic using sliding window
    EXTRACTCTRLWINDOWEDPEAKS(CALCULATECTRLWINDOWDIFFS.out.tsv, peak_window)

    // Tabix all stat TSVs
    TABIX_MERGEDSTARTS(MERGEREADSTARTS.out.tsv)
    TABIX_WINDOWDIFFS(CALCULATEWINDOWDIFFS.out.tsv)
    TABIX_WINDOWEDPEAKS(EXTRACTWINDOWEDPEAKS.out.tsv)
    TABIX_WINDOWSUMS(CALCULATEWINDOWSUMS.out.tsv)
    TABIX_CTRLWINDOWSUMS(CALCULATECTRLWINDOWSUMS.out.tsv)
    TABIX_CTRLWIDEWINDOWSUMS(CALCULATECTRLWIDEWINDOWSUMS.out.tsv)
    ch_versions = ch_versions.mix(TABIX_MERGEDSTARTS.out.versions)

    // Gather inputs for site calling
    ch_test_window_diffs = CALCULATEWINDOWDIFFS.out.tsv
        .join(
            TABIX_WINDOWDIFFS.out.tbi,
            failOnDuplicate: true,
            failOnMismatch: true
        )
    ch_test_window_sums = CALCULATEWINDOWSUMS.out.tsv
        .join(
            TABIX_WINDOWSUMS.out.tbi,
            failOnDuplicate: true,
            failOnMismatch: true
        )
    ch_ctrl_window_sums = CALCULATECTRLWINDOWSUMS.out.tsv
        .join(
            TABIX_CTRLWINDOWSUMS.out.tbi,
            failOnDuplicate: true,
            failOnMismatch: true
        )
    ch_ctrl_wide_window_sums = CALCULATECTRLWIDEWINDOWSUMS.out.tsv
        .join(
            TABIX_CTRLWIDEWINDOWSUMS.out.tbi,
            failOnDuplicate: true,
            failOnMismatch: true
        )
    ch_site_calling_inputs = EXTRACTWINDOWEDPEAKS.out.tsv
        .join(
            ch_test_window_sums,
            failOnDuplicate: true,
            failOnMismatch: true
        )
        .map { meta, test_peaks_tsv, test_sums_tsv, test_sums_tbi ->
            [ meta.control_group, meta, test_peaks_tsv, test_sums_tsv, test_sums_tbi ]
        }
        .combine(
            ch_ctrl_window_sums.map { meta, tsv, tbi ->
                [ meta.sample_group, tsv, tbi ]
            },
            by: 0
        )
        .combine(
            ch_ctrl_wide_window_sums.map { meta, tsv, tbi ->
                [ meta.sample_group, tsv, tbi ]
            },
            by: 0
        )
        .combine(
            CALCULATECTRLWIDEWINDOWSUMS.out.dist.map { meta, tsv ->
                [ meta.sample_group, tsv ]
            },
            by: 0
        )
        .combine(
            EXTRACTCTRLWINDOWEDPEAKS.out.dist.map { meta, tsv ->
                [ meta.sample_group, tsv ]
            },
            by: 0
        )
        .map { vals -> vals[1..-1] }

    // Report/call sites based on test signal strength relative to control signal distribution
    // and background control distribution
    REPORTANDCALLSITES(
        ch_site_calling_inputs,
        test_pval_reporting_thresh,
        test_pval_calling_thresh,
        ctrl_noise_pval_calling_thresh,
        include_noise_pval_fails,
        ch_fasta,
        ch_fai,
        motif_search_factor,
        motif_match_score,
        motif_mismatch_pen,
        motif_gap_pen
    )

    // Create sample groups (groups of replicate groups, i.e. all samples within a sample group)
    // Figure out how many replicate groups there are per sample group (allows groups to emit as they're ready)
    ch_sample_group_sizes = getSampleGroupSizes(ch_meta)
    ch_test_sites = REPORTANDCALLSITES.out.site_tsv
        .join(
            ch_test_window_diffs,
            failOnDuplicate: true,
            failOnMismatch: true
        )
        .join(
            ch_test_window_sums,
            failOnDuplicate: true,
            failOnMismatch: true
        )
    ch_agg_site_inputs = createSampleGroups(ch_test_sites, ch_sample_group_sizes)
        .map { vals -> [ vals[0].control_group ] + vals }
        .combine(
            ch_ctrl_window_sums.map { meta, tsv, tbi ->
                [ meta.sample_group, tsv, tbi ]
            },
            by: 0
        )
        .combine(
            ch_ctrl_wide_window_sums.map { meta, tsv, tbi ->
                [ meta.sample_group, tsv, tbi ]
            },
            by: 0
        )
        .combine(
            CALCULATECTRLWIDEWINDOWSUMS.out.dist.map { meta, tsv ->
                [ meta.sample_group, tsv ]
            },
            by: 0
        )
        .combine(
            EXTRACTCTRLWINDOWEDPEAKS.out.dist.map { meta, tsv ->
                [ meta.sample_group, tsv ]
            },
            by: 0
        )
        .map { vals -> vals[1..-1] }

    // Aggregate reported/called sites across RNP concentrations
    AGGREGATESITES(
        ch_agg_site_inputs,
        test_pval_reporting_thresh,
        test_pval_calling_thresh,
        ctrl_noise_pval_calling_thresh,
        include_noise_pval_fails,
        agg_peak_window
    )

    emit:
    count_reads_json        = COUNTREADSTARTS.out.stats_json       // channel: [ val(meta), [ json   ] ]
    site_tsv                = REPORTANDCALLSITES.out.site_tsv      // channel: [ val(meta), [ tsv    ] ]
    agg_site_tsv            = AGGREGATESITES.out.agg_site_tsv      // channel: [ val(meta), [ tsv    ] ]
    window_diff_bigwig      = TSVTOBIGWIGDIFFS.out.bigwig          // channel: [ val(meta), [ bigwig ] ]
    peaks_dist_tsv          = EXTRACTWINDOWEDPEAKS.out.dist        // channel: [ val(meta), [ tsv    ] ]
    ctrl_peaks_dist_tsv     = EXTRACTCTRLWINDOWEDPEAKS.out.dist    // channel: [ val(meta), [ tsv    ] ]
    ctrl_wide_sums_dist_tsv = CALCULATECTRLWIDEWINDOWSUMS.out.dist // channel: [ val(meta), [ tsv    ] ]

    versions = ch_versions                                         // channel: [ versions.yml ]
}

