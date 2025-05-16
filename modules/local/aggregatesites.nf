
process AGGREGATESITES {
    tag "$meta.id"
    label 'process_single'
    label 'site_seq_toolkit'

    input:
    tuple(
        val(meta),
        path(called_site_tsvs),
        path(test_window_diffs_tsvs),
        path(test_window_diffs_tbis),
        path(test_window_sums_tsvs),
        path(test_window_sums_tbis),
        path(control_window_sums_tsv),
        path(control_window_sums_tbi),
        path(control_wide_window_sums_tsv),
        path(control_wide_window_sums_tbi),
        path(ctrl_wide_window_sum_dist_tsv),
        path(ctrl_peaks_dist_tsv)
    )
    val(test_pval_reporting_thresh)
    val(test_pval_calling_thresh)
    val(control_pval_calling_thresh)
    val(include_control_pval_fails)
    val(collapse_window)

    output:
    tuple val(meta), path("*.tsv"), emit: agg_site_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def include_control_fails = include_control_pval_fails ? "--include_control_pval_fails" : ""
    """
    aggregate_sites.py \\
        --identifier ${meta.id} \\
        --output_tsv ${prefix}.aggregated_sites.tsv \\
        --reported_site_tsvs ${called_site_tsvs.join(' ')} \\
        --test_window_diffs_tsvs ${test_window_diffs_tsvs.join(' ')} \\
        --test_window_sums_tsvs ${test_window_sums_tsvs.join(' ')} \\
        --control_window_sums_tsv $control_window_sums_tsv \\
        --control_wide_window_sums_tsv $control_wide_window_sums_tsv \\
        --ctrl_peaks_dist_tsv $ctrl_peaks_dist_tsv \\
        --ctrl_wide_window_sum_dist_tsv $ctrl_wide_window_sum_dist_tsv \\
        --test_pval_reporting_thresh $test_pval_reporting_thresh \\
        --test_pval_calling_thresh $test_pval_calling_thresh \\
        --control_pval_calling_thresh $control_pval_calling_thresh \\
        --site_collapse_window $collapse_window \\
        $include_control_fails
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aggregated_sites.tsv
    """
}
