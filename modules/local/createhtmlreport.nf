
process CREATEHTMLREPORT {
    tag "$meta.id"
    label 'process_single'
    label 'site_seq_toolkit'

    input:
    tuple(
        val(meta),
        path(agg_site_tsv),
        path(peaks_dist_tsvs),
        path(ctrl_peaks_dist_tsv),
        path(ctrl_wide_window_sum_dist_tsv)
    )
    path(read_stats_tsv)
    val(test_pval_reporting_thresh)
    val(test_pval_calling_thresh)
    val(background_pval_calling_thresh)

    output:
    tuple val(meta), path("*.html"), emit: html

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    create_html_report.py \\
        --sample_group ${meta.sample_group} \\
        --control_group ${meta.control_group} \\
        --on_target_motif ${meta.on_target_motif} \\
        --on_target_location ${meta.on_target_location} \\
        --read_stats_tsv $read_stats_tsv \\
        --agg_site_tsv $agg_site_tsv \\
        --peaks_dist_tsvs ${peaks_dist_tsvs.join(' ')} \\
        --ctrl_peaks_dist_tsv $ctrl_peaks_dist_tsv \\
        --ctrl_wide_window_sum_dist_tsv $ctrl_wide_window_sum_dist_tsv \\
        --test_pval_reporting_thresh $test_pval_reporting_thresh \\
        --test_pval_calling_thresh $test_pval_calling_thresh \\
        --background_pval_calling_thresh $background_pval_calling_thresh \\
        --out_path ${prefix}_report.html \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_report.html
    """
}
