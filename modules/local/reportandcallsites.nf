
process REPORTANDCALLSITES {
    tag "$meta.id"
    label 'process_single'
    label 'site_seq_toolkit'

    input:
    tuple(
        val(meta),
        path(test_peaks_tsv),
        path(test_window_sums_tsv),
        path(test_window_sums_tbi),
        path(control_window_sums_tsv),
        path(control_window_sums_tbi),
        path(control_wide_window_sums_tsv),
        path(control_wide_window_sums_tbi),
        path(ctrl_wide_window_sum_dist_tsv, stageAs: "wide_window_dist.tsv"),
        path(ctrl_peaks_dist_tsv, stageAs: "ctrl_peaks_dist.tsv")
    )
    val(test_pval_reporting_thresh)
    val(test_pval_calling_thresh)
    val(control_pval_calling_thresh)
    val(include_control_pval_fails)
    tuple val(fa_meta), path(fasta)
    tuple val(fai_meta), path(fai)
    val(search_factor)
    val(motif_match_score)
    val(motif_mismatch_pen)
    val(motif_gap_pen)

    output:
    tuple val(meta), path("*.tsv"), emit: site_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def include_control_fails = include_control_pval_fails ? "--include_control_pval_fails" : ""
    """
    report_and_call_sites.py \\
        --identifier ${meta.id} \\
        --conc ${meta.conc} \\
        --output_tsv ${prefix}.reported_sites.tsv \\
        --test_peaks_tsv $test_peaks_tsv \\
        --test_window_sums_tsv $test_window_sums_tsv \\
        --control_window_sums_tsv $control_window_sums_tsv \\
        --control_wide_window_sums_tsv $control_wide_window_sums_tsv \\
        --ctrl_peaks_dist_tsv $ctrl_peaks_dist_tsv \\
        --ctrl_wide_window_sum_dist_tsv $ctrl_wide_window_sum_dist_tsv \\
        --test_pval_reporting_thresh $test_pval_reporting_thresh \\
        --test_pval_calling_thresh $test_pval_calling_thresh \\
        --control_pval_calling_thresh $control_pval_calling_thresh \\
        $include_control_fails \\
        --fasta_path $fasta \\
        --search_motif ${meta.on_target_motif} \\
        --on_target_location ${meta.on_target_location} \\
        --search_factor $search_factor \\
        --motif_match_score $motif_match_score \\
        --motif_mismatch_pen $motif_mismatch_pen \\
        --motif_gap_pen $motif_gap_pen
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reported_sites.tsv
    """
}
