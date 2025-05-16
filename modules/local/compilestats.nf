
process COMPILESTATS {
    tag "$meta.id"
    label 'process_single'
    label 'site_seq_toolkit'

    input:
    tuple(
        val(meta),
        path(fastp_json),
        path(umi_filter_report),
        path(bowtie2_log),
        path(umitools_dedup_log),
        path(count_reads_json)
    )

    output:
    tuple val(meta), path("*.tsv"), emit: stats_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def umi_filter_report_arg = umi_filter_report ? "--umi_filter_report $umi_filter_report" : ""
    def umitools_dedup_arg = umitools_dedup_log ? "--umitools_dedup_log $umitools_dedup_log" : ""

    """
    compile_stats.py \\
        --identifier ${meta.id} \\
        --sample_group ${meta.sample_group} \\
        --rep ${meta.rep} \\
        --conc ${meta.conc} \\
        --fastp_json $fastp_json \\
        $umi_filter_report_arg \\
        --bowtie2_log $bowtie2_log \\
        $umitools_dedup_arg \\
        --count_reads_json $count_reads_json \\
        --out_path ${prefix}_stats.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_stats.tsv
    """
}
