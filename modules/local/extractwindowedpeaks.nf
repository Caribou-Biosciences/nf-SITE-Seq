
process EXTRACTWINDOWEDPEAKS {
    tag "$meta.id"
    label 'process_single'
    label 'site_seq_toolkit'

    input:
    tuple val(meta), path(window_diffs_tsv)
    val(peak_window)

    output:
    tuple val(meta), path("*.tsv.gz"), emit: tsv
    tuple val(meta), path("*.dist.tsv"), emit: dist

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_windowed_peaks \\
        -i $window_diffs_tsv \\
        -d ${prefix}.peaks.dist.tsv \\
        -w $peak_window | bgzip -c > ${prefix}.peaks.tsv.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.peaks.tsv.gz
    touch ${prefix}.peaks.dist.tsv
    """
}
