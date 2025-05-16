
process CALCULATEWINDOWSUMS {
    tag "$meta.id"
    label 'process_single'
    label 'site_seq_toolkit'

    input:
    tuple val(meta), path(merged_tsv)
    tuple val(genome_meta), path(fai)
    val(window)

    output:
    tuple val(meta), path("*.tsv.gz"), emit: tsv
    tuple val(meta), path("*.tsv")   , emit: dist

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    calculate_window_sums \\
        -i $merged_tsv \\
        -f $fai \\
        -d ${prefix}.window_sums.W${window}.dist.tsv \\
        -w $window | bgzip -c > ${prefix}.window_sums.W${window}.tsv.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.window_sums.W${window}.tsv.gz
    touch ${prefix}.window_sums.W${window}.dist.tsv
    """
}
