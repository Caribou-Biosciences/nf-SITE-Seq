
process COUNTREADSTARTS {
    tag "$meta.id"
    label 'process_single'
    label 'site_seq_toolkit'

    input:
    tuple val(meta), path(bam)
    val(apply_normalization)
    val(min_aln_length)
    val(min_map_qual)
    val(max_clipping)

    output:
    tuple val(meta), path("*.tsv.gz"), emit: tsv
    tuple val(meta), path("*.json")  , emit: stats_json

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def apply_norm = apply_normalization ? "--apply-norm" : ""
    """
    count_read_starts \\
        --bam-path $bam \\
        --stats-path ${prefix}.stats.json \\
        --min-aln-length $min_aln_length \\
        --min-map-qual $min_map_qual \\
        --max-clipping $max_clipping \\
        $apply_norm \\
        | bgzip -c > ${prefix}.read_starts.tsv.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.read_starts.tsv.gz
    touch ${prefix}.stats.json
    """
}
