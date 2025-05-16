process MERGEREADSTARTS {
    tag "$meta.id"
    label 'process_single'
    label 'site_seq_toolkit'

    input:
    tuple val(meta), path(tsvs)

    output:
    tuple val(meta), path("*.tsv.gz"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_read_starts \\
        --bam-paths ${tsvs.join(' ')} \\
        | bgzip -c > ${prefix}.tsv.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv.gz
    """
}
