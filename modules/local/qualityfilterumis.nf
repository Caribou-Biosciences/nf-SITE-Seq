
process QUALITYFILTERUMIS {
    tag "$meta.id"
    label 'process_single'
    label 'site_seq_toolkit'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.json")    , emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    quality_filter_umis \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o ${prefix}.umi-filtered-1.fastq.gz \\
        -O ${prefix}.umi-filtered-2.fastq.gz \\
        -r ${prefix}.umi-filter-report.json \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi-filtered-1.fastq.gz
    touch ${prefix}.umi-filtered-2.fastq.gz
    touch ${prefix}.umi-filter-report.json
    """
}
