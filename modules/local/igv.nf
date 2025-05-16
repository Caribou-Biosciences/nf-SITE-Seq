process IGV {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3':
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), val(panel_paths)
    tuple val(fa_meta), val(fasta_path)

    output:
    tuple val(meta), path(igv_xml), emit: xml
    path "versions.yml"           , emit: versions

    script:
    igv_xml = "${meta.id}_igv_session.xml"
    def panel_args = panel_paths.collect { paths ->
        "--panel ${paths.flatten().join(' ')}"
    }
    """
    create_igv_session.py \\
        --out-path $igv_xml \\
        --genome-ref $fasta_path \\
        ${panel_args.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    igv_xml = "${meta.id}_igv_session.xml"
    """
    touch $igv_xml
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
