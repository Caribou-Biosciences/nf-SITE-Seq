process TSVTOBIGWIG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:472--h9b8f530_1':
        'quay.io/biocontainers/ucsc-bedgraphtobigwig:472--h9b8f530_1' }"

    input:
    tuple val(meta), path(tsv)
    tuple val(genome_meta), path(fai)

    output:
    tuple val(meta), path("*.bw"), emit: bigwig

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    replicates=(${meta.reps.join(" ")})
    if [[ "$tsv" =~ \\.gz\$ ]]; then
        read_cmd="zcat"
    else
        read_cmd="cat"
    fi

    num_columns=\$(awk -F'\\t' '{print NF; exit}' <( \$read_cmd "$tsv" | grep -v '^#' ))
    if [[ -z "\$num_columns" || "\$num_columns" -le 3 ]]; then
        echo "Error: Not enough columns in '$tsv'"
        exit 1
    fi

    \$read_cmd "$tsv" | grep -v '^#' | awk -F'\\t' -v OFS='\\t' -v total_cols="\$num_columns" '
    BEGIN {
        for (i = 0; i < total_cols-3; i++) {
            file = "tmp_column_" i ".tsv"
            file_handles[i] = file
            print "Generating " file
        }
    }
    {
        for (col = 4; col <= total_cols; col++) {
            if (\$col > 0) {
                print \$1, \$2, \$3, \$col > file_handles[col-4]
            }
        }
    }
    END {
        # Close all output file handles
        for (col in file_handles) {
            close(file_handles[col])
        }
    }
    '

    end=\$((\$num_columns - 4))
    for idx in \$(seq 0 \$end); do
        rep=\${replicates[idx]}
        bedGraphToBigWig tmp_column_\$idx.tsv $fai ${prefix}_R\$rep.bw
    done
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    replicates=(${meta.reps.join(" ")})
    for rep in "\${replicates[@]}"; do
        touch ${prefix}_R\$rep.bw
    done
    """
}
