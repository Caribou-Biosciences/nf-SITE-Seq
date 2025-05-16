
process CALCULATEWINDOWDIFFS {
    tag "$test_meta.id"
    label 'process_single'
    label 'site_seq_toolkit'

    input:
    tuple val(test_meta), path(test_merged_tsv), val(control_meta), path(control_merged_tsv)
    tuple val(genome_meta), path(fai)
    val(stat_window)
    val(self_comparison)

    output:
    tuple val(test_meta), path("*.tsv.gz"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${test_meta.id}"
    def test = "-t $test_merged_tsv"
    def control = control_meta ? "-c $control_merged_tsv" : ""
    def self_comp = self_comparison ? "-s" : ""
    """
    calculate_window_diffs \\
        $test \\
        $control \\
        -w $stat_window \\
        -f $fai \\
        $self_comp \\
        | bgzip -c > ${prefix}.window_diffs.W${stat_window}.tsv.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${test_meta.id}"
    """
    touch ${prefix}.window_diffs.W${stat_window}.tsv.gz
    """
}
