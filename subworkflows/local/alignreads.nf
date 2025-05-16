include { BOWTIE2_ALIGN                          } from '../../modules/nf-core/bowtie2/align/main'
include { UMITOOLS_DEDUP                         } from '../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX                         } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILT  } from '../../modules/nf-core/samtools/index/main'

workflow ALIGNREADS {

    take:
    ch_reads              // channel: [ val(meta), [ fastq ] ]
    ch_fasta              // channel: [ val(meta), [ fasta ] ]
    ch_bowtie2_index      // channel: [ val(meta), [ /path/to/bowtie2/index/ ] ]
    use_umis              // value: whether to use UMIs to deduplicate alignments

    main:

    ch_versions = Channel.empty()

    //
    // Align reads to genome
    //
    BOWTIE2_ALIGN(
        ch_reads,
        ch_bowtie2_index,
        ch_fasta,
        false,
        true
    )
    ch_bam = BOWTIE2_ALIGN.out.bam
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

    SAMTOOLS_INDEX(
        BOWTIE2_ALIGN.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    if (use_umis) {
        UMITOOLS_DEDUP(
            BOWTIE2_ALIGN.out.bam
                .join(SAMTOOLS_INDEX.out.bai, by: [0]),
            true
        )
        ch_bam = UMITOOLS_DEDUP.out.bam
        ch_umitools_dedup_log = UMITOOLS_DEDUP.out.log
        ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions)

        SAMTOOLS_INDEX_DEDUP(
            ch_bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_DEDUP.out.versions)
    } else {
        ch_umitools_dedup_log = ch_bam.map { meta, _bam -> [ meta, [] ]}
    }

    emit:
    bam                = ch_bam                   // channel: [ val(meta), [ bam ] ]
    bowtie2_log        = BOWTIE2_ALIGN.out.log    // channel: [ val(meta), [ log ] ]
    umitools_dedup_log = ch_umitools_dedup_log    // channel: [ val(meta), [ log ] ]

    versions           = ch_versions              // channel: [ versions.yml ]
}

