include { UNTAR as UNTAR_INDEX   } from '../../modules/nf-core/untar/main'
include { UNZIP as UNZIP_INDEX   } from '../../modules/nf-core/unzip/main'
include { GUNZIP as GUNZIP_FASTA } from '../../modules/nf-core/gunzip/main'
include { BOWTIE2_BUILD          } from '../../modules/nf-core/bowtie2/build/main'
include { SAMTOOLS_FAIDX         } from '../../modules/nf-core/samtools/faidx/main'


process PUBLISHGENOMEFASTA {
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta)

    script:
    """
    echo "Publishing file: $fasta"
    """
}


workflow PREPAREGENOME {

    take:
    fasta              // file: /path/to/genome.fasta
    bowtie2_index      // file: /path/to/genome-bt2idx/

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    ch_fasta = Channel.empty()
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA([ [:], fasta ]).gunzip
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = [ [:], file(fasta, checkIfExists: true) ]
    }

    //
    // Uncompress Bowtie2 index or generate from scratch if required
    //
    if (bowtie2_index) {
        if (bowtie2_index.endsWith('.tar.gz')) {
            ch_bowtie2_index = UNTAR_INDEX([ [:], bowtie2_index ] ).untar
            ch_versions  = ch_versions.mix(UNTAR_INDEX.out.versions)
        } else if (bowtie2_index.endsWith('.zip')) {
            ch_bowtie2_index = UNZIP_INDEX([ [:], bowtie2_index ]).unzipped_archive
            ch_versions  = ch_versions.mix(UNZIP_INDEX.out.versions)
        } else {
            ch_bowtie2_index = [ [:], file(bowtie2_index, checkIfExists: true) ]
        }
    } else {
        ch_bowtie2_index = BOWTIE2_BUILD(ch_fasta).index
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    }

    SAMTOOLS_FAIDX(ch_fasta, [[], []])
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // Publish input FASTA to output folder
    PUBLISHGENOMEFASTA(ch_fasta)

    emit:
    fasta         = ch_fasta                    // channel: [ val(meta), [ fasta         ] ]
    fai           = SAMTOOLS_FAIDX.out.fai      // channel: [ val(meta), [ fai           ] ]
    bowtie2_index = ch_bowtie2_index            // channel: [ val(meta), [ bowtie2_index ] ]

    versions      = ch_versions                 // channel: [ versions.yml                 ]
}

