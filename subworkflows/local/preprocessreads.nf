include { CAT_FASTQ         } from '../../modules/nf-core/cat/fastq/main'
include { FASTP             } from '../../modules/nf-core/fastp/main'
include { QUALITYFILTERUMIS } from '../../modules/local/qualityfilterumis'
include { UMITOOLS_EXTRACT  } from '../../modules/nf-core/umitools/extract/main'
include { SEQTK_SAMPLE      } from '../../modules/nf-core/seqtk/sample/main'

workflow PREPROCESSREADS {

    take:
    ch_fastq                    // channel: [ val(meta), [ fastq(s) ] ]
    subsample_fastqs            // val: whether to subsample input FASTQs
    subsample_target_reads      // val: target number of reads to subsample to

    main:

    ch_versions = Channel.empty()

    // Hack: populate single_end field with custom has_r2 field so that read processing tools
    // will include R2 when it is present
    ch_fastq = ch_fastq.map { meta, fastqs ->
        def meta_clone = meta.clone()
        meta_clone.single_end = !meta_clone.has_r2
        return [ meta_clone, fastqs ]
    }

    if (subsample_fastqs) {
        // subsample FASTQs if requested
        ch_fastq_subsample = ch_fastq
            .map{[ it[0], it[1], subsample_target_reads ]}
        SEQTK_SAMPLE(
            ch_fastq_subsample
        )
        ch_fastq = SEQTK_SAMPLE.out.reads
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)
    }

    // Concatenate FASTQs for samples with multiple FASTQS (e.g. multiple lanes, multiple flow cells)
    ch_fastq
        .groupTuple(by: [0])
        .branch { meta, fastqs ->
            single: fastqs.size() == 1
                return [ meta, fastqs.flatten() ]
            multiple: fastqs.size() > 1
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_fastq_branch }
    CAT_FASTQ(
        ch_fastq_branch.multiple
    )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    ch_cat_fastq = CAT_FASTQ.out.reads
        .mix(ch_fastq_branch.single)

    // Run Fastp for quality/length filtering
    FASTP(
        ch_cat_fastq, [], false, false, false
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    if (params.use_umis) {
        // Apply custom UMI quality filtering that allows a variable number of bases below threshold
        QUALITYFILTERUMIS(
            FASTP.out.reads
        )
        ch_umi_filter_report = QUALITYFILTERUMIS.out.report

        // Extract UMIs using UMI tools
        // UMI tools v1.1.5 (latest version on biocontainers) does not directly support read2-only UMIs
        // so read1/read2 must be swapped.
        // TODO: UMI tools v1.1.6 does support read2-only UMIs, update accordingly when container is available
        ch_reads_rev = QUALITYFILTERUMIS.out.reads.map { meta, reads -> [ meta, reads.reverse() ] }
        UMITOOLS_EXTRACT(
            ch_reads_rev
        )
        ch_reads = UMITOOLS_EXTRACT.out.reads.map { meta, reads -> [ meta, reads.reverse() ] }
        ch_umitools_extract_log = UMITOOLS_EXTRACT.out.log
    } else {
        ch_reads = FASTP.out.reads
        ch_umi_filter_report = ch_fastq.map { meta, _fastqs -> [ meta, [] ]}
        ch_umitools_extract_log = ch_fastq.map { meta, _fastqs -> [ meta, [] ]}
    }

    // Hack: remove R2 and set single_end=true so downstream tools will only use R1
    ch_reads = ch_reads.map { meta, reads ->
        def meta_clone = meta.clone()
        meta_clone.single_end = true
        return [ meta_clone, reads[0] ]
    }

    emit:
    reads                = ch_reads                // channel: [ val(meta), [ fastq(s) ] ]
    fastp_json           = FASTP.out.json          // channel: [ val(meta), [ json     ] ]
    umi_filter_report    = ch_umi_filter_report    // channel: [ val(meta), [ json     ] ]
    umitools_extract_log = ch_umitools_extract_log // channel: [ val(meta), [ log      ] ]

    versions   = ch_versions                       // channel: [ versions.yml            ]
}

