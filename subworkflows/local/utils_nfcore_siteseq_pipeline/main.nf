//
// Subworkflow with functionality specific to the Caribou-Biosciences/nf-SITE-Seq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                def meta_clone = [ id: "${meta.sample_group}_${meta.conc}nM_R${meta.rep}" ] + meta.clone()
                if (!fastq_2) {
                    if (params.use_umis) {
                        error(
                            "--use_umis is set but sample sheet row has no fastq_2"
                        )
                    }
                    return [ meta_clone.id, meta_clone + [ single_end:true, has_r2:false ], [ fastq_1 ] ]
                } else {
                    // Despite having 2 FASTQs, we set single_end:true because R2 is only used to
                    // extract UMIs, i.e. only R1 is used during genome alignment
                    // We use the has_r2 field to track whether R2 is present and overwrite single_end
                    // downstream to "trick" other tools into analyzing R1/R2 based on need
                    return [ meta_clone.id, meta_clone + [ single_end:true, has_r2:true ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { sample_rows ->
            validateReplicateSamples(sample_rows)
        }
        .map {
            meta, fastqs ->
                return [ meta: meta, fastqs: fastqs.flatten() ]
        }
        .toSortedList { a, b -> (b.meta.conc == 0) <=> (a.meta.conc == 0) } // Sort so that controls are first
        .flatten()
        .map { it -> [ it.meta, it.fastqs ] }
        .set { ch_samples }


    // Do a test run through of creating replicate groups (groups of samples with the same sample
    // group and RNP concentration), replicate group combinations (matched test and control groups),
    // and sample groups (groups of samples with the same sample group). Validation is baked into
    // these functions so the pipeline will fail before doing data processing if there are metadata
    // issues.
    ch_meta = ch_samples.map { meta, _fastqs -> meta }
    ch_rep_group_sizes = getReplicateGroupSizes(ch_meta)
    ch_control_group_limits = getControlReplicateGroupLimits(ch_meta)
    ch_replicate_groups = createReplicateGroups(ch_samples, ch_rep_group_sizes, ch_control_group_limits, false)
    ch_replicate_groups
    _ch_rep_group_combos = createReplicateGroupCombos(ch_replicate_groups)
    ch_sample_group_sizes = getSampleGroupSizes(ch_meta)
    _ch_sample_groups = createSampleGroups(ch_replicate_groups, ch_sample_group_sizes)


    emit:
    samplesheet = ch_samples
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_report.toList()
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//
// Return a string identifier for the given sample group and concentration
//
def getReplicateGroupId(sample_group, conc) {
    return "${sample_group}_${conc}nM"
}

//
// Return a channel containing each replicate group ID and how many samples
// are in that group
//
def getReplicateGroupSizes(ch_meta) {
    def ch_group_size = ch_meta
        .map { meta ->
            return [ getReplicateGroupId(meta.sample_group, meta.conc), meta ]
        }
        .groupTuple()
        .map { key, metas ->
            return [ key, metas.size() ]
        }
    return ch_group_size
}

//
// Return a channel containing each control sample group and how many samples
// are in that group
//
def getControlReplicateGroupSizes(ch_meta) {
    def ch_control_group_sizes = ch_meta
        .filter { meta -> meta.conc == 0 }
        .map { meta ->
            return [ meta.sample_group, meta ]
        }
        .groupTuple()
        .map { key, metas ->
            return [ key, metas.size() ]
        }
    return ch_control_group_sizes
}

//
// Return a channel containing each control sample group and how many samples are
// in each test replicate group that references that control group
// This is used to limit the number of control samples are included in test + control
// analyses
//
def getControlReplicateGroupLimits(ch_meta) {
    def ch_control_group_limit = ch_meta
        .filter { meta -> meta.conc > 0 }
        .map { meta ->
            return [ getReplicateGroupId(meta.sample_group, meta.conc), meta ]
        }
        .groupTuple()
        .map { _key, metas ->
            return [ metas[0].control_group, metas.size() ]
        }
        .groupTuple()
        .map { key, counts ->
            def unique_counts = counts.unique().collect();
            if (unique_counts.size() != 1) {
                error(
                    "Control group ${key} is referenced by test replicate groups (i.e. a test "+
                    "sample group at a particular concentration) with differing numbers of " +
                    "replicates (${unique_counts.join(', ')}). All test replicate groups that " +
                    "reference the same control group must have the same number of replicates"
                )
            }
            return [ key, counts[0] ]
        }
    // Add dummy test sample group records (allows use of .combine in createReplicateGroups)
    def ch_full_control_group_limit = ch_control_group_limit
        .concat(
            ch_meta.map { meta -> [ meta.sample_group ] }
            .groupTuple()
            .join(ch_control_group_limit, by: 0, remainder: true)
            .filter { _sample_group, matched -> !matched }
        )
    return ch_full_control_group_limit
}

//
// Return a channel of replicate groups from a channel of files - these are
// samples with the same sample group and RNP concentration
// Uses ch_group_sizes to know how many samples there are in each group - this
// allows groups to be emitted as soon as their samples are available
// Uses ch_control_group_limits to limit the number of control samples analyzed
//
def createReplicateGroups(ch_files, ch_group_sizes, ch_control_group_limits, limit_controls = true) {
    def replicate_groups = ch_files
        .map { meta, file ->
            def key = getReplicateGroupId(meta.sample_group, meta.conc)
            return [ meta.sample_group, key, meta, file ]
        }
        .combine(ch_control_group_limits, by: 0)
        .map { _sample_group, repGroupKey, meta, file, controlGroupLimit ->
            [ repGroupKey, meta, file, controlGroupLimit ]
        }
        .combine(ch_group_sizes, by: 0)
        .map { key, meta, file, controlGroupLimit, groupSize ->
            return tuple( groupKey(key, groupSize), [ meta, file, controlGroupLimit ] )
        }
        .groupTuple(by: 0, sort: { s1, s2 -> s1[0].rep <=> s2[0].rep })
        .map { groupId, tups ->
            def metas = tups.collect { meta, _file, _lim -> meta }
            validateReplicateGroup(metas)
            def ref_meta = metas[0]
            def reps = tups.collect { meta, _file, _lim -> meta.rep }
            def files = tups.collect { _meta, file, _lim -> file }
            if (limit_controls && ref_meta.conc == 0) {
                // Truncate control replicates to the number of test group replicates or 2, whichever is greater
                def controlGroupLimit = Math.max(tups[0][2], 2)
                def endIdx = controlGroupLimit - 1
                reps = reps[0..endIdx]
                files = files[0..endIdx]
            }
            def group_meta = [
                id: groupId,
                sample_group: ref_meta.sample_group,
                control_group: ref_meta.control_group,
                conc: ref_meta.conc,
                reps: reps,
                is_control: ref_meta.conc == 0,
                single_end: ref_meta.single_end,
                has_r2: ref_meta.has_r2,
                on_target_motif: ref_meta.on_target_motif,
                on_target_location: ref_meta.on_target_location
            ]
            return [ group_meta, files ]
        }
    return replicate_groups
}

//
// Get number of concentrations (i.e. replicate groups) per sample group
//
def getSampleGroupSizes(ch_meta) {
    def ch_group_size = ch_meta
        .map { meta ->
            return [ meta.sample_group, meta ]
        }
        .groupTuple()
        .map { key, metas ->
            return [ key, metas.collect{ meta -> meta.conc }.unique().size() ]
        }
    return ch_group_size
}

//
// Create a channel of sample groups from a channel of files
// Uses ch_group_sizes to know how many samples there are in each group - this
// allows groups to be emitted as soon as their samples are available
//
def createSampleGroups(ch_files, ch_group_sizes) {
    def sample_groups = ch_files
        .map { vals ->
            return tuple( vals[0].sample_group, vals[0], vals[1..-1] )
        }
        .combine(ch_group_sizes, by: 0)
        .map { key, meta, files, groupSize ->
            return tuple( groupKey(key, groupSize), [ meta, files ] )
        }
        .groupTuple(by: 0, sort: { s1, s2 -> s1[0].conc <=> s2[0].conc ?: s1[0].rep <=> s2[0].rep })
        .map { groupId, tups ->
            def metas = tups.collect { meta, _file -> meta }
            validateSampleGroup(metas)
            def ref_meta = metas[0]
            def concs = tups.collect { meta, _file -> meta.conc }
            def files = tups.collect { _meta, file -> file }
            def group_meta = [
                id: groupId,
                sample_group: ref_meta.sample_group,
                control_group: ref_meta.control_group,
                concs: concs,
                is_control: ref_meta.is_control,
                single_end: ref_meta.single_end,
                has_r2: ref_meta.has_r2,
                on_target_motif: ref_meta.on_target_motif,
                on_target_location: ref_meta.on_target_location
            ]
            return [ group_meta ] + files.transpose()
        }
    return sample_groups
}

//
// Create a channel of test replicate groups combined with their matched control replicate
// group
//
def createReplicateGroupCombos(ch_replicate_groups) {
    def ch_replicate_groups_branch = ch_replicate_groups
        .branch { meta, _file ->
            control: meta.is_control
            test: !meta.is_control
        }

    def ch_rep_group_combos = ch_replicate_groups_branch.test
        .map { meta, file ->
            return [ meta.control_group, meta, file ]
        }
        .combine(
            ch_replicate_groups_branch.control.map { meta, file ->
                return [ meta.sample_group, meta, file ]
            },
            by: 0
        )
        .map { _control_group, test_meta, test_file, control_meta, control_file ->
            return [ test_meta, test_file, control_meta, control_file ]
        }

    validateControlGroupReferences(ch_replicate_groups, ch_rep_group_combos)
    return ch_rep_group_combos
}

//
// Validate channels from input samplesheet
//
def validateReplicateSamples(sample_rows) {
    def (metas, fastqs) = sample_rows[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.has_r2 }.unique().size == 1
    if (!endedness_ok) {
        error(
            "Please check input samplesheet -> Multiple runs of a sample must be of the same " +
            "datatype i.e. single-end or paired-end (${metas[0].sample_group}, replicate ${metas[0].rep})"
        )
    }

    metas.collect { meta ->
        if (meta.conc < 0) {
            error("Please check input samplesheet -> Concentration must be ≥0")
        } else if (meta.conc == 0) {
            if (meta.on_target_motif) {
                error(
                    "Please check input samplesheet -> on_target_motif should be blank for " +
                    "control samples (${meta.sample_group}, replicate ${meta.rep})"
                )
            } else if (meta.on_target_location) {
                error(
                    "Please check input samplesheet -> on_target_location should be blank for " +
                    "control samples (${meta.sample_group}, replicate ${meta.rep})"
                )
            } else if (meta.control_group) {
                error(
                    "Please check input samplesheet -> control_group should be blank for " +
                    "control samples (${meta.sample_group}, replicate ${meta.rep})"
                )
            }
        } else {
            if (!meta.on_target_motif) {
                error(
                    "Please check input samplesheet -> on_target_motif must be defined for all " +
                    "test samples (${meta.sample_group}, conc. ${meta.conc}, replicate ${meta.rep})"
                )
            } else if (!meta.on_target_location) {
                error(
                    "Please check input samplesheet -> on_target_location must be defined for " +
                    "all test samples (${meta.sample_group}, conc. ${meta.conc}, replicate ${meta.rep})"
                )
            } else if (!meta.control_group) {
                error(
                    "Please check input samplesheet -> control_group must be defined for all test " +
                    "samples (${meta.sample_group}, conc. ${meta.conc}, replicate ${meta.rep})"
                )
            }
        }
    }

    return [ metas[0], fastqs ]
}

//
// Validate a collection of meta maps that have the same sample group and RNP concentration
//
def validateReplicateGroup(metas) {
    def sample_group = metas[0].sample_group
    def conc = metas[0].conc
    if (conc == 0) {
        if (metas.size < 2) {
            error(
                "Control groups must have at least 2 replicates (${sample_group})"
            )
        } else if (metas.collect{ meta -> meta.control_group }.any()) {
            log.warn(
                "The control_group column for control sample groups should be blank"
            )
        } else if (metas.collect{ meta -> meta.on_target_motif }.any()) {
            log.warn(
                "The on_target_motif column for control sample groups should be blank"
            )
        } else if (metas.collect{ meta -> meta.on_target_location }.any()) {
            log.warn(
                "The on_target_location column for control sample groups should be blank"
            )
        }

    } else {
        def control_groups_ok = metas.collect{ meta -> meta.control_group }.unique().size() == 1
        if (!control_groups_ok) {
            error(
                "Please check input samplesheet -> All samples in a sample group must have the " +
                "same control_group: ${sample_group}"
            )
        }

        def motifs_ok = metas.collect{ meta -> meta.on_target_motif }.unique().size == 1
        if (!motifs_ok) {
            error(
                "Please check input samplesheet -> All samples in a sample group must have the " +
                "same on_target_motif: ${sample_group}"
            )
        }

        def motif_locs_ok = metas.collect{ meta -> meta.on_target_location }.unique().size == 1
        if (!motif_locs_ok) {
            error(
                "Please check input samplesheet -> All samples in a sample group must have the " +
                "same on_target_location: ${sample_group}"
            )
        }
    }
}

//
// Validate control group references
//
def validateControlGroupReferences(replicate_groups, replicate_group_combos) {
    replicate_group_combos.map { test_meta, _test_file, control_meta, _control_file ->
        validateReplicateGroupCombo(test_meta, control_meta)
    }

    def replicate_groups_branch = replicate_groups
        .branch { meta, _file ->
            control: meta.is_control
            test: !meta.is_control
        }

    // Validate that all test groups have a control group
    def combo_test_groups = replicate_group_combos
        .map { test_meta, test_file, _control_meta, _control_file -> [ test_meta, test_file ] }
    replicate_groups_branch.test
        .join(combo_test_groups, by: [0], remainder: true)
        .map { meta, _file, match ->
            if (!match) {
                error(
                    "Could not identify any samples matching control_group ${meta.control_group} " +
                    "for sample group ${meta.sample_group}, concentration ${meta.conc}"
                )
            }
        }

    // Validate that all control groups are referenced by a test group
    def combo_control_groups = replicate_group_combos
        .map { _test_meta, _test_file, control_meta, control_file -> [ control_meta, control_file ]}
    replicate_groups_branch.control
        .join(combo_control_groups, by: [0], remainder: true)
        .map { meta, _file, match ->
            if (!match) {
                log.warn(
                    "Control sample group ${meta.sample_group} is not referenced by any test samples, " +
                    "these control samples will not be used in off-target analysis"
                )
            }
        }
}


def validateReplicateGroupCombo(test_group_meta, control_group_meta) {
    if (test_group_meta.control_group != control_group_meta.sample_group) {
        error(
            "Error during samplesheet parsing, control_group for sample group " +
            "${test_group_meta.sample_group} is ${test_group_meta.control_group} but given " +
            "control group ${control_group_meta.sample_group}"
        )
    } else if (test_group_meta.reps.size > control_group_meta.reps.size) {
        error(
            "Sample group ${test_group_meta.sample_group} has more replicates " +
            "(${test_group_meta.reps.size}) than control group " +
            "${control_group_meta.sample_group} (${control_group_meta.reps.size})"
        )
    } else if (test_group_meta.reps.size < control_group_meta.reps.size && test_group_meta.reps.size > 1) {
        log.warn(
            "Sample group ${test_group_meta.sample_group} has fewer replicates " +
            "(${test_group_meta.reps.size}) than control group " +
            "${control_group_meta.sample_group} (${control_group_meta.reps.size}), extra " +
            "control replicates won't be used when analyzing ${test_group_meta.sample_group}"
        )
    } else if (test_group_meta.has_r2 != control_group_meta.has_r2) {
        error(
            "Sample group ${test_group_meta.sample_group} datatype (i.e. single-end or " +
            "paired-end) does not match control group datatype"
        )
    }
}


def validateSampleGroup(metas) {
    def sample_group = metas[0].sample_groups

    def control_groups_ok = metas.collect{ meta -> meta.control_group }.unique().size == 1
    if (!control_groups_ok) {
        error(
            "Please check input samplesheet -> Sample group ${sample_group} references multiple " +
            "control groups"
        )
    }

    def motifs_ok = metas.collect{ meta -> meta.on_target_motif }.unique().size == 1
    if (!motifs_ok) {
        error(
            "Please check input samplesheet -> ${sample_group} has more than 1 unique on_target_motif"
        )
    }

    def motif_locs_ok = metas.collect{ meta -> meta.on_target_location }.unique().size == 1
    if (!motif_locs_ok) {
        error(
            "Please check input samplesheet -> ${sample_group} has more than 1 unique on_target_location"
        )
    }

    def conc_counts = [:]
    metas.collect { meta ->
        meta.reps.collect { _rep ->
            if (!conc_counts[meta.conc]) {
                conc_counts[meta.conc] = 1
            } else {
                conc_counts[meta.conc] += 1
            }
        }
    }

    if (conc_counts.values().flatten().unique().size != 1) {
        def count_str = conc_counts.collect { conc, count ->
            "${conc}: ${count} sample(s)"
        }.join(", ")
        log.warn(
            "Concentration groups for sample group ${sample_group} do not all have the same " +
            "number of samples (${count_str})"
        )
    }
}


//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
            "Tools used in the workflow included:",
            "FastP (Chen 2023),",
            params["use_umis"] ? "UMI-tools (Smith et al. 2017)" : "",
            "Bowtie2 (Langmead et al. 2012)",
            "MultiQC (Ewels et al. 2016)",
            "attrs (Schlawack)",
            "pandas",
            "NumPy (Harris et al. 2020)",
            "Parasail (Daily 2016)",
            "plotly.py (Kruchten et al. 2024)",
            "pysam",
            "samtools (Danecek et al. 2021)",
            "HTSLib (Bonfield et al. 2021)",
            "nlohmann/json (Lohmann 2023)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            "<li>Shifu Chen. 2023. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta 2: e107. https://doi.org/10.1002/imt2.107.</li>",
            params["use_umis"] ? "<li>Smith T, Heger A, Sudbery I. UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Res. 2017 Mar;27(3):491-499. doi: 10.1101/gr.209601.116. Epub 2017 Jan 18. PMID: 28100584; PMCID: PMC5340976.</li>" : "",
            "<li>Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>",
            "<li>Schlawack, H. attrs [Computer software]. https://doi.org/10.5281/zenodo.6925130</li>",
            "<li>The pandas development team. pandas-dev/pandas: Pandas [Computer software]. https://doi.org/10.5281/zenodo.3509134</li>",
            "<li>Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2.</li>",
            "<li>Daily, Jeff. (2016). Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments. BMC Bioinformatics, 17(1), 1-11. doi:10.1186/s12859-016-0930-z</li>",
            "<li>Kruchten, N., Seier, A., & Parmer, C. (2024). An interactive, open-source, and browser-based graphing library for Python (Version 5.24.1) [Computer software]. https://doi.org/10.5281/zenodo.14503524</li>",
            "<li>Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li, Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008</li>",
            "<li>James K Bonfield, John Marshall, Petr Danecek, Heng Li, Valeriu Ohan, Andrew Whitwham, Thomas Keane, Robert M Davies, HTSlib: C library for reading/writing high-throughput sequencing data, GigaScience, Volume 10, Issue 2, February 2021, giab007, https://doi.org/10.1093/gigascience/giab007</li>",
            "<li>Lohmann, N. (2023). JSON for Modern C++ (Version 3.11.3) [Computer software]. https://github.com/nlohmann</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

