# Caribou-Biosciences/nf-SITE-Seq: Usage

## Introduction

This document describes the requirements and options available to run the SITE-Seq® assay pipeline.

## Sample requirements

To run the SITE-Seq® assay analysis pipeline, at least one test sample (i.e. RNP concentration > 0 nM) and two control samples (i.e. RNP concentration = 0 nM) must be included in the samplesheet. If multiple test sample replicates are included, at least that many control sample replicates must also be included. For example, if 2 test samples are included, at least 2 control samples must be included. **We recommend running 3 test sample replicates per guide and RNP concentration and 3 control sample replicates.**

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyze before running the pipeline. Use the `--input` parameter to specify its location. See the examples below.

```bash
--input '[path to samplesheet file]'
```

### Full samplesheet

> [!WARNING]
> While the pipeline supports single-end sequencing data, this is incompatible with UMI-based analysis which must be disabled using `--use_umis false`

A full samplesheet testing two targets, each at two RNP concentrations, is given below. Single-end sequencing data may be used by omitting the `fastq_2` column, however, this is incompatible with UMI-based analysis which must be disabled using `--use_umis false`.

```csv title="samplesheet.csv"
sample_group,control_group,concentration,replicate,on_target_motif,on_target_location,fastq_1,fastq_2
CONTROL,,0,1,,,/data/CONTROL_Rep1_R1.fastq.gz,/data/CONTROL_Rep1_R2.fastq.gz
CONTROL,,0,2,,,/data/CONTROL_Rep2_R1.fastq.gz,/data/CONTROL_Rep2_R2.fastq.gz
CONTROL,,0,3,,,/data/CONTROL_Rep3_R1.fastq.gz,/data/CONTROL_Rep3_R2.fastq.gz
AAVS1,CONTROL,16,1,GGGGCCACTAGGGACAGGATNGG,chr19:55115749-55115771,/data/AAVS1_16_Rep1_R1.fastq.gz,/data/AAVS1_16_Rep1_R2.fastq.gz
AAVS1,CONTROL,16,2,GGGGCCACTAGGGACAGGATNGG,chr19:55115749-55115771,/data/AAVS1_16_Rep2_R1.fastq.gz,/data/AAVS1_16_Rep2_R2.fastq.gz
AAVS1,CONTROL,16,3,GGGGCCACTAGGGACAGGATNGG,chr19:55115749-55115771,/data/AAVS1_16_Rep3_R1.fastq.gz,/data/AAVS1_16_Rep3_R2.fastq.gz
AAVS1,CONTROL,128,1,GGGGCCACTAGGGACAGGATNGG,chr19:55115749-55115771,/data/AAVS1_128_Rep1_R1.fastq.gz,/data/AAVS1_128_Rep1_R2.fastq.gz
AAVS1,CONTROL,128,2,GGGGCCACTAGGGACAGGATNGG,chr19:55115749-55115771,/data/AAVS1_128_Rep2_R1.fastq.gz,/data/AAVS1_128_Rep2_R2.fastq.gz
AAVS1,CONTROL,128,3,GGGGCCACTAGGGACAGGATNGG,chr19:55115749-55115771,/data/AAVS1_128_Rep3_R1.fastq.gz,/data/AAVS1_128_Rep3_R2.fastq.gz
FANCF,CONTROL,16,1,GGAATCCCTTCTGCAGCACCNGG,chr11:22625786-22625808,/data/FANCF_16_Rep1_R1.fastq.gz,/data/FANCF_16_Rep1_R2.fastq.gz
FANCF,CONTROL,16,2,GGAATCCCTTCTGCAGCACCNGG,chr11:22625786-22625808,/data/FANCF_16_Rep2_R1.fastq.gz,/data/FANCF_16_Rep2_R2.fastq.gz
FANCF,CONTROL,16,3,GGAATCCCTTCTGCAGCACCNGG,chr11:22625786-22625808,/data/FANCF_16_Rep3_R1.fastq.gz,/data/FANCF_16_Rep3_R2.fastq.gz
FANCF,CONTROL,128,1,GGAATCCCTTCTGCAGCACCNGG,chr11:22625786-22625808,/data/FANCF_128_Rep1_R1.fastq.gz,/data/FANCF_128_Rep1_R2.fastq.gz
FANCF,CONTROL,128,2,GGAATCCCTTCTGCAGCACCNGG,chr11:22625786-22625808,/data/FANCF_128_Rep2_R1.fastq.gz,/data/FANCF_128_Rep2_R2.fastq.gz
FANCF,CONTROL,128,3,GGAATCCCTTCTGCAGCACCNGG,chr11:22625786-22625808,/data/FANCF_128_Rep3_R1.fastq.gz,/data/FANCF_128_Rep3_R2.fastq.gz
```

| Column               | Required | Description                                                                                                                                                                                                          |
| -------------------- | -------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample_group`       | Yes      | An identifier used to group samples. Typically this is the name of a particular CRISPR guide but may include other information on the conditions/treatment of the sample, e.g. a modification to the CRISPR protein. |
| `control_group`      | Yes      | The `sample_group` of the corresponding control samples. For control samples, this column should be empty. The same `control_group` may be used for multiple test sample groups.                                     |
| `concentration`      | Yes      | The RNP concentration used in the sample. For control samples, this column should be 0.                                                                                                                              |
| `replicate`          | Yes      | A unique numeric identifier that differentiates samples with the same `sample_group` and `concentration`.                                                                                                            |
| `on_target_motif`    | Yes      | The on-target motif sequence for the target; includes the spacer sequence and the PAM sequence. Degenerate PAM bases should be included, e.g. "NGG" for Cas9. For control samples, this column should be empty.      |
| `on_target_location` | Yes      | The genomic coordinates of the `on_target_motif` sequence in the reference genome used. For control samples, this column should be empty.                                                                            |
| `fastq_1`            | Yes      | An absolute path to the read 1 FASTQ file for this sample.                                                                                                                                                           |
| `fastq_2`            | No       | An absolute path to the read 2 FASTQ file for this sample. Optional, but UMI analysis must be disabled (`--use_umis false`) if this column is omitted or empty.                                                      |

An [example samplesheet](../assets/test_data/samplesheet.csv) has been provided with the pipeline.

### Multiple runs of the same sample

When a sample has been sequenced more than once, e.g. to increase sequencing depth, it can be included in multiple rows of the samplesheet as long as the `sample_group`, `concentration`, and `replicate` columns are identical—these columns are used to identify the sample. The pipeline will then concatenate the raw reads before performing any downstream analysis. Below is an example in which each sample has been sequenced across 2 lanes:

```csv title="samplesheet.csv"
sample_group,control_group,concentration,replicate,on_target_motif,on_target_location,fastq_1,fastq_2
CONTROL,,0,1,,,/data/CONTROL_Rep1_R1_L001_001.fastq.gz,/data/CONTROL_Rep1_R2_L001_001.fastq.gz
CONTROL,,0,1,,,/data/CONTROL_Rep1_R1_L002_001.fastq.gz,/data/CONTROL_Rep1_R2_L002_001.fastq.gz
CONTROL,,0,2,,,/data/CONTROL_Rep2_R1_L001_001.fastq.gz,/data/CONTROL_Rep2_R2_L001_001.fastq.gz
CONTROL,,0,2,,,/data/CONTROL_Rep2_R1_L002_001.fastq.gz,/data/CONTROL_Rep2_R2_L002_001.fastq.gz
CONTROL,,0,3,,,/data/CONTROL_Rep3_R1_L001_001.fastq.gz,/data/CONTROL_Rep3_R2_L001_001.fastq.gz
CONTROL,,0,3,,,/data/CONTROL_Rep3_R1_L002_001.fastq.gz,/data/CONTROL_Rep3_R2_L002_001.fastq.gz
AAVS1,CONTROL,16,1,GGGGCCACTAGGGACAGGATNGG,chr19:55115749-55115771,/data/AAVS1_16_Rep1_R1_L001_001.fastq.gz,/data/AAVS1_16_Rep1_R2_L001_001.fastq.gz
AAVS1,CONTROL,16,1,GGGGCCACTAGGGACAGGATNGG,chr19:55115749-55115771,/data/AAVS1_16_Rep1_R1_L002_001.fastq.gz,/data/AAVS1_16_Rep1_R2_L002_001.fastq.gz
AAVS1,CONTROL,16,2,GGGGCCACTAGGGACAGGATNGG,chr19:55115749-55115771,/data/AAVS1_16_Rep2_R1_L001_001.fastq.gz,/data/AAVS1_16_Rep2_R2_L001_001.fastq.gz
AAVS1,CONTROL,16,2,GGGGCCACTAGGGACAGGATNGG,chr19:55115749-55115771,/data/AAVS1_16_Rep2_R1_L002_001.fastq.gz,/data/AAVS1_16_Rep2_R2_L002_001.fastq.gz
AAVS1,CONTROL,16,3,GGGGCCACTAGGGACAGGATNGG,chr19:55115749-55115771,/data/AAVS1_16_Rep3_R1_L001_001.fastq.gz,/data/AAVS1_16_Rep3_R2_L001_001.fastq.gz
AAVS1,CONTROL,16,3,GGGGCCACTAGGGACAGGATNGG,chr19:55115749-55115771,/data/AAVS1_16_Rep3_R1_L002_001.fastq.gz,/data/AAVS1_16_Rep3_R2_L002_001.fastq.gz
```

## Dependencies

The following software are required to run the SITE-Seq assay analysis pipeline:

- [Nextflow](https://www.nextflow.io/) version ≥24.04.2
- [Docker](https://www.docker.com/) version ≥20

Additional software packages are installed automatically during build. See [CITATIONS.md](CITATIONS.md) for a complete list.

## Running the pipeline

To run the SITE-Seq® assay analysis pipeline you will first need to clone this repository. Because the SITE-Seq® assay pipeline Docker image currently isn't hosted in any public registry, you will also need to build the Docker image once before running the pipeline for the first time and once for every new release of the pipeline. Here are commands to clone the repo, build the Docker image, and run the pipeline using the test profile:

```bash
git clone https://github.com/Caribou-Biosciences/nf-SITE-Seq.git && cd nf-SITE-Seq/toolkits
bash build_image.sh && cd ..
nextflow run main.nf -profile test,docker --outdir example_results
```

The test profile runs the pipeline on a subset of sequencing reads from a SITE-Seq assay run used to detect on- and off-target cleavage using Cas9 with an AAVS1 targeting guide. Please note that the majority of reads were removed from these FASTQ files to minimize download size, but the resulting report is similar to what would be seen in a normal experiment.

The typical command for running the pipeline on your own data is as follows:

```bash
nextflow run main.nf --input ./samplesheet.csv --outdir ./results --fasta /path/to/genome.fasta -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run main.nf -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
fasta: '/path/to/genome.fasta'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## Pipeline parameters

Below are descriptions of the pipeline parameters.

### Required parameters

These parameters are required for every run:

- `--input (String, required)`: Path to the samplesheet containing sample information.
- `--outdir (String, required)`: Path to the directory where outputs will be published.
- `--fasta (String, required)`: Path to the genome reference FASTA file to align input reads against.

### Utility parameters

These parameters save computation time or create more useful outputs:

- `--bowtie2_index (String, optional)`: Path to a directory containing a pre-built Bowtie2 index. Should be built from the same reference genome as `--fasta`. If left empty, an index will be generated from `--fasta`.
- `--save_bowtie2_index (Boolean, optional)`: Whether to save the Bowtie2 index created during pipeline execution to the output folder. Only relevant when `--bowtie2_index` is not specified. Destination folder is `<OUTDIR>/bowtie2_index`.
- `--use_relative_igv_paths (Boolean, optional)`: Whether to use relative paths when defining file paths in IGV session files. This parameter is ignored if `--file_server_url` is set.
- `--file_server_url (String, optional)`: The URL to a file server that serves files present in `--outdir`. Used to generate URLs pointing to output files within IGV sessions.
- `--file_server_prefixes (String, optional)`: A comma-separated string of path prefixes to replace in output file paths with `--file_server_url`. Used to generate URLs pointing to output files within IGV sessions. Example: `--file_server_prefixes /data/prefix1,/data/prefix2`.

### Configurable analysis parameters suitable for customization

These parameters modify overall analysis strategy and site calling/reporting criteria. The default settings should work well in most cases:

- `--test_pval_calling_thresh (Float, optional, default: 1e-5)`: The p-value threshold below which sites will be called.
- `--test_pval_reporting_thresh (Float, optional, default: 1e-4)`: The p-value threshold below which sites will be reported, i.e. included in the final outputs.
- `--ctrl_noise_pval_calling_thresh (Float, optional, default: 0.005)`: The control background p-value threshold below which sites will not be called due to high background, regardless of the test p-value.
- `--include_noise_pval_fails (Boolean, optional, default: false)`: Whether to include sites that fail the high-background control filter in the output.
- `--apply_normalization (Boolean, optional, default: true)`: Whether to use normalization when counting read pileups. Normalization is performed per-million.
- `--use_umis (Boolean, optional, default: true)`: Whether to use UMIs to deduplicate alignments. Requires that `fastq_2` be included in the input samplesheet.

### Critical analysis parameters you probably don't want to touch

These parameters affect the basic building blocks of the algorithm. Modification is not recommended:

- `--umi_length (Int, optional, default: 7)`: The length of the sequence to take from the start of read 2 as the UMI.
- `--umi_base_qual_thresh (Int, optional, default: 20)`: The Q-score threshold to use when filtering UMIs. At most `--max_lowq_umi_bases` UMI bases can be below this threshold, otherwise the read is filtered out.
- `--max_lowq_umi_bases (Int, optional, default: 1)`: The maximum number of UMI bases with a quality score below `--umi_base_qual_thresh` before filtering out a read.
- `--stat_window (Int, optional, default: 5)`: The window size (in base pairs) used to calculate signal statistics. You may want to increase this parameter if your nuclease of choice produces staggered cuts with a large distance between the two nicks. The default setting is sufficient for the staggered cuts produced by Cas12a.
- `--peak_window (Int, optional, default: 5)`: The window size (in base pairs) used to identify peaks of the signal statistic.
- `--control_noise_window (Int, optional, default: 150)`: The window size (in base pairs) used to calculate high-background signals in the control samples.
- `--agg_peak_window (Int, optional, default: 10)`: The window size (in base pairs) used to aggregate peaks across RNP concentrations.
- `--min_aln_length (Integer, optional, default: 20)`: The minimum reference length that an alignment must have to be included in the analysis.
- `--max_clipping (Integer, optional, default: 0)`: The maximum allowed clipping at the start of the read for an alignment to be included in the analysis.
- `--min_map_qual (Integer, optional, default: 5)`: The minimum mapping quality that an alignment must have to be included in the analysis.
- `--motif_search_factor (Integer, optional, default: 1)`: Used to define the size of the window around a signal peak to search for motif sequences. This value is multiplied by the length of the on-target motif which determines how far upstream and downstream from a peak to search for motif sequences.
- `--motif_match_score (Integer, optional, default: 1)`: The score for a matched base when searching for motif sequences.
- `--motif_mismatch_pen (Integer, optional, default: 1)`: The penalty for a mismatched base when searching for motif sequences. Should be positive.
- `--motif_gap_pen (Integer, optional, default: 8)`: The penalty for a gap when searching for motif sequences. Should be positive.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment. Below are the profiles currently available in the pipeline:

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/). This is currently the only supported containerization strategy for this pipeline.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run Nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
