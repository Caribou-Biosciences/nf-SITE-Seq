# Caribou-Biosciences/nf-SITE-Seq: Site-Calling Algorithm

## Introduction

This document describes the core site-calling algorithm of the SITE-Seq® assay analysis pipeline. The algorithm does not include the upstream read preprocessing, alignment, and deduplication.

## The SITE-Seq® assay site-calling algorithm

The input to the algorithm is a BAM file containing alignments of reads to the reference genome that has been deduplicated to remove PCR duplicates.

### Count read starts

The first step of the algorithm is to parse the input BAM file and count how many reads start at each position in the reference genome. Under default parameters, these counts are normalized by the per-million observations in that sample (normalization can be disabled using `--apply_normalization false`). For a read to be counted, its alignment must pass three filters:

- No more than `--max_clipping` (default: `0`) clipped bases at the start of the alignment
- Have a mapping quality no less than `--min_map_qual` (default: `5`)
- Span at least `--min_aln_length` (default: `20`) reference bases

### Merge read starts

Next, the read start counts for all replicates of the same sample group and RNP concentration are merged into a single file. No manipulations are made to the counts themselves, but it's good to highlight that subsequent steps are performed on all replicates at once.

### Calculate corrected signal

Next, each test replicate (i.e., RNP concentration > 0 nM) is matched with a control replicate (i.e., RNP concentration = 0 nM) to calculate a corrected signal value at each position in the reference genome. For each position $x$ in the genome, this signal is calculated as:

$$
Signal_x =
\begin{cases}
  \sum\limits_{i=x-w}^{x+w} T_i - \sum\limits_{i=x-w}^{x+w} C_i, & \text{if } \sum\limits_{i=x-w}^{x+w} T_i - \sum\limits_{i=x-w}^{x+w} C_i \ge 0 \\
  0, & \text{if } \sum\limits_{i=x-w}^{x+w} T_i - \sum\limits_{i=x-w}^{x+w} C_i \lt 0
\end{cases}
$$

where $T_i$ and $C_i$ are the number of reads counted at position $i$ in the test and control replicate respectively, and $w$ is the window size (controlled by `--stat_window`, `5` by default).

### Calculate control vs control signal

To contextualize the values observed in the [test vs control signal calculation](#calculate-corrected-signal), the same signal calculation is performed on paired control replicates. If 3 control replicates are included, $C_1$ is compared to $C_2$, $C_2$ is compared to $C_3$, and $C_3$ is compared to $C_1$. This calculation is why at least 2 control replicates are required to run the pipeline.

### Calculate background signal

The signal calculation is sensitive to high-backround regions of the genome. To avoid calling false positives in these regions, a simple sliding window sum is performed on the read start counts of the control samples using a wider window (controlled by `--control_noise_window`, `150` by default). These background signal values feed into the site-calling step downstream.

### Identify signal peaks

The signal values of both the [test vs control](#calculate-corrected-signal) and [control vs control](#calculate-control-vs-control-signal) signal calculations are processed to identify signal peaks. A simple sliding window is passed over the signal values and any continuous chain of windows containing a non-zero signal are collapsed into the single position with the maximum value (ties are given to the lesser coordinate).

### Report and call sites

We're finally ready to call some sites! First, two eCDFs are constructed, one from the [control vs control signal peaks](#calculate-control-vs-control-signal) and one from the [control background signal values](#calculate-background-signal). Then, each [test vs control signal peak](#calculate-corrected-signal) is evaluated to determine whether it should be "reported", i.e. included in the output, and "called", i.e. labeled as a potential off-target site. A p-value is calculated by comparing the mean peak signal strength across all replicates to the eCDF of control vs control signal peaks. An important note is that if a signal value exceeds the maximum value observed in the eCDF, the p-value is set to $\frac{1}{N}$ where $N$ is the total number of observations in the control vs control distribution. This p-value is compared to the site reporting (`--test_pval_reporting_thresh`, `1e-4` by default) and site calling (`--test_pval_calling_thresh`, `1e-5` by default) thresholds; if the p-value is below these thresholds, the site has passed the first step of the site reporting and site calling criteria respectively. Next, the same genomic position is evaluated in the control to determine whether it lies in a high-background region. The mean of the [background signal values](#calculate-background-signal) at that location is calculated and compared to the eCDF of control background signal values. If the value of the background signal exceeds the point at which the inverse eCDF of the control background signal reaches `--ctrl_noise_pval_calling_thresh` (`5e-3` by default), that site is labeled as high-background and is excluded from the output, regardless of signal strength in the test samples. To include sites that fail the high-background check, use the `--include_noise_pval_fails` parameter. Once the sites passing the reporting criteria are identified, they are sorted to produce a ranked list. The on-target (if present in the list) is given a rank of 0, and remaining sites are ranked by their mean corrected signal strength.

The output of this step is a ranked list of sites that have passed the reporting criteria, some, all, or none of which have also passed the calling criteria. One list of sites is produced for each combination of sample group and RNP concentration, regardless of the number of replicates included.

### Identify off-target motifs

For each site passing the reporting criteria, the sequence motif best-matching the on-target motif within close proximity of the site is identified. The value of the `--motif_search_factor` (default: `1`) is multiplied by the length of the on-target motif to define the number of base pairs upstream and downstream of the site to search within. An alignment between this sequence and the on-target motif is performed using the `--motif_match_score` (default: `1`), `--motif_mismatch_pen` (default: `1`), `--motif_gap_pen` (default: `8`) to set the match score, mismatch penalty, and gap penalty respectively. The reference sequence of the alignment with the highest alignment score is chosen. It's worth highlighting that all motif information is purely annotative and is not used to select or filter sites.

### Aggregate sites

We've called sites in each RNP concentration in isolation, but the true power of the SITE-Seq® assay is in combining signal information across multiple RNP concentrations. To accomplish this, the algorithm aggregates sites reported at individual RNP concentrations into a single list of ranked sites for each sample group. Since a site may have slightly different coordinates across different RNP concentrations, or may be missing entirely from a given RNP concentration, the algorithm picks a consensus location for each site and recomputes statistic values to ensure they are accurate at the consensus position across all RNP concentrations. All sites reported by all RNP concentrations are grouped using a sliding window of length `--agg_peak_window` (default: `10`) and each unique location among the sites in a group is considered as a candidate for the consensus location. For each candidate location, signal statistics are extracted for every RNP concentration and the signal values for each candidate are compared to identify the location with the strongest "aggregate" signal. A tuple is constructed for each candidate and is populated with a signal value for each concentration in ascending order; this value is the mean corrected signal at that RNP concentration if the site was called at that concentration and all higher concentrations and is 0 otherwise. The final item in the tuple is the total corrected signal across all RNP concentrations, regardless of site calls. A component-wise comparison is performed between tuples such that the first element with differing values determines which tuple is of greater value. The candidate location with the highest value tuple is chosen as the consensus location for the site. Finally, aggregated sites are compared to determine the final site ranking. The on-target, if included in the reported sites, is given a rank of 0, and the remaining sites are sorted according to the same tuple comparison process that was used to pick consensus locations.
