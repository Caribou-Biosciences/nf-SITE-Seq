#!/usr/bin/env python3

#########################################################################
#########################################################################
## Parse read and signal distributions for one or more replicates of   ##
## and RNP concentration and compare peak signal strengths to the      ##
## control distributions in order to determine which sites should be   ##
## reported and which sites should be "called". Output a TSV of all    ##
## reported/called sites containing statistic values and best-match    ##
## motif information.                                                  ##
#########################################################################
#########################################################################


import argparse
import attrs
import pandas as pd

from site_seq_utils.site import Site
from site_seq_utils.coordinate import Zerodinate, Onedinate
from site_seq_utils.dist_tsv_parser import DistTSVParser
from site_seq_utils.motif_finder import MotifFinder
from site_seq_utils.site_calculator import SiteCalculator


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--conc", type=int, required=True)
    parser.add_argument("--identifier", type=str, required=True)
    parser.add_argument("--output_tsv", type=str, required=True)
    parser.add_argument("--test_peaks_tsv", type=str, required=True)
    parser.add_argument("--test_window_sums_tsv", type=str, required=True)
    parser.add_argument("--control_window_sums_tsv", type=str, required=True)
    parser.add_argument("--control_wide_window_sums_tsv", type=str, required=True)
    parser.add_argument("--ctrl_peaks_dist_tsv", type=str, required=True)
    parser.add_argument("--ctrl_wide_window_sum_dist_tsv", type=str, required=True)
    parser.add_argument("--test_pval_reporting_thresh", type=float, required=True)
    parser.add_argument("--test_pval_calling_thresh", type=float, required=True)
    parser.add_argument("--control_pval_calling_thresh", type=float, required=True)
    parser.add_argument("--include_control_pval_fails", action="store_true")
    parser.add_argument("--fasta_path", type=str, required=True)
    parser.add_argument("--search_motif", type=str, required=True)
    parser.add_argument("--on_target_location", type=str, required=True)
    parser.add_argument("--search_factor", type=float, required=True)
    parser.add_argument("--motif_match_score", type=int, required=True)
    parser.add_argument("--motif_mismatch_pen", type=int, required=True)
    parser.add_argument("--motif_gap_pen", type=int, required=True)
    args = parser.parse_args()

    # Helper object for finding the best-match sequence to the on-target motif near
    # a particular location
    motif_finder = MotifFinder(
        args.search_motif,
        args.fasta_path,
        args.search_factor,
        args.motif_match_score,
        args.motif_mismatch_pen,
        args.motif_gap_pen,
    )
    # Helper object for calculating site information based on test/control stat values
    site_calculator = SiteCalculator(
        control_peak_dist_path=args.ctrl_peaks_dist_tsv,
        wide_control_peak_dist_path=args.ctrl_wide_window_sum_dist_tsv,
        peak_calling_thresh=args.test_pval_calling_thresh,
        background_filter_thresh=args.control_pval_calling_thresh,
        peak_reporting_thresh=args.test_pval_reporting_thresh,
        include_background_peaks=args.include_control_pval_fails,
        motif_finder=motif_finder,
        on_target_loc=Onedinate(args.on_target_location),
        terminate_early=True,
    )

    # Objects for parsing tabix-indexed TSV files of stat distributions
    peak_parser = DistTSVParser(args.test_peaks_tsv)
    test_window_sum_parser = DistTSVParser(args.test_window_sums_tsv)
    control_window_sum_parser = DistTSVParser(args.control_window_sums_tsv)
    control_wide_window_sum_parser = DistTSVParser(args.control_wide_window_sums_tsv)

    # Validate that the samples in the input distributions match each other
    comparison_names = peak_parser.columns[2:]
    test_names = test_window_sum_parser.columns[2:]
    control_names = control_window_sum_parser.columns[2:]
    ww_control_names = control_wide_window_sum_parser.columns[2:]
    if len(comparison_names) != len(test_names):
        raise ValueError(
            f"Column count mismatch(es) in files {args.test_peaks_tsv}, "
            f"{args.test_window_sums_tsv}"
        )
    elif len(test_names) > len(control_names):
        raise ValueError(
            f"More columns in {args.test_peaks_tsv} and {args.test_window_sums_tsv} "
            f"than in {args.control_window_sums_tsv}"
        )
    elif control_names != ww_control_names:
        raise ValueError(
            f"Columns don't match in control files {args.control_window_sums_tsv}, "
            f"{args.control_wide_window_sums_tsv}"
        )

    # Initialize iterators for all records in the distribution TSVs
    peak_iter = peak_parser.fetch()
    test_window_sum_iter = test_window_sum_parser.fetch()
    control_window_sum_iter = control_window_sum_parser.fetch()
    control_wide_window_sum_iter = control_wide_window_sum_parser.fetch()

    # Get the first records for the distribution TSVs (other than the test peak TSV)
    try:
        cur_test_window_rec = next(test_window_sum_iter)
    except StopIteration:
        cur_test_window_rec = None
    try:
        cur_control_window_rec = next(control_window_sum_iter)
    except StopIteration:
        cur_control_window_rec = None
    try:
        cur_control_wide_window_rec = next(control_wide_window_sum_iter)
    except StopIteration:
        cur_control_wide_window_rec = None

    # Iterate through every record in the test peak distribution TSV, get the overlapping
    # record for each other distribution TSV (if any), and use the distribution values
    # at each position to determine whether the peak should become a reported/called site.
    # Since each distribution TSV is sorted, we can parse each record from the peak
    # iterator and advance each other iterator until we reach or pass the location of the
    # peak, thus we only have to parse each TSV once.
    sites = []
    max_peak_vals = [0] * peak_parser.num_samples
    on_target_found = False
    while True:
        try:
            peak_record = next(peak_iter)
        except StopIteration:
            break

        test_signals = [peak_record[n] for n in peak_parser.sample_names]
        if sum(test_signals) < 1:
            continue

        peak_loc = Zerodinate(
            peak_record["chrom"], peak_record["start"], peak_record["end"]
        )

        if cur_test_window_rec:
            # Advance this iterator until we reach or pass the current peak location
            while peak_record["chrom"] > cur_test_window_rec["chrom"] or (
                peak_record["chrom"] == cur_test_window_rec["chrom"]
                and peak_record["start"] >= cur_test_window_rec["end"]
            ):
                try:
                    cur_test_window_rec = next(test_window_sum_iter)
                except StopIteration:
                    cur_test_window_rec = None
                    break

            test_window_loc = (
                Zerodinate(
                    cur_test_window_rec["chrom"],
                    cur_test_window_rec["start"],
                    cur_test_window_rec["end"],
                )
                if cur_test_window_rec
                else None
            )

        if test_window_loc and peak_loc.is_overlapping_with(test_window_loc):
            # Record overlaps the peak record, use its stat values
            test_window_vals = [
                cur_test_window_rec[n] for n in test_window_sum_parser.sample_names
            ]
        else:
            # No overlapping record, default to zeros
            test_window_vals = [0] * test_window_sum_parser.num_samples

        if cur_control_window_rec:
            # Advance this iterator until we reach or pass the current peak location
            while peak_record["chrom"] > cur_control_window_rec["chrom"] or (
                peak_record["chrom"] == cur_control_window_rec["chrom"]
                and peak_record["start"] >= cur_control_window_rec["end"]
            ):
                try:
                    cur_control_window_rec = next(control_window_sum_iter)
                except StopIteration:
                    cur_control_window_rec = None
                    break

            control_window_loc = (
                Zerodinate(
                    cur_control_window_rec["chrom"],
                    cur_control_window_rec["start"],
                    cur_control_window_rec["end"],
                )
                if cur_control_window_rec
                else None
            )

        if control_window_loc and peak_loc.is_overlapping_with(control_window_loc):
            # Record overlaps the peak record, use its stat values
            control_window_vals = [
                cur_control_window_rec[n]
                for n in control_window_sum_parser.sample_names
            ]
        else:
            # No overlapping record, default to zeros
            control_window_vals = [0] * control_window_sum_parser.num_samples

        if cur_control_wide_window_rec:
            # Advance this iterator until we reach or pass the current peak location
            while peak_record["chrom"] > cur_control_wide_window_rec["chrom"] or (
                peak_record["chrom"] == cur_control_wide_window_rec["chrom"]
                and peak_record["start"] >= cur_control_wide_window_rec["end"]
            ):
                try:
                    cur_control_wide_window_rec = next(control_wide_window_sum_iter)
                except StopIteration:
                    cur_control_wide_window_rec = None
                    break

            control_wide_window_loc = (
                Zerodinate(
                    cur_control_wide_window_rec["chrom"],
                    cur_control_wide_window_rec["start"],
                    cur_control_wide_window_rec["end"],
                )
                if cur_control_wide_window_rec
                else None
            )

        if control_wide_window_loc and peak_loc.is_overlapping_with(
            control_wide_window_loc
        ):
            # Record overlaps the peak record, use its stat values
            control_wide_window_vals = [
                cur_control_wide_window_rec[n]
                for n in control_wide_window_sum_parser.sample_names
            ]
        else:
            # No overlapping record, default to zeros
            control_wide_window_vals = [0] * control_wide_window_sum_parser.num_samples

        # Now that we have the values for every overlapping record (or empty values if
        # there is none), we can calculate the site values for this peak
        site = site_calculator.calculate_site(
            peak_record["chrom"],
            peak_record["start"],
            test_signals,
            test_window_vals,
            control_window_vals,
            control_wide_window_vals,
        )
        if site is not None:
            sites.append(site)
            if site.site_called:
                for idx in range(peak_parser.num_samples):
                    max_peak_vals[idx] = max(max_peak_vals[idx], site.signal_vals[idx])
            if site.is_on_target:
                on_target_found = True

    # Sort sites. Sorting criteria are represented in Site class
    sites = sorted(sites, reverse=True)

    # Calculate site ranks. We make the on-target rank 0 if it was reported
    ranks = (
        list(range(len(sites))) if on_target_found else list(range(1, len(sites) + 1))
    )

    # Calculate % of max peak statistic
    for idx, site in enumerate(sites):
        pct_max_peak_vals = [
            100 * wd / mwd if mwd > 0 else None
            for wd, mwd in zip(site.signal_vals, max_peak_vals)
        ]
        updated_site = attrs.evolve(
            site, rank=ranks[idx], pct_max_peak_vals=pct_max_peak_vals
        )
        sites[idx] = updated_site

    # Define output columns and ordering
    fields = attrs.fields(Site)
    columns = [
        fields.rank.name,
        fields.chrom.name,
        fields.position.name,
        fields.position_1b.name,
        fields.is_on_target.name,
        fields.motif_location.name,
        fields.closest_motif.name,
        fields.motif_score.name,
        fields.num_subs.name,
        fields.num_gaps.name,
        fields.site_called.name,
        fields.digested.name,
        fields.high_background.name,
        fields.test_pval.name,
        fields.background_pval.name,
        fields.total_signal.name,
        fields.mean_signal.name,
        fields.signal_vals.name,
        fields.total_test_reads.name,
        fields.total_control_reads.name,
        fields.test_read_vals.name,
        fields.control_read_vals.name,
        fields.total_wide_control_reads.name,
        fields.mean_wide_control_reads.name,
        fields.wide_control_read_vals.name,
        fields.pct_max_peak_vals.name,
        fields.mean_pct_max_peak.name,
    ]

    # Convert to DataFrame
    df_output = pd.DataFrame([attrs.asdict(site) for site in sites], columns=columns)

    # Create output TSV
    with open(args.output_tsv, "w") as fh:
        fh.write(f"#identifier={args.identifier}\n")
        fh.write(f"#concentration={args.conc}\n")
        df_output.to_csv(fh, sep="\t", index=False)


if __name__ == "__main__":
    main()
