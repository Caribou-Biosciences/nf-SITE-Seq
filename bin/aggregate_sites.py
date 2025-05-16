#!/usr/bin/env python3

#########################################################################
#########################################################################
## Aggregate reported sites among different RNP concentrations in the  ##
## same sample group                                                   ##
#########################################################################
#########################################################################


import argparse
import os
import attrs
import numpy as np
import pandas as pd
from typing import Tuple, Dict, List

from site_seq_utils.site import Site, AggregateSite
from site_seq_utils.dist_tsv_parser import DistTSVParser
from site_seq_utils.site_calculator import SiteCalculator


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--identifier", type=str, required=True)
    parser.add_argument("--output_tsv", type=str, required=True)
    parser.add_argument("--reported_site_tsvs", type=str, nargs="+", required=True)
    parser.add_argument("--test_window_sums_tsvs", type=str, nargs="+", required=True)
    parser.add_argument("--test_window_diffs_tsvs", type=str, nargs="+", required=True)
    parser.add_argument("--control_window_sums_tsv", type=str, required=True)
    parser.add_argument("--control_wide_window_sums_tsv", type=str, required=True)
    parser.add_argument("--ctrl_peaks_dist_tsv", type=str, required=True)
    parser.add_argument("--ctrl_wide_window_sum_dist_tsv", type=str, required=True)
    parser.add_argument("--site_collapse_window", type=int, required=True)
    parser.add_argument("--test_pval_reporting_thresh", type=float, required=True)
    parser.add_argument("--test_pval_calling_thresh", type=float, required=True)
    parser.add_argument("--control_pval_calling_thresh", type=float, required=True)
    parser.add_argument("--include_control_pval_fails", action="store_true")
    args = parser.parse_args()

    (
        sites_by_conc,
        test_signal_parsers_by_conc,
        test_window_sums_parsers_by_conc,
    ) = parse_and_validate_input_files(
        args.reported_site_tsvs,
        args.test_window_sums_tsvs,
        args.test_window_diffs_tsvs,
    )
    control_window_sum_parser = DistTSVParser(args.control_window_sums_tsv)
    control_wide_window_sum_parser = DistTSVParser(args.control_wide_window_sums_tsv)

    site_calculator = SiteCalculator(
        control_peak_dist_path=args.ctrl_peaks_dist_tsv,
        wide_control_peak_dist_path=args.ctrl_wide_window_sum_dist_tsv,
        peak_calling_thresh=args.test_pval_calling_thresh,
        background_filter_thresh=args.control_pval_calling_thresh,
        peak_reporting_thresh=args.test_pval_reporting_thresh,
        include_background_peaks=args.include_control_pval_fails,
        terminate_early=False,
    )

    # Create a mapping between every site location across all concentrations and the
    # reported site at each concentration (if any)
    concs = sorted(sites_by_conc)
    location_site_map = {}
    for conc, sites in sites_by_conc.items():
        for conc_site in sites:
            loc_key = (conc_site.chrom, conc_site.position)
            location_site_map.setdefault(loc_key, {})[conc] = conc_site

    # Go through each site location reported across all concentrations and create
    # an AggregateSite object. Calculate site values for concentrations with no
    # reported site at that location
    aggregate_site_map = {}
    for (chrom, pos), conc_site_map in location_site_map.items():
        control_window_vals = control_window_sum_parser.get_values_at_position(
            chrom, pos, default_val=0.0
        )
        control_wide_window_vals = (
            control_wide_window_sum_parser.get_values_at_position(
                chrom, pos, default_val=0.0
            )
        )
        for conc in concs:
            if conc not in conc_site_map:
                # This concentration does not have a reported site at this location,
                # calculate the site values
                signals = test_signal_parsers_by_conc[conc].get_values_at_position(
                    chrom, pos, default_val=0.0
                )
                test_window_vals = test_window_sums_parsers_by_conc[
                    conc
                ].get_values_at_position(chrom, pos, default_val=0.0)
                conc_site = site_calculator.calculate_site(
                    chrom,
                    pos,
                    signals,
                    test_window_vals,
                    control_window_vals,
                    control_wide_window_vals,
                )
                conc_site_map[conc] = conc_site

        aggregate_site = AggregateSite(site_dict=conc_site_map)
        aggregate_site_map[(chrom, pos)] = aggregate_site

    # Now that we have an aggregate site at each reported site position,
    # we use a sliding window to collapse nearby sites into a single site call
    # (different concentrations might differ slightly on the cut-site for a given site).
    # The aggregate site with the highest value in each window is chosen. Value criteria
    # are encoded in the AggregateSite class.
    all_locs = sorted(aggregate_site_map)
    best_aggregated_sites = []
    prev_chrom = prev_pos = None
    best_agg_site = None
    for cur_chrom, cur_pos in all_locs:
        cur_agg_site = aggregate_site_map[(cur_chrom, cur_pos)]
        if prev_chrom != cur_chrom or prev_pos + args.site_collapse_window < cur_pos:
            if best_agg_site:
                best_aggregated_sites.append(best_agg_site)
            best_agg_site = cur_agg_site
        elif best_agg_site is None or cur_agg_site > best_agg_site:
            best_agg_site = cur_agg_site

        prev_chrom, prev_pos = cur_chrom, cur_pos

    if best_agg_site:
        best_aggregated_sites.append(best_agg_site)

    # Sort aggregated sites. Comparison criteria encoded in AggregateSite class.
    best_aggregated_sites = sorted(
        best_aggregated_sites,
        reverse=True,
    )

    # Add site ranks and recalulate % of max peak stat since the max peak for a given
    # concentration might have shifted slightly
    on_target_found = any(site.is_on_target for site in best_aggregated_sites)
    ranks = (
        list(range(len(best_aggregated_sites)))
        if on_target_found
        else list(range(1, len(best_aggregated_sites) + 1))
    )
    max_signals_by_conc = {
        conc: [0.0] * test_signal_parsers_by_conc[conc].num_samples for conc in concs
    }
    for agg_site in best_aggregated_sites:
        for conc in concs:
            cur_maxes = max_signals_by_conc[conc]
            for idx in range(test_signal_parsers_by_conc[conc].num_samples):
                cur_maxes[idx] = max(
                    agg_site.site_dict[conc].signal_vals[idx], cur_maxes[idx]
                )

    for idx, agg_site in enumerate(best_aggregated_sites):
        agg_site = attrs.evolve(agg_site, rank=ranks[idx])
        best_aggregated_sites[idx] = agg_site
        for conc in concs:
            conc_site = agg_site.site_dict[conc]
            pct_max_peak_vals = [
                100 * sig / max_sig if max_sig else np.nan
                for sig, max_sig in zip(
                    conc_site.signal_vals, max_signals_by_conc[conc]
                )
            ]
            agg_site.site_dict[conc] = attrs.evolve(
                conc_site, pct_max_peak_vals=pct_max_peak_vals
            )

    # Define output columns and ordering
    site_fields = attrs.fields(Site)
    agg_site_fields = attrs.fields(AggregateSite)
    meta_cols = [
        site_fields.chrom.name,
        site_fields.position.name,
        site_fields.position_1b.name,
        site_fields.closest_motif.name,
        site_fields.is_on_target.name,
        site_fields.motif_location.name,
        site_fields.motif_score.name,
        site_fields.num_subs.name,
        site_fields.num_gaps.name,
    ]
    site_data_fields = [
        agg_site_fields.total_signal,
        site_fields.site_called,
        site_fields.test_pval,
        site_fields.background_pval,
        site_fields.mean_pct_max_peak,
        site_fields.pct_max_peak_vals,
        site_fields.mean_signal,
        site_fields.signal_vals,
        site_fields.total_test_reads,
        site_fields.total_control_reads,
        site_fields.total_wide_control_reads,
        site_fields.mean_wide_control_reads,
        site_fields.test_read_vals,
        site_fields.control_read_vals,
        site_fields.wide_control_read_vals,
        site_fields.total_signal,
        site_fields.digested,
        site_fields.high_background,
    ]
    site_data_cols = [
        col
        for data_field in site_data_fields
        for col in AggregateSite.get_field_keys(data_field, concs)
    ]
    columns = [
        agg_site_fields.rank.name,
        *meta_cols,
        *site_data_cols,
    ]

    # Convert to DataFrame
    df_output = pd.DataFrame(
        [site.to_dict() for site in best_aggregated_sites], columns=columns
    )
    with open(args.output_tsv, "w") as fh:
        fh.write(f"#identifier={args.identifier}\n")
        fh.write(f"#concentrations=[{', '.join(map(str, sorted(concs)))}]\n")
        df_output.to_csv(fh, sep="\t", index=False)


def parse_and_validate_input_files(
    reported_site_tsvs: List[str],
    test_window_sums_tsvs: List[str],
    test_window_diffs_tsvs: List[str],
) -> Tuple[Dict[int, List[Site]], Dict[int, DistTSVParser], Dict[int, DistTSVParser]]:
    sites_by_conc, site_paths_by_conc = {}, {}
    test_window_sums_parsers_by_conc = {}
    test_signal_parsers_by_conc = {}
    for site_tsv_path in reported_site_tsvs:
        with open(site_tsv_path, "r") as fh:
            identifier = conc = None
            line = fh.readline()
            while line.startswith("#"):
                key, val = line.strip("#\n").split("=")
                if key == "concentration":
                    conc = int(val)
                elif key == "identifier":
                    identifier = val

                line = fh.readline()

            if identifier is None or conc is None:
                raise Exception(
                    f"Could not parse metadata from reported site TSV {site_tsv_path}"
                )
            elif conc in sites_by_conc:
                raise Exception(
                    f"More than one site TSV found for concentration {conc}: "
                    f"{site_tsv_path}, {site_paths_by_conc[conc]}"
                )

            columns = line.strip().split("\t")
            df = pd.read_csv(fh, sep="\t", names=columns)
            sites = [Site.load(row) for _, row in df.iterrows()]
            sites_by_conc[conc] = sites
            site_paths_by_conc[conc] = site_tsv_path

            conc_key = f"_{conc}nM"
            test_window_sum_matches = [
                p for p in test_window_sums_tsvs if conc_key in os.path.basename(p)
            ]
            if len(test_window_sum_matches) > 1:
                raise Exception(
                    f"Multiple matching files for {conc_key}: {', '.join(test_window_sum_matches)}"
                )
            elif not test_window_sum_matches:
                raise Exception(f"No matching test window sum TSV for {conc_key}")
            else:
                test_window_sum_tsv = test_window_sum_matches[0]
                test_window_sum_parser = DistTSVParser(test_window_sum_tsv)
                test_window_sums_parsers_by_conc[conc] = test_window_sum_parser

            test_window_diff_matches = [
                p for p in test_window_diffs_tsvs if conc_key in os.path.basename(p)
            ]
            if len(test_window_diff_matches) > 1:
                raise Exception(
                    f"Multiple matching files for {conc_key}: {', '.join(test_window_diff_matches)}"
                )
            elif not test_window_diff_matches:
                raise Exception(f"No matching test window diff TSV for {conc_key}")
            else:
                test_window_diff_tsv = test_window_diff_matches[0]
                test_window_diff_parser = DistTSVParser(test_window_diff_tsv)
                test_signal_parsers_by_conc[conc] = test_window_diff_parser

    return (
        sites_by_conc,
        test_signal_parsers_by_conc,
        test_window_sums_parsers_by_conc,
    )


if __name__ == "__main__":
    main()
