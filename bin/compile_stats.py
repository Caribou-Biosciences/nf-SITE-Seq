#!/usr/bin/env python3

############################################################################
############################################################################
## Parse various logs and reports for a single sample and output a        ##
## single-row TSV containing read/UMI/alignment statistics for the sample ##
## calling summary, a site threshold summary, and a table of site calls   ##
############################################################################
############################################################################

import argparse
import json
from typing import Dict, Optional
import csv

from site_seq_utils.constants import SITESeqQCColumn


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--conc", type=int, required=True)
    parser.add_argument("--identifier", type=str, required=True)
    parser.add_argument("--sample_group", type=str, required=True)
    parser.add_argument("--rep", type=int, required=True)
    parser.add_argument("--out_path", type=str, required=True)
    parser.add_argument("--fastp_json", type=str, required=True)
    parser.add_argument("--bowtie2_log", type=str, required=True)
    parser.add_argument("--count_reads_json", type=str, required=True)
    parser.add_argument("--umi_filter_report", type=str, required=False)
    parser.add_argument("--umitools_dedup_log", type=str, required=False)
    args = parser.parse_args()

    metadata = {
        SITESeqQCColumn.IDENTIFIER: args.identifier,
        SITESeqQCColumn.SAMPLE_GROUP: args.sample_group,
        SITESeqQCColumn.REPLICATE: args.rep,
        SITESeqQCColumn.CONCENTRATION: args.conc,
    }
    stats = {}
    stats.update(parse_fastp_json(args.fastp_json))
    stats.update(parse_umi_filter_report(args.umi_filter_report))
    stats.update(parse_bowtie2_log(args.bowtie2_log))
    stats.update(parse_umitools_dedup_log(args.umitools_dedup_log))
    stats.update(parse_count_reads_json(args.count_reads_json))

    # Calculate derived stats
    stats.update(
        {
            SITESeqQCColumn.PCT_READS_PF: (
                100
                * stats[SITESeqQCColumn.READS_PF]
                / stats[SITESeqQCColumn.TOTAL_READS]
                if stats[SITESeqQCColumn.TOTAL_READS]
                else None
            ),
        }
    )

    column_order = [i.value for i in SITESeqQCColumn]
    sorted_stats = {
        k: v for k, v in sorted(stats.items(), key=lambda t: column_order.index(t[0]))
    }
    row = {**metadata, **sorted_stats}
    with open(args.out_path, "w") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def parse_fastp_json(json_path: str) -> Dict:
    """Parse a Fastp JSON report"""
    with open(json_path, "r") as fh:
        json_data = json.load(fh)

    stats = {
        SITESeqQCColumn.TOTAL_READS: json_data["summary"]["before_filtering"][
            "total_reads"
        ]
        // 2,
        SITESeqQCColumn.LOWQ_READS: json_data["filtering_result"]["low_quality_reads"]
        // 2,
        SITESeqQCColumn.TOO_MANY_N_READS: json_data["filtering_result"][
            "too_many_N_reads"
        ]
        // 2,
        SITESeqQCColumn.SHORT_READS: json_data["filtering_result"]["too_short_reads"]
        // 2,
        SITESeqQCColumn.READS_PF: json_data["filtering_result"]["passed_filter_reads"]
        // 2,
        SITESeqQCColumn.PCT_Q30_READS_PF: 100
        * json_data["summary"]["before_filtering"]["q30_rate"],
    }
    return stats


def parse_bowtie2_log(log_path: str) -> Dict:
    """Parse a Bowtie2 alignment log"""
    total_reads = unaligned_reads = None
    with open(log_path, "r") as fh:
        for line in fh.readlines():
            if line.strip().endswith(" reads; of these:"):
                total_reads = int(line.split(maxsplit=1)[0])
            elif line.strip().endswith(" aligned 0 times"):
                unaligned_reads = int(line.strip().split(maxsplit=1)[0])

    if total_reads is None or unaligned_reads is None:
        raise RuntimeError(f"Could not parse alignment stats from {log_path}")

    aligned_reads = total_reads - unaligned_reads
    stats = {
        SITESeqQCColumn.ALIGNED_READS: aligned_reads,
        SITESeqQCColumn.PCT_ALIGNED: 100 * aligned_reads / total_reads,
    }
    return stats


def parse_umitools_extract_log(log_path: Optional[str]) -> Dict:
    """Parse a UMI-tools extract log"""
    if log_path is None:
        stats = {SITESeqQCColumn.LOWQ_UMIS: "N/A"}
        return stats
    lowq_umis = None
    with open(log_path, "r") as fh:
        for line in fh.readlines():
            if "Reads output:" in line:
                output_reads = int(line.strip().split("Reads output: ")[-1])
            elif "filtered: umi quality:" in line:
                lowq_umis = int(line.strip().split("filtered: umi quality: ")[-1])

    if lowq_umis is None:
        raise RuntimeError(f"Could not parse dedup stats from {log_path}")
    stats = {
        SITESeqQCColumn.LOWQ_UMIS: lowq_umis,
        SITESeqQCColumn.READS_PF: output_reads,
    }
    return stats


def parse_umitools_dedup_log(log_path: Optional[str]) -> Dict:
    """Parse a UMI-tools dedup log"""
    if log_path is None:
        stats = {
            SITESeqQCColumn.TOTAL_DEDUP_UMIS: "N/A",
            SITESeqQCColumn.PCT_DUPLICATION: "N/A",
        }
        return stats
    total_umis = total_dedup_umis = None
    with open(log_path, "r") as fh:
        for line in fh.readlines():
            if "Reads: Input Reads:" in line:
                total_umis = int(line.strip().split("Reads: Input Reads: ")[-1])
            if "Number of reads out:" in line:
                total_dedup_umis = int(line.strip().split("Number of reads out: ")[-1])

    if total_umis is None or total_dedup_umis is None:
        raise RuntimeError(f"Could not parse dedup stats from {log_path}")
    pct_duplication = (
        100 * (total_umis - total_dedup_umis) / total_umis if total_umis else None
    )
    stats = {
        SITESeqQCColumn.TOTAL_DEDUP_UMIS: total_dedup_umis,
        SITESeqQCColumn.PCT_DUPLICATION: pct_duplication,
    }
    return stats


def parse_umi_filter_report(json_path: Optional[str]) -> Dict:
    """Parse a UMI filtering report"""
    if json_path is None:
        stats = {SITESeqQCColumn.LOWQ_UMIS: "N/A"}
        return stats

    with open(json_path, "r") as fh:
        json_data = json.load(fh)

    stats = {
        SITESeqQCColumn.LOWQ_UMIS: json_data["failing_reads"],
        SITESeqQCColumn.READS_PF: json_data["passing_reads"],
    }
    return stats


def parse_count_reads_json(json_path: str) -> Dict:
    """Parse a count-read-starts JSON report"""
    with open(json_path, "r") as fh:
        json_data = json.load(fh)

    stats = {
        SITESeqQCColumn.CLIPPED_ALIGNMENTS: json_data["filters"]["above_clipping"],
        SITESeqQCColumn.LOW_MAPQ_ALIGNMENTS: json_data["filters"][
            "low_mapping_quality"
        ],
        SITESeqQCColumn.ALIGNMENTS_PF: json_data["passing_reads"],
        SITESeqQCColumn.TOTAL_UNIQUE_SITES: json_data["total_sites"],
    }
    return stats


if __name__ == "__main__":
    main()
