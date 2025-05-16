#!/usr/bin/env python3

##########################################################################
##########################################################################
## Create an HTML report summarizing various information and results    ##
## from a particular sample group, include sample metadata, a site      ##
## calling summary, a site threshold summary, and a table of site calls ##
##########################################################################
##########################################################################

import argparse
from typing import List, Callable, Dict, Tuple
import attrs
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from site_seq_utils.site import Site, AggregateSite
from site_seq_utils.formatters import (
    format_integer,
    format_decimal,
    format_pvalue,
    format_percentage,
    format_dna_sequence,
)
from site_seq_utils.site_calculator import SiteCalculator
from site_seq_utils.constants import SITESeqQCColumn

PLOTLY_THEME = "seaborn"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_group", type=str, required=True)
    parser.add_argument("--control_group", type=str, required=True)
    parser.add_argument("--on_target_motif", type=str, required=True)
    parser.add_argument("--on_target_location", type=str, required=True)
    parser.add_argument("--read_stats_tsv", type=str, required=True)
    parser.add_argument("--agg_site_tsv", type=str, required=True)
    parser.add_argument(
        "--peaks_dist_tsvs",
        required=True,
        nargs="+",
    )
    parser.add_argument("--ctrl_peaks_dist_tsv", type=str, required=True)
    parser.add_argument("--ctrl_wide_window_sum_dist_tsv", type=str, required=True)
    parser.add_argument("--test_pval_reporting_thresh", type=float, required=True)
    parser.add_argument("--test_pval_calling_thresh", type=float, required=True)
    parser.add_argument("--background_pval_calling_thresh", type=float, required=True)
    parser.add_argument("--pipeline_version", type=str, required=True)
    parser.add_argument("--out_path", type=str, required=True)

    args = parser.parse_args()
    create_html_report(
        args.sample_group,
        args.control_group,
        args.on_target_motif,
        args.on_target_location,
        args.read_stats_tsv,
        args.agg_site_tsv,
        args.peaks_dist_tsvs,
        args.ctrl_peaks_dist_tsv,
        args.ctrl_wide_window_sum_dist_tsv,
        args.test_pval_reporting_thresh,
        args.test_pval_calling_thresh,
        args.background_pval_calling_thresh,
        args.pipeline_version,
        args.out_path,
    )


def create_html_report(
    sample_group: str,
    control_group: str,
    on_target_motif: str,
    on_target_location: str,
    read_stats_tsv_path: str,
    agg_site_tsv_path: str,
    peaks_dist_tsvs: List[str],
    ctrl_peaks_dist_tsv: str,
    ctrl_wide_window_sum_dist_tsv: str,
    test_pval_reporting_thresh: float,
    test_pval_calling_thresh: float,
    background_pval_calling_thresh: float,
    pipeline_version: str,
    out_path: str,
) -> str:
    df_threshold, threshold_formatter_mapper = get_threshold_summary(
        test_pval_calling_thresh,
        test_pval_reporting_thresh,
        background_pval_calling_thresh,
        ctrl_peaks_dist_tsv,
        ctrl_wide_window_sum_dist_tsv,
    )
    threshold_table_str = df_threshold.to_html(
        table_id="threshold-table",
        formatters=threshold_formatter_mapper,
        index=False,
        classes="table table-bordered table-striped table-hover",
    )

    df_read_stats, stats_formatter_mapper, concs = get_read_stats_summary(
        read_stats_tsv_path, sample_group, control_group
    )
    test_concs = [c for c in concs if c != 0]
    stats_table_str = df_read_stats.to_html(
        table_id="stats-table",
        formatters=stats_formatter_mapper,
        index=False,
        classes="table table-bordered table-striped table-hover",
    )

    df_sites = pd.read_csv(agg_site_tsv_path, sep="\t", comment="#")
    df_num_called_sites, site_calling_formatter_mapper = get_site_calling_breakdown(
        df_sites, test_concs
    )

    agg_site_fields = attrs.fields(AggregateSite)
    site_fields = attrs.fields(Site)
    output_fields = [
        agg_site_fields.rank,
        site_fields.chrom,
        site_fields.position_1b,
        site_fields.closest_motif,
        site_fields.motif_location,
        site_fields.motif_score,
        site_fields.num_subs,
        site_fields.num_gaps,
        site_fields.site_called,
        site_fields.test_pval,
        site_fields.background_pval,
        site_fields.mean_pct_max_peak,
        site_fields.pct_max_peak_vals,
        site_fields.mean_signal,
        site_fields.signal_vals,
        site_fields.total_test_reads,
        site_fields.total_control_reads,
        site_fields.mean_wide_control_reads,
        site_fields.test_read_vals,
        site_fields.control_read_vals,
        site_fields.wide_control_read_vals,
        site_fields.total_signal,
        site_fields.digested,
        site_fields.high_background,
    ]
    field_keys, title_mapper, formatter_mapper = [], {}, {}
    for field in output_fields:
        keys = AggregateSite.get_field_keys(field, test_concs)
        titles = AggregateSite.get_field_titles(field, test_concs)
        for key, title in zip(keys, titles):
            field_keys.append(key)
            title_mapper[key] = title
            if "formatter" in field.metadata:
                formatter_mapper[title] = field.metadata["formatter"]

    df_sites_titles = df_sites[field_keys].rename(columns=title_mapper)
    sites_table_str = df_sites_titles.to_html(
        table_id="sites-table",
        formatters=formatter_mapper,
        escape=False,
        index=False,
        classes="table table-bordered table-striped table-hover",
    )

    num_called_sites_table_str = df_num_called_sites.to_html(
        table_id="num-called-sites-table",
        formatters=site_calling_formatter_mapper,
        index=False,
        classes="table table-bordered table-striped table-hover",
    )

    peak_signal_plot = get_peak_signal_plot([ctrl_peaks_dist_tsv] + peaks_dist_tsvs)

    html_str = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Interactive Table</title>
        <!-- Bootstrap CSS -->
        <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
        <!-- DataTables CSS -->
        <link rel="stylesheet" href="https://cdn.datatables.net/2.1.8/css/dataTables.dataTables.min.css">
        <!-- jQuery and DataTables JS -->
        <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
        <script src="https://cdn.datatables.net/2.1.8/js/dataTables.min.js"></script>

        <style>

            /*
            table.dataTable thead>tr>th.dt-orderable-asc span.dt-column-order:before,
            table.dataTable thead>tr>th.dt-ordering-asc span.dt-column-order:before,
            table.dataTable thead>tr>td.dt-orderable-asc span.dt-column-order:before,
            table.dataTable thead>tr>td.dt-ordering-asc span.dt-column-order:before {{
                bottom: 18% !important;
            }}

            table.dataTable thead>tr>th.dt-orderable-desc span.dt-column-order:after,
            table.dataTable thead>tr>th.dt-ordering-desc span.dt-column-order:after,
            table.dataTable thead>tr>td.dt-orderable-desc span.dt-column-order:after,
            table.dataTable thead>tr>td.dt-ordering-desc span.dt-column-order:after {{
                top: 82% !important;
            }}
            */

            table.dataTable thead>tr>th.dt-orderable-asc span.dt-column-order,
            table.dataTable thead>tr>th.dt-orderable-desc span.dt-column-order,
            table.dataTable thead>tr>th.dt-ordering-asc span.dt-column-order,
            table.dataTable thead>tr>th.dt-ordering-desc span.dt-column-order,
            table.dataTable thead>tr>td.dt-orderable-asc span.dt-column-order,
            table.dataTable thead>tr>td.dt-orderable-desc span.dt-column-order,
            table.dataTable thead>tr>td.dt-ordering-asc span.dt-column-order,
            table.dataTable thead>tr>td.dt-ordering-desc span.dt-column-order {{
                position: absolute;
                right: 3px;
                top: 0;
                bottom: 0;
                width: 12px;
            }}


            div.dt-container div.dt-layout-row {{
                margin: 0 !important;
            }}

            .dt-search {{
                margin-bottom: 0.25rem !important;
            }}

            table.dataTable tbody tr {{
                line-height: 1.2;
                height: 30px;
            }}

            table.dataTable tbody td {{
                padding: 4px 8px;
                white-space: nowrap;
            }}

            table.dataTable thead th {{
                padding: 6px 20px 6px 10px !important;
            }}

            .card {{
                border-color: rgba(0, 0, 0, 0.40);
                background-color: #f7fafc;
            }}

            .meta-label-cell {{
                font-weight: bold;
                white-space: nowrap;
                vertical-align: top;
                padding-right: 8px;
            }}

            div.dt-scroll-head {{
                background-color: #055081;
                border-top-left-radius: 5px;
                border-top-right-radius: 5px;
            }}

            .dt-scroll-body {{
                background: white !important;
                border: 0 !important;
            }}

            .dataTable > thead > tr {{
                background-color: #055081;
                background-clip: padding-box;
            }}

            .dataTable > thead > tr > th {{
                background-color: #055081;
                color: #ffffff;
                /*font-weight: normal !important;*/
                border-bottom: 0px !important;
                white-space: pre-wrap !important;
            }}

            .table td, .table th {{
                border: 1px solid #dee2e6;
            }}

            .dt-scroll {{
                border-radius: 4px;
                overflow: hidden;
            }}

            .dt-layout-table {{
                border: 1px #055081 solid;
                background-color: #055081;
                border-radius: 5px;
            }}

            table.table-bordered.dataTable th:last-child,
            table.table-bordered.dataTable th:last-child,
            table.table-bordered.dataTable td:last-child,
            table.table-bordered.dataTable td:last-child {{
                border-right-width: 0;
            }}

            table.dataTable.table-bordered th {{
                border-top-width: 0;
            }}

            .dataTable > thead > tr {{
                background-color: #055081;
                background-clip: padding-box;
            }}

            .dataTable tr {{
                border-width: 0 !important;
            }}

            .dataTable > thead {{
                border-bottom: 0;
            }}

            table.table-bordered.dataTable th, table.table-bordered.dataTable td {{
                border-left-width: 0;
            }}

            table.dataTable {{
                border: 0px;
                font-size: 10pt;
            }}

            .dt-length > select {{
                margin-right: 0.25rem;
            }}

            .flaggedCell {{
                background-color: rgba(220, 53, 69, 0.8) !important;
                font-weight: bold;
            }}

        </style>
    </head>
    <body>
        <div class="container-fluid" style="padding: 0.5rem 2rem">
            <div class="row mb-1">
                <div class="col">
                    <h2 class="mb-0" style="margin-left: 0.25rem">The SITE-Seq® Assay Report</h2>
                </div>
            </div>
            <div class="row mb-2">
                <div class="col-5" style="padding-right: 0.25rem !important">
                    <div class="card" style="padding: 0.5rem; height: 100%">
                        <h4 class="mb-1">Metadata</h4>
                        <table id="meta-table">
                            <div>
                                <tr>
                                    <td class="meta-label-cell">Sample Group:</td>
                                    <td>{sample_group}</td>
                                </tr>
                                <tr>
                                    <td class="meta-label-cell">Control Group:</td>
                                    <td>{control_group}</td>
                                </tr>
                                <tr>
                                    <td class="meta-label-cell">Concentrations:</td>
                                    <td>{", ".join(f"{c}nM" for c in concs)}</td>
                                </tr>
                                <tr>
                                    <td class="meta-label-cell">On-Target Motif:</td>
                                    <td>{format_dna_sequence(on_target_motif)}</td>
                                </tr>
                                <tr>
                                    <td class="meta-label-cell">On-Target Location:</td>
                                    <td>{on_target_location}</td>
                                </tr>
                                <tr>
                                    <td class="meta-label-cell">Pipeline Version:</td>
                                    <td>{pipeline_version}</td>
                                </tr>
                            </div>
                        </table>
                    </div>
                </div>
                <div class="col-7" style="padding-left: 0.25rem !important">
                    <div class="card" style="padding: 0.5rem">
                        <h4 class="mb-1" style="margin-left: 0.25rem">Called Site Summary</h4>
                        {num_called_sites_table_str}
                    </div>
                </div>
            </div>
            <div class="row mb-2">
                <div class="col">
                    <div class="card" style="padding: 0.5rem">
                        <h4 class="mb-1" style="margin-left: 0.25rem">Read & Alignments Stats</h4>
                        {stats_table_str}
                    </div>
                </div>
            </div>
            <div class="row mb-2">
                <div class="col">
                    <div class="card" style="padding: 0.5rem">
                        <h4 class="mb-1" style="margin-left: 0.25rem">Thresholds</h4>
                        {threshold_table_str}
                    </div>
                </div>
            </div>
            <div class="row mb-2">
                <div class="col">
                    <div class="card" style="padding: 0.5rem">
                        <h4 class="mb-1" style="margin-left: 0.25rem">Peak Signal Distributions</h4>
                        {peak_signal_plot.to_html(full_html=False, include_plotlyjs="cdn")}
                    </div>
                </div>
            </div>
            <div class="row mb-2">
                <div class="col">
                    <div class="card" style="padding: 0.5rem">
                        {sites_table_str}
                    </div>
                </div>
            </div>
        </div>
        <script>
            // Initialize DataTables on the table
            $(document).ready(function() {{
                $('#threshold-table').DataTable({{
                    scrollCollapse: true,
                    scrollY: false,
                    scrollX: false,
                    paging: false,
                    searching: false,
                    info: false,
                    ordering: false,
                    layout: {{
                        topStart: null,
                        topEnd: null,
                        bottomStart: null,
                        bottomEnd: null
                    }}
                }});
                $('#num-called-sites-table').DataTable({{
                    scrollCollapse: true,
                    scrollY: false,
                    scrollX: false,
                    paging: false,
                    searching: false,
                    info: false,
                    layout: {{
                        topStart: null,
                        topEnd: null
                    }},
                    rowCallback: function (row, data) {{
                        if (data[2] == 'False') {{
                            $('td:eq(2)', row).addClass('flaggedCell');
                        }}
                    }}
                }});
                $('#stats-table').DataTable({{
                    scrollCollapse: true,
                    scrollY: '75vh',
                    scrollX: true,
                    paging: false,
                    searching: false,
                    info: false,
                    layout: {{
                        topStart: null,
                        topEnd: null
                    }}
                }});
                $('#sites-table').DataTable({{
                    scrollCollapse: true,
                    scrollY: '70vh',
                    scrollX: true,
                    pageLength: 200,
                    lengthMenu: [50, 100, 200, 500, 2000],
                    layout: {{
                        topStart: () => '<h4 class="mb-0">Called/Reported Sites</h4>',
                        bottomStart: ["pageLength", "info"]
                    }}
                }});
            }});
        </script>
    </body>
    </html>
    """

    with open(out_path, "w") as fh:
        fh.write(html_str)


###########################################
##   Report component creation methods   ##
###########################################


def get_read_stats_summary(
    read_stats_tsv_path: str, sample_group: str, control_group: str
) -> Tuple[pd.DataFrame, Dict[str, Callable], List[int]]:
    title_mapper = {
        SITESeqQCColumn.IDENTIFIER: "Sample",
        SITESeqQCColumn.TOTAL_READS: "Total Reads",
        SITESeqQCColumn.LOWQ_READS: "Low-Q Reads",
        SITESeqQCColumn.TOO_MANY_N_READS: "N Reads",
        SITESeqQCColumn.SHORT_READS: "Short Reads",
        SITESeqQCColumn.LOWQ_UMIS: "Low-Q UMIs",
        SITESeqQCColumn.READS_PF: "Reads PF",
        SITESeqQCColumn.PCT_READS_PF: "% Reads PF",
        SITESeqQCColumn.PCT_Q30_READS_PF: "% Q30 (PF)",
        SITESeqQCColumn.ALIGNED_READS: "Aligned Reads",
        SITESeqQCColumn.PCT_ALIGNED: "% Aligned",
        SITESeqQCColumn.TOTAL_DEDUP_UMIS: "Total Dedup. Alignments",
        SITESeqQCColumn.PCT_DUPLICATION: "% Duplication",
        SITESeqQCColumn.CLIPPED_ALIGNMENTS: "Clipped Alignments",
        SITESeqQCColumn.LOW_MAPQ_ALIGNMENTS: "Low-Q Alignments",
        SITESeqQCColumn.ALIGNMENTS_PF: "Alignments PF",
        SITESeqQCColumn.TOTAL_UNIQUE_SITES: "Unique Sites",
    }
    formatter_mapper = {
        SITESeqQCColumn.TOTAL_READS: format_integer,
        SITESeqQCColumn.LOWQ_READS: format_integer,
        SITESeqQCColumn.TOO_MANY_N_READS: format_integer,
        SITESeqQCColumn.SHORT_READS: format_integer,
        SITESeqQCColumn.READS_PF: format_integer,
        SITESeqQCColumn.PCT_READS_PF: format_percentage,
        SITESeqQCColumn.PCT_Q30_READS_PF: format_percentage,
        SITESeqQCColumn.LOWQ_UMIS: format_integer,
        SITESeqQCColumn.ALIGNED_READS: format_integer,
        SITESeqQCColumn.PCT_ALIGNED: format_percentage,
        SITESeqQCColumn.TOTAL_DEDUP_UMIS: format_integer,
        SITESeqQCColumn.PCT_DUPLICATION: format_percentage,
        SITESeqQCColumn.ALIGNMENTS_PF: format_integer,
        SITESeqQCColumn.CLIPPED_ALIGNMENTS: format_integer,
        SITESeqQCColumn.LOW_MAPQ_ALIGNMENTS: format_integer,
        SITESeqQCColumn.TOTAL_UNIQUE_SITES: format_integer,
    }
    title_formatter_mapper = {
        title_mapper[key]: formatter for key, formatter in formatter_mapper.items()
    }
    df = pd.read_csv(read_stats_tsv_path, sep="\t")
    mask = df.sample_group.isin((sample_group, control_group))
    sub_df_titled = df[mask][list(title_mapper)].rename(columns=title_mapper)
    return (
        sub_df_titled,
        title_formatter_mapper,
        sorted(df[mask][SITESeqQCColumn.CONCENTRATION].unique()),
    )


def get_site_calling_breakdown(
    df_sites: pd.DataFrame, concs: List[int]
) -> Tuple[pd.DataFrame, Dict[str, Callable]]:
    rows = []
    site_fields = attrs.fields(Site)
    is_on_target_key = site_fields.is_on_target.name
    site_called_keys = AggregateSite.get_field_keys(site_fields.site_called, concs)
    mean_signal_keys = AggregateSite.get_field_keys(site_fields.mean_signal, concs)
    for conc, site_called_key, mean_signal_key in zip(
        concs, site_called_keys, mean_signal_keys
    ):
        on_target_called = (
            df_sites[is_on_target_key] & (df_sites[site_called_key])
        ).any()
        max_called_signal = df_sites[mean_signal_key][df_sites[site_called_key]].max()
        rows.append(
            {
                "Concentration": conc,
                "Sites Called": df_sites[site_called_key].sum(),
                "On-Target Called": on_target_called,
                "Max Called Signal": (
                    max_called_signal if not pd.isnull(max_called_signal) else None
                ),
            }
        )
    formatter_mapper = {
        "Sites Called": format_integer,
        "Max Called Signal": format_decimal,
    }
    return pd.DataFrame(rows), formatter_mapper


def get_threshold_summary(
    peak_calling_thresh: float,
    peak_reporting_thresh: float,
    background_filter_thresh: float,
    control_peak_dist_path: str,
    control_wide_sums_dist_path: str,
) -> Tuple[pd.DataFrame, Dict[str, Callable]]:
    site_calculator = SiteCalculator(
        control_peak_dist_path,
        control_wide_sums_dist_path,
        peak_calling_thresh,
        background_filter_thresh,
        peak_reporting_thresh,
    )
    site_fields = attrs.fields(Site)
    rows = [
        {
            "Name": "Calling threshold",
            "Stat": site_fields.mean_signal.metadata["title"],
            "Threshold": peak_calling_thresh,
            "Thresh. Value": site_calculator.test_pval_calc.get_threshold_for_pvalue(
                peak_calling_thresh
            ),
        },
        {
            "Name": "Reporting threshold",
            "Stat": site_fields.mean_signal.metadata["title"],
            "Threshold": peak_reporting_thresh,
            "Thresh. Value": site_calculator.test_pval_calc.get_threshold_for_pvalue(
                peak_reporting_thresh
            ),
        },
        {
            "Name": "Background filter threshold",
            "Stat": site_fields.mean_wide_control_reads.metadata["title"],
            "Threshold": background_filter_thresh,
            "Thresh. Value": site_calculator.background_pval_calc.get_threshold_for_pvalue(
                peak_reporting_thresh
            ),
        },
    ]
    formatter_mapper = {
        "Threshold": format_pvalue,
        "Thresh. Value": format_decimal,
    }
    return pd.DataFrame(rows), formatter_mapper


def get_peak_signal_plot(peaks_dist_tsvs: List[str]) -> go.Figure:
    num_sites_col = "Num. Sites"
    conc_col = "Conc."
    stat_col = attrs.fields(Site).mean_signal.metadata["title"]
    rows = []
    for peaks_dist_tsv in peaks_dist_tsvs:
        conc = f"{int(peaks_dist_tsv.split('nM')[0].split('_')[-1])}nM"
        stat_map = SiteCalculator.parse_dist_tsv(peaks_dist_tsv)
        rows.extend(
            {conc_col: conc, num_sites_col: count, stat_col: stat}
            for stat, count in stat_map.items()
        )

    df = pd.DataFrame(rows)
    bin_width = 1
    nbins = min(int(df[stat_col].max() / bin_width) * 2, 2000)
    fig = px.histogram(
        df,
        x=stat_col,
        y=num_sites_col,
        color=conc_col,
        log_y=True,
        nbins=nbins,
        barmode="overlay",
        marginal="rug",
        template=PLOTLY_THEME,
        title=None,
        height=300,
    )
    updatemenus = [
        dict(
            type="buttons",
            direction="down",
            buttons=list(
                [
                    dict(
                        args=[{"yaxis.type": "log"}],
                        label="Log Scale",
                        method="relayout",
                    ),
                    dict(
                        args=[{"yaxis.type": "linear"}],
                        label="Linear Scale",
                        method="relayout",
                    ),
                ]
            ),
        ),
    ]
    fig.update_layout(
        yaxis_title=num_sites_col,
        updatemenus=updatemenus,
        margin=dict(l=20, r=20, t=20, b=20),
        paper_bgcolor="rgba(0, 0, 0, 0)",
    )
    return fig


if __name__ == "__main__":
    main()
