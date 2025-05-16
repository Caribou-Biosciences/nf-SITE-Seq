#########################################################################
#########################################################################
##      Class for calculating site values based on signal inputs       ##
#########################################################################
#########################################################################

from typing import Union, Dict, Optional, List
import pandas as pd
from .site import Site
from .pval_calculator import PValCalculator
from .motif_finder import MotifFinder
from .coordinate import Onedinate

Number = Union[int, float]


class SiteCalculator:
    """Class for calculating site values based on signal inputs"""

    def __init__(
        self,
        control_peak_dist_path: str,
        wide_control_peak_dist_path: str,
        peak_calling_thresh: float,
        background_filter_thresh: float,
        peak_reporting_thresh: float,
        include_background_peaks: bool = False,
        terminate_early: bool = False,
        motif_finder: Optional[MotifFinder] = None,
        on_target_loc: Optional[Onedinate] = None,
    ) -> None:
        self.test_pval_calc = SiteCalculator.get_pval_calculator(control_peak_dist_path)
        self.background_pval_calc = SiteCalculator.get_pval_calculator(
            wide_control_peak_dist_path
        )
        self.peak_calling_thresh = peak_calling_thresh
        self.background_filter_thresh = background_filter_thresh
        self.peak_reporting_thresh = peak_reporting_thresh
        self.include_background_peaks = include_background_peaks
        self.terminate_early = terminate_early
        self.motif_finder = motif_finder
        self.on_target_loc = on_target_loc
        if self.motif_finder is not None and self.on_target_loc is None:
            raise ValueError(
                "on_target_loc must be provided when motif_finder is provided"
            )

    @staticmethod
    def parse_dist_tsv(dist_tsv_path: str) -> Dict[Number, int]:
        """Parse a TSV representing a distribution of statistics"""
        df = pd.read_csv(dist_tsv_path, sep="\t", comment="#")
        return {row["mean"]: int(row["count"]) for i, row in df.iterrows()}

    @staticmethod
    def get_pval_calculator(dist_tsv_path: str) -> PValCalculator:
        """
        Return a PValCalculator objects that calculates p-values from the given distribution file
        """
        dist = SiteCalculator.parse_dist_tsv(dist_tsv_path)
        return PValCalculator(dist)

    def calculate_site(
        self,
        chrom: str,
        position: int,
        peak_stat_vals: List[Number],
        test_window_sum_vals: List[Number],
        control_window_sum_vals: List[Number],
        control_wide_window_sum_vals: List[Number],
    ) -> Optional[Dict]:
        """
        Return a dictionary of site values for the given genomic location with the given signal
        values
        """
        test_stat_mean = sum(peak_stat_vals) / len(peak_stat_vals)
        test_pval = self.test_pval_calc(test_stat_mean)
        if self.terminate_early and test_pval > self.peak_reporting_thresh:
            return None

        mean_control_wide_window = sum(control_wide_window_sum_vals) / len(
            control_wide_window_sum_vals
        )
        background_pval = self.background_pval_calc(mean_control_wide_window)
        if (
            self.terminate_early
            and not self.include_background_peaks
            and background_pval < self.background_filter_thresh
        ):
            return None

        digested = test_pval <= self.peak_calling_thresh
        high_background = background_pval <= self.background_filter_thresh
        site_called = digested and not high_background
        if self.motif_finder is not None:
            closest_motif, motif_location, motif_score, num_subs, num_gaps = (
                self.motif_finder.find_motif(chrom, position)
            )
            on_target_inter = Onedinate(motif_location).intersection(
                self.on_target_loc, ignore_strand=True
            )
            is_on_target = (
                bool(on_target_inter)
                and on_target_inter.length / self.on_target_loc.length > 0.9
            )
        else:
            closest_motif = motif_location = motif_score = num_subs = num_gaps = (
                is_on_target
            ) = None

        return Site(
            chrom=chrom,
            position=position,
            position_1b=position + 1,
            is_on_target=is_on_target,
            motif_location=motif_location,
            closest_motif=closest_motif,
            motif_score=motif_score,
            num_subs=num_subs,
            num_gaps=num_gaps,
            site_called=site_called,
            digested=digested,
            high_background=high_background,
            test_pval=test_pval,
            background_pval=background_pval,
            signal_vals=peak_stat_vals,
            test_read_vals=test_window_sum_vals,
            control_read_vals=control_window_sum_vals,
            wide_control_read_vals=control_wide_window_sum_vals,
        )
