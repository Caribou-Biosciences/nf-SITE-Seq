#########################################################################
#########################################################################
##           Class for calculating p-values based on eCDF              ##
#########################################################################
#########################################################################


from typing import Union, Dict
import bisect
import numpy as np

Number = Union[int, float]


class PValCalculator:
    """
    Class for calculating p-values from an eCDF with linear interpolation
    """

    def __init__(
        self,
        stat_count_map: Dict[Number, int],
        precision: int = 6,
        pval_precision: int = 16,
        min_p: bool = True,
    ) -> None:
        self.stat_count_map = stat_count_map
        self.stat_vals = sorted(stat_count_map)
        self.stat_counts = [stat_count_map[v] for v in self.stat_vals]
        self.precision = precision
        self.pval_precision = pval_precision
        self._min_val = min(stat_count_map)
        self._max_val = max(stat_count_map)
        self._total_count = sum(stat_count_map.values())
        self._min_pval = (
            round(float(1 / self._total_count), self.pval_precision) if min_p else 0
        )
        self._max_pval = 1.0
        self._memo = {}

    def get_threshold_for_pvalue(self, pval: float) -> Number:
        """Return the test value that would produce the given p-value"""
        pval = round(pval, self.pval_precision)
        if pval <= self._min_pval:
            return self._max_val

        target_count = pval * self._total_count
        cumsum = np.cumsum(self.stat_counts[::-1])[::-1]
        target_gte_index = np.argwhere(cumsum >= target_count).max()
        if target_gte_index == len(self.stat_vals) - 1:
            return self._max_val

        gte_total = cumsum[target_gte_index]
        lt_total = cumsum[target_gte_index + 1]
        gte_stat = self.stat_vals[target_gte_index]
        lt_stat = self.stat_vals[target_gte_index + 1]

        threshold = gte_stat + (gte_total - target_count) / (gte_total - lt_total) * (
            lt_stat - gte_stat
        )
        return threshold

    def p_value(self, test_stat: Number) -> float:
        """
        Calculate p-value for the given test statistic value as the sum of the PDF
        for values greater than or equal to the test value. If the test value lies
        between values in the null distribution, linearly interpolate between the
        surrounding values.
        """
        test_stat = float(test_stat)
        if test_stat in self._memo:
            return self._memo[test_stat]

        test_stat = round(test_stat, self.precision)

        # Return 1.0 if given value is less than tne minimum of the null distribution
        if test_stat <= self._min_val:
            p_val = self._max_pval
        # Return 0 or 1 / n if given value is larger than the maximum of the null distribution
        elif test_stat > self._max_val:
            p_val = self._min_pval
        # If the given value exists in null distribution, return p value by summing of the counts
        # within the range from given value to the right of the null distribution.
        #                         [      p value     ]
        #                         ------------------->
        # ------------------------|-------------------
        # 0  1  3  4  5  6  7  8  9  10  11  12  13  14
        #                         v
        elif test_stat in self.stat_count_map:
            total_gte_count = sum(
                count for val, count in self.stat_count_map.items() if val >= test_stat
            )
            p_val = total_gte_count / self._total_count
        # If the given value does not exist in the null distribution,
        # sum the counts from the next larger value to the right, and
        # then get the addition counts from next less value.
        #                       [       p value      ]
        #                       -->------------------>
        # ----------------------|-|-------------------
        # 0  1  3  4  5  6  7  8  9  10  11  12  13  14
        #                       v
        else:
            index = bisect.bisect_left(self.stat_vals, test_stat)
            lt_val = self.stat_vals[index - 1]
            gt_val = self.stat_vals[index]
            lt_count = self.stat_counts[index - 1]
            ratio = abs(gt_val - test_stat) / abs(gt_val - lt_val)
            total_gte_count = sum(self.stat_counts[index:])
            interp_count = ratio * lt_count
            p_val = (total_gte_count + interp_count) / self._total_count

        if p_val < 0 or p_val > 1.0:
            raise RuntimeError(p_val)

        p_val = round(p_val, self.pval_precision)
        self._memo[test_stat] = p_val
        return p_val

    def __call__(self, test_stat: Number) -> float:
        return self.p_value(test_stat)
