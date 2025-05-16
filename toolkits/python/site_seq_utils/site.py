#########################################################################
#########################################################################
##                   Classes for representing sites                    ##
#########################################################################
#########################################################################

from __future__ import annotations
import attrs
import re
import numpy as np
import pandas as pd

from typing import List, Dict, Any, Union, Tuple, Optional
from functools import cached_property

from .formatters import (
    format_dna_sequence,
    format_coordinate,
    format_number_list,
    format_percentage_list,
    format_integer,
    format_decimal,
    format_pvalue,
    format_percentage,
)


def float_list_tuple_converter(
    value: Union[List[float], str, None],
) -> Optional[List[float]]:
    if value is None:
        return None
    elif isinstance(value, tuple):
        return value
    elif isinstance(value, list):
        return tuple(value)
    elif isinstance(value, str):
        eval_result = eval(value)
        if not isinstance(eval_result, (list, tuple)):
            raise ValueError(
                f'Did not get a list or tuple when evaluating "{value}" (got {type(eval_result)})'
            )
        return tuple(eval_result)
    else:
        raise ValueError(f"Expected a tuple, list, or str, got {type(value)}")


def list_to_tuple_converter(value: Optional[List[Any]]) -> Optional[Tuple[Any]]:
    if value is None:
        return None
    elif isinstance(value, list):
        return tuple(value)
    else:
        raise ValueError(f"Expected a list, got {tuple(value)}")


@attrs.define(frozen=True)
class Site:
    """
    Class for representing a genomic site that may or may not have evidence of biochemical
    digestion. Stores genomic coordinate, sequence motif information, and various measurements
    taken from test and control samples at that site. Measurements are taken from one or more
    replicates of a single RNP concentration.

    ...
    Attributes
    ----------
    chrom : str
        The chromosome of the site in the reference genome
    position : int
        The zero-based position of the site in the reference genome
    position_1b : int
        The one-based position of the site in the reference genome
    is_on_target : bool
        Whether the site corresponds to the on-target location
    site_called : bool
        Whether the site has been "called", i.e. considered to be the product of true
        CRISPR-induced biochemical digestion
    digested : bool
        Whether the site has evidence of biochemical digestion before accounting for regions
        with high background digestion
    high_background : bool
        Whether the site is in a region with high background digestion
    test_pval : float
        The p-value of the statistical test comparing the signal at this site to the control
        vs control signal distribution
    background_pval : float
        The p-value of the statistical test comparing the background digestion to the
        distribution of background digestion
    signal_vals : Tuple[float]
        The signal measurements in the test replicate(s)
    test_read_vals : Tuple[float]
        The number of reads measured in each test replicate
    control_read_vals : Tuple[float]
        The number of reads measured in each control replicate
    wide_control_read_vals : Tuple[float]
        The number of reads measured in each control replicate within a wider window
    rank : Optional[int]
        The rank of the site among all other sites in the reference genome
    motif_location : Optional[str]
        The location of the sequence motif best-matching the on-target motif within close
        proximity of the site
    closest_motif : Optional[str]
        The motif sequence best-matching the on-target motif within close proximity of the
        site
    motif_score : Optional[int]
        A score representing the similarity between closest_motif and the on-target motif
    num_subs : Optional[int]
        The number of mismatched bases in the alignment between closest_motif and the on-target
        motif
    num_gaps : Optional[int]
        The number of gaps in the alignment between closest_motif and the on-target motif
    pct_max_peak_vals : Optional[Tuple[Optional[float]]]
        The signal strength of the site as a percentage of the maximum signal strength measured
        in all called sites in each test replicate
    total_signal : float
        The sum of the signal measurements in the test replicates
    mean_signal : float
        The mean of the signal measurements in the test replicates
    total_test_reads : float
        The total number of reads measured in the test replicates
    total_control_reads : float
        The total number of reads measured in the control replicates
    total_wide_control_reads : float
        The total number of reads measured in the control replicates within a wider winder
    mean_wide_control_reads : float
        The mean of the number of reads measured in the control replicates within a wider window
    mean_pct_max_peak : Optional[float]
        The mean of the signal strength of the site as a percentage of the maximum signal strength
        measured in all called sites in each test replicate
    comparison_tuple : Tuple
        A tuple representing the "value" of this site. Used to compare sites

    Methods
    -------
    load(cls, Union[Dict, pd.Series]) -> Site
        Create a Site object from a dictionary or pd.Series object

    """

    chrom: str = attrs.field(
        metadata=dict(
            title="Chrom",
            data_field=False,
            control_data=False,
            motif_field=False,
        )
    )
    position: int = attrs.field(
        metadata=dict(
            title="Position",
            data_field=False,
            control_data=False,
            motif_field=False,
            formatter=format_integer,
        )
    )
    position_1b: int = attrs.field(
        metadata=dict(
            title="Position (1-based)",
            data_field=False,
            control_data=False,
            motif_field=False,
            formatter=format_integer,
        )
    )
    is_on_target: bool = attrs.field(
        metadata=dict(
            title="Is On-Target",
            data_field=False,
            control_data=False,
            motif_field=False,
        )
    )
    site_called: bool = attrs.field(
        metadata=dict(
            title="Called",
            data_field=True,
            control_data=False,
            motif_field=False,
        )
    )
    digested: bool = attrs.field(
        metadata=dict(
            title="Digested",
            data_field=True,
            control_data=False,
            motif_field=False,
        )
    )
    high_background: bool = attrs.field(
        metadata=dict(
            title="High-background",
            data_field=True,
            control_data=True,
            motif_field=False,
        )
    )
    test_pval: float = attrs.field(
        metadata=dict(
            title="Test p-value",
            data_field=True,
            control_data=False,
            motif_field=False,
            formatter=format_pvalue,
        )
    )
    background_pval: float = attrs.field(
        metadata=dict(
            title="Background p-value",
            data_field=True,
            control_data=True,
            motif_field=False,
            formatter=format_pvalue,
        )
    )
    signal_vals: Tuple[float] = attrs.field(
        converter=float_list_tuple_converter,
        metadata=dict(
            title="Signals",
            data_field=True,
            control_data=False,
            motif_field=False,
            formatter=format_number_list,
        ),
    )
    test_read_vals: Tuple[float] = attrs.field(
        converter=float_list_tuple_converter,
        metadata=dict(
            title="Reads",
            data_field=True,
            control_data=False,
            motif_field=False,
            formatter=format_number_list,
        ),
    )
    control_read_vals: Tuple[float] = attrs.field(
        converter=float_list_tuple_converter,
        metadata=dict(
            title="Control Reads",
            data_field=True,
            control_data=True,
            motif_field=False,
            formatter=format_number_list,
        ),
    )
    wide_control_read_vals: Tuple[float] = attrs.field(
        converter=float_list_tuple_converter,
        metadata=dict(
            title="Control Reads (Wide)",
            data_field=True,
            control_data=True,
            motif_field=False,
            formatter=format_number_list,
        ),
    )

    rank: Optional[int] = attrs.field(
        default=None,
        metadata=dict(
            title="Rank",
            data_field=True,
            control_data=False,
            motif_field=False,
            formatter=format_integer,
        ),
    )
    motif_location: Optional[str] = attrs.field(
        default=None,
        metadata=dict(
            title="Motif Loc.",
            data_field=False,
            control_data=False,
            motif_field=True,
            formatter=format_coordinate,
        ),
    )
    closest_motif: Optional[str] = attrs.field(
        default=None,
        metadata=dict(
            title="Closest Motif",
            data_field=False,
            control_data=False,
            motif_field=True,
            formatter=format_dna_sequence,
        ),
    )
    motif_score: Optional[int] = attrs.field(
        default=None,
        metadata=dict(
            title="Motif Score",
            data_field=False,
            control_data=False,
            motif_field=True,
            formatter=format_integer,
        ),
    )
    num_subs: Optional[int] = attrs.field(
        default=None,
        metadata=dict(
            title="Num. Subs",
            data_field=False,
            control_data=False,
            motif_field=True,
            formatter=format_integer,
        ),
    )
    num_gaps: Optional[int] = attrs.field(
        default=None,
        metadata=dict(
            title="Num. Gaps",
            data_field=False,
            control_data=False,
            motif_field=True,
            formatter=format_integer,
        ),
    )
    pct_max_peak_vals: Optional[Tuple[Optional[float]]] = attrs.field(
        converter=float_list_tuple_converter,
        default=None,
        metadata=dict(
            title="% of Max Peak",
            data_field=True,
            control_data=False,
            motif_field=False,
            formatter=format_percentage_list,
        ),
    )

    # Derived fields
    total_signal: float = attrs.field(
        init=False,
        metadata=dict(
            title="Total Signal",
            data_field=True,
            control_data=False,
            motif_field=False,
            formatter=format_decimal,
        ),
    )
    mean_signal: float = attrs.field(
        init=False,
        metadata=dict(
            title="Mean Signal",
            data_field=True,
            control_data=False,
            motif_field=False,
            formatter=format_decimal,
        ),
    )
    total_test_reads: float = attrs.field(
        init=False,
        metadata=dict(
            title="Total Reads",
            data_field=True,
            control_data=False,
            motif_field=False,
            formatter=format_decimal,
        ),
    )
    total_control_reads: float = attrs.field(
        init=False,
        metadata=dict(
            title="Total Control Reads",
            data_field=True,
            control_data=True,
            motif_field=False,
            formatter=format_decimal,
        ),
    )
    total_wide_control_reads: float = attrs.field(
        init=False,
        metadata=dict(
            title="Total Control Reads (Wide)",
            data_field=True,
            control_data=True,
            motif_field=False,
            formatter=format_decimal,
        ),
    )
    mean_wide_control_reads: float = attrs.field(
        init=False,
        metadata=dict(
            title="Mean Control Reads (Wide)",
            data_field=True,
            control_data=True,
            motif_field=False,
            formatter=format_decimal,
        ),
    )
    mean_pct_max_peak: Optional[float] = attrs.field(
        init=False,
        metadata=dict(
            title="Mean % of Max Peak",
            data_field=True,
            control_data=False,
            motif_field=False,
            formatter=format_percentage,
        ),
    )

    #####################################
    # Derived-field calculation methods #
    #####################################

    @total_signal.default
    def _total_signal(self) -> float:
        return sum(self.signal_vals)

    @mean_signal.default
    def _mean_signal(self) -> float:
        return np.mean(self.signal_vals)

    @total_test_reads.default
    def _total_test_reads(self) -> float:
        return sum(self.test_read_vals)

    @total_control_reads.default
    def _total_control_reads(self) -> float:
        return sum(self.control_read_vals)

    @total_wide_control_reads.default
    def _total_wide_control_reads(self) -> float:
        return sum(self.control_read_vals)

    @mean_wide_control_reads.default
    def _mean_wide_control_reads(self) -> float:
        return np.mean(self.control_read_vals)

    @mean_pct_max_peak.default
    def _mean_pct_max_peak(self) -> Optional[float]:
        return (
            np.mean(self.pct_max_peak_vals)
            if self.pct_max_peak_vals
            and all(v is not None for v in self.pct_max_peak_vals)
            else None
        )

    @classmethod
    def load(cls, record: Union[Dict, pd.Series]) -> Site:
        return cls(
            **{
                field.name: getattr(record, field.name)
                for field in attrs.fields(cls)
                if field.init
            }
        )

    @property
    def comparison_tuple(self) -> Tuple:
        return (
            self.is_on_target,
            self.site_called,
            self.mean_signal,
        )

    def __lt__(self, other: Site) -> bool:
        return self.comparison_tuple < other.comparison_tuple


@attrs.define(frozen=True)
class AggregateSite:
    """
    Class for representing a collection of Site objects obtained from multiple RNP concentrations.

    Attributes
    ----------
    site_dict : Dict[int, Site]
        A dictionary mapping RNP concentrations to the Site object obtained at that concentration
    rank : Optional[int]
        The rank of the site relative to all other sites measured across all RNP concentrations
    total_signal : float
        The total signal of the site across all test replicates at all RNP concentrations

    Methods
    -------
    to_dict(self) -> Dict
        Helper method to convert this object to a dictionary
    load(cls, Union[Dict, pd.Series]) -> Site
        Create an AggregateSite object from a dictionary or pd.Series object
    is_on_target(self) -> bool
        Whether this site is the on-target
    comparison_tuple(self) -> Tuple
        A tuple representing the "value" of this site. Used to compare aggregate sites
    get_field_keys(cls, attrs.Attibute, List[int]) -> List[str]
        Return the keys for the given Site class field for the given RNP concentrations. The
        concentrations are used to prefix the keys
    get_field_titles(cls, attrs.Attibute, List[int]) -> List[str]
        Return the titles for the given Site class field for the given RNP concentrations. The
        concentrations are used to suffix the keys
    get_field_by_key(cls, str) -> attrs.Attribute
        Return the field of the Site class corresponding to the given key
    """

    site_dict: Dict[int, Site] = attrs.field()
    rank: Optional[int] = attrs.field(
        default=None,
        metadata=dict(
            title="Rank",
            formatter=format_integer,
        ),
    )
    total_signal: float = attrs.field(
        init=False,
        metadata=dict(
            title="Total Signal",
            formatter=format_decimal,
        ),
    )

    _conc_field_pattern = re.compile("^(?P<conc>\d+)nM_(?P<field>\S+)$")

    @total_signal.default
    def _total_signal(self) -> float:
        return sum(site.total_signal for site in self.site_dict.values())

    @site_dict.validator
    def validate_site_dict_format(self, attribute: Any, value: Any) -> None:
        if not isinstance(value, dict):
            raise ValueError(f"{attribute.name} is not a dict")
        elif len(value) == 0:
            raise ValueError(f"{attribute.name} is empty")

        for conc, site in value.items():
            if not isinstance(conc, int):
                raise ValueError(
                    f'{attribute.name} has a key "{conc}" that is not an int (type {type(conc)})'
                )
            elif not isinstance(site, Site):
                raise ValueError(
                    f'{attribute.name} has a value "{conc}" that is not a Site (type {type(site)})'
                )

    @site_dict.validator
    def validate_common_control_data(self, attribute: Any, value: Any) -> None:
        control_data_fields = [
            field for field in attrs.fields(Site) if field.metadata["control_data"]
        ]
        control_data_vals = {}
        for site in value.values():
            for field in control_data_fields:
                control_data_vals.setdefault(field.name, set()).add(
                    getattr(site, field.name)
                )

        for field_name, field_vals in control_data_vals.items():
            if len(field_vals) != 1:
                raise ValueError(
                    f'{attribute.name} has multiple values for control data field "{field_name}": '
                    f"{', '.join(map(str, field_vals))}"
                )

    @site_dict.validator
    def validate_common_non_data_fields(self, attribute: Any, value: Any) -> None:
        non_data_fields = [
            field for field in attrs.fields(Site) if not field.metadata["data_field"]
        ]
        non_data_vals = {}
        for site in value.values():
            for field in non_data_fields:
                field_value = getattr(site, field.name)
                if field_value is not None:
                    non_data_vals.setdefault(field.name, set()).add(
                        getattr(site, field.name)
                    )

        for field_name, field_vals in non_data_vals.items():
            if len(field_vals) > 1:
                raise ValueError(
                    f'{attribute.name} has multiple values for control data field "{field_name}": '
                    f"{', '.join(map(str, field_vals))}"
                )

    @site_dict.validator
    def validate_non_null_motif_fields(self, attribute: Any, value: Any) -> None:
        motif_fields = [
            field for field in attrs.fields(Site) if field.metadata["motif_field"]
        ]
        for field in motif_fields:
            if not any(
                getattr(site, field.name) is not None for site in value.values()
            ):
                raise ValueError(
                    f'{attribute.name} has all-None values for motif field "{field.name}": '
                )

    def to_dict(self) -> Dict[str, Any]:
        data = attrs.asdict(self, filter=attrs.filters.exclude(AggregateSite.site_dict))
        for field in attrs.fields(Site):
            if field.metadata["data_field"] and not field.metadata["control_data"]:
                # Test data field, add for all concentrations, prefix key with concentration
                for conc, site in self.site_dict.items():
                    data[f"{conc}nM_{field.name}"] = getattr(site, field.name)
            else:
                # Find the first non-null value
                value = None
                for site in self.site_dict.values():
                    value = getattr(site, field.name)
                    if value is not None:
                        break
                data[field.name] = value
        return data

    @classmethod
    def load(cls, record: Union[Dict, pd.Series]) -> Site:
        agg_data, conc_site_data_map = {}, {}
        common_site_keys = []
        record_data = record.to_dict() if isinstance(record.pd.Series) else record
        for key, value in record_data.items():
            conc_field_match = re.match(cls._conc_field_pattern, key)
            if conc_field_match:
                conc = int(conc_field_match.group("conc"))
                site_field = conc_field_match.group("field")
                if not hasattr(Site, site_field):
                    raise ValueError(
                        f"Unrecognized concentration-prefixed field: {key}"
                    )
                conc_site_data_map.setdefault(conc, {})[site_field] = value
            elif hasattr(cls, key):
                agg_data[key] = value
            elif hasattr(Site, key):
                common_site_keys.append(key)
            else:
                raise ValueError(f"Unrecognized field: {key}")

        for key in common_site_keys:
            for conc, conc_site_data in conc_site_data_map.items():
                conc_site_data[key] = record_data[key]

        site_dict = {
            conc: Site.load(conc_site_data)
            for conc, conc_site_data in conc_site_data_map.items()
        }
        return cls(
            site_dict=site_dict,
            **agg_data,
        )

    @cached_property
    def is_on_target(self) -> bool:
        for site in self.site_dict.values():
            if site.is_on_target is not None:
                return site.is_on_target
        # Should be caught by validator in theory
        raise ValueError(f"No non-null value for {Site.is_on_target.name}")

    @cached_property
    def comparison_tuple(self) -> Tuple:
        min_all_called_conc = float("inf")
        for conc in sorted(self.site_dict, reverse=True):
            if self.site_dict[conc].site_called:
                min_all_called_conc = conc
            else:
                break

        return (
            self.is_on_target,
            *(
                (
                    site.mean_signal
                    if site.site_called and conc >= min_all_called_conc
                    else 0
                )
                for conc, site in sorted(self.site_dict.items(), key=lambda t: t[0])
            ),
            self.total_signal,
        )

    def __lt__(self, other: Site) -> bool:
        return self.comparison_tuple < other.comparison_tuple

    @classmethod
    def get_field_keys(cls, field: attrs.Attibute, concs: List[int]) -> List[str]:
        """
        Return the keys for the given Site class field for the given RNP concentrations. The
        concentrations are used to prefix the keys
        """
        if not field.metadata.get("data_field") or field.metadata.get("control_data"):
            return [field.name]
        else:
            return [f"{conc}nM_{field.name}" for conc in concs]

    @classmethod
    def get_field_titles(cls, field: attrs.Attibute, concs: List[int]) -> List[str]:
        """
        Return the titles for the given Site class field for the given RNP concentrations. The
        concentrations are used to suffix the keys
        """
        if not field.metadata.get("data_field") or field.metadata.get("control_data"):
            return [field.metadata["title"]]
        else:
            return [f"{field.metadata['title']} ({conc}nM)" for conc in concs]

    @classmethod
    def get_field_by_key(cls, key: str) -> attrs.Attribute:
        """
        Return the field of the Site class corresponding to the given key
        """
        conc_field_match = re.match(cls._conc_field_pattern, key)
        if conc_field_match:
            field_key = conc_field_match.group("field")
            if not hasattr(Site, field_key):
                raise ValueError(f"Unrecognized concentration-prefixed field: {key}")
            return getattr(attrs.fields(Site), field_key)
        elif hasattr(cls, key):
            return getattr(attrs.fields(cls), key)
        elif hasattr(Site, key):
            return getattr(attrs.fields(Site), key)
        else:
            raise ValueError(f"Unrecognized field: {key}")
