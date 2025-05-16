#########################################################################
#########################################################################
##                Class for parsing distribution TSVs                  ##
#########################################################################
#########################################################################

import pysam
import gzip
from typing import Generator, Dict, List, Optional


class DistTSVParser:
    """
    Class for parsing/querying distribution TSVs (potentially bgzipped and Tabix indexed)
    """

    def __init__(self, tsv_path: str) -> None:
        self.tsv_path = tsv_path
        try:
            self.tbq = pysam.TabixFile(tsv_path)
            self.indexed = True
        except OSError as err:
            err_str = str(err)
            if "index" in err_str and "not found" in err_str:
                self.indexed = False
            else:
                raise
        self._columns = None
        self._memo = {}

    @property
    def columns(self) -> List[str]:
        """Return column names in this TSV"""
        if self._columns is None:
            with gzip.open(self.tsv_path, "rt") as fh:
                self._columns = fh.readline().strip().split("\t")
        return self._columns

    @property
    def sample_names(self) -> List[str]:
        """
        Return names of samples described in this TSV. Assumed to the the column names after the
        chromosome, start, and end columns
        """
        return self.columns[3:]

    @property
    def num_samples(self) -> int:
        """
        Return the number of samples (aka measurements) stored in this TSV
        """
        return len(self.columns) - 3

    def get_values_at_position(
        self, chrom: str, position: int, default_val: float = 0.0
    ) -> List[float]:
        """
        Return the values at the given position, or the given default value if no record is found
        """
        key = (chrom, position)
        if key in self._memo:
            return self._memo[key]

        for record in self.fetch(chrom, position, position + 1):
            values = [record[name] for name in self.sample_names]
            self._memo[key] = values
            if record["start"] <= position < record["end"]:
                return values

        default_vals = [default_val] * self.num_samples
        self._memo[key] = default_vals
        return default_vals

    def fetch(
        self,
        chrom: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
    ) -> Generator[Dict, None, None]:
        """
        Fetch all of the records stored at the given genomic coordinate range
        """
        if self.indexed:
            try:
                for record_str in self.tbq.fetch(chrom, start, end):
                    yield self.parse_record_str(record_str)
            except ValueError as err:
                if "could not create iterator for region" in str(err):
                    return
                else:
                    raise
        else:
            if any(v is not None for v in (chrom, start, end)):
                raise Exception("file is not indexed, cannot query by region")
            with gzip.open(self.tsv_path, "rt") as fh:
                for idx, line in enumerate(fh.readlines()):
                    if idx == 0:
                        continue
                    yield self.parse_record_str(line.strip())

    def parse_record_str(self, record_str: str) -> Dict:
        """
        Parse the given record string and convert it to a dictionary
        """
        record_vals = record_str.split("\t")
        data = {
            "chrom": record_vals[0],
            "start": int(record_vals[1]),
            "end": int(record_vals[2]),
        }
        for idx in range(3, len(record_vals)):
            data[self.columns[idx]] = float(record_vals[idx])
        return data
