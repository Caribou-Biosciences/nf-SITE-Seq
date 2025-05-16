#########################################################################
#########################################################################
##           Constant values used in SITE-Seq assay analysis           ##
#########################################################################
#########################################################################

from enum import Enum


class SITESeqQCColumn(str, Enum):
    TOTAL_READS = "total_reads"
    LOWQ_READS = "lowq_reads"
    TOO_MANY_N_READS = "too_many_n_reads"
    SHORT_READS = "short_reads"
    LOWQ_UMIS = "lowq_umis"
    HIGHQ_UMIS = "highq_umis"
    READS_PF = "reads_passing_filter"
    PCT_READS_PF = "pct_reads_passing_filter"
    PCT_Q30_READS_PF = "pct_q30_reads_pf"
    ALIGNED_READS = "aligned_reads"
    PCT_ALIGNED = "pct_aligned"
    TOTAL_DEDUP_UMIS = "total_dedup_umis"
    PCT_DUPLICATION = "pct_duplication"
    CLIPPED_ALIGNMENTS = "clipped_alignments"
    LOW_MAPQ_ALIGNMENTS = "low_mapq_alignments"
    ALIGNMENTS_PF = "alignments_pf"
    TOTAL_UNIQUE_SITES = "total_unique_sites"
    SAMPLE_GROUP = "sample_group"
    REPLICATE = "replicate"
    CONCENTRATION = "concentration"
    IDENTIFIER = "identifier"
