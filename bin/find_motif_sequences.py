#!/usr/bin/env python3


#########################################################################
#########################################################################
## Helper methods for determining the best-match sequence to the       ##
## on-target motif within close proximity of a reported site           ##
#########################################################################
#########################################################################


from typing import Dict, Tuple, Callable
import itertools
import parasail
import pysam

from site_seq_utils.coordinate import Zerodinate


def get_motif_search_fn(
    search_motif: str,
    fasta_path: str,
    search_factor: float,
    match_score: int,
    mismatch_pen: int,
    gap_pen: int,
) -> Callable:
    """
    Return a function that will find the best-match sequence for a given genomic
    position using the given search motif (i.e. on-target motif), genome FASTA,
    search factor (how far around the genomic position to search), and alignment
    scoring parameters
    """
    fasta = pysam.FastaFile(fasta_path)
    scoring_matrix, base_index_map, base_match_map = get_parasail_scoring_matrix(
        match_score, -mismatch_pen
    )
    rev_comp_motif = reverse_complement(search_motif)
    search_window = search_factor * len(search_motif) // 1

    def find_best_motif(chrom: str, position: int) -> Tuple[str, str, int, int, int]:
        window_start = max(position - search_window, 0)
        window_end = position + search_window
        sequence = fasta.fetch(chrom, window_start, window_end)

        aln_fwd = parasail.sg_dx_trace(
            search_motif, sequence, gap_pen, gap_pen, scoring_matrix
        )
        aln_rev = parasail.sg_dx_trace(
            rev_comp_motif, sequence, gap_pen, gap_pen, scoring_matrix
        )
        if aln_fwd.score >= aln_rev.score:
            best_aln = aln_fwd
            best_query = search_motif
            aln_is_reverse = False
        else:
            best_aln = aln_rev
            best_query = rev_comp_motif
            aln_is_reverse = True

        aln_score = best_aln.score
        aln_ref_start, aln_ref_end, formatted_seq, num_subs, num_gaps = (
            format_and_count_alignment(
                best_query,
                sequence,
                best_aln.cigar,
                base_match_map,
            )
        )
        ref_start = window_start + aln_ref_start
        ref_end = window_start + aln_ref_end
        if aln_is_reverse:
            formatted_seq = reverse_complement(formatted_seq)
        motif_location = (
            Zerodinate(chrom, ref_start, ref_end, aln_is_reverse)
            .to_onedinate()
            .string_value
        )
        return formatted_seq, motif_location, aln_score, num_subs, num_gaps

    return find_best_motif


def get_parasail_scoring_matrix(
    match_score: int,
    mismatch_pen: int,
) -> parasail.Matrix:
    """
    Return a parasail scoring matrix that allows for degenerate base matching
    """
    match_map = {
        "A": "A",
        "G": "G",
        "C": "C",
        "T": "T",
        "N": "ACGT",
        "R": "AG",
        "Y": "CT",
        "S": "GC",
        "W": "AT",
        "K": "GT",
        "M": "AC",
        "V": "ACG",
        "H": "ACT",
        "D": "AGT",
        "B": "CGT",
    }
    uppercase_alphabet = "".join(match_map)
    alphabet = uppercase_alphabet + uppercase_alphabet.lower()
    scoring_matrix = parasail.matrix_create(alphabet, match_score, mismatch_pen)
    base_idx_map = {base: idx for idx, base in enumerate(alphabet)}
    for base1, base2 in itertools.combinations(match_map, 2):
        uidx1, uidx2 = base_idx_map[base1], base_idx_map[base2]
        lidx1, lidx2 = base_idx_map[base1.lower()], base_idx_map[base2.lower()]
        if set(match_map[base1]).intersection(match_map[base2]):
            val = match_score
        else:
            val = mismatch_pen

        for uidx, lidx in ((uidx1, lidx1), (uidx2, lidx2)):
            scoring_matrix[uidx, lidx] = match_score
            scoring_matrix[uidx, lidx] = match_score

        for idx1, idx2 in itertools.product((uidx1, lidx1), (uidx2, lidx2)):
            scoring_matrix[idx1, idx2] = val
            scoring_matrix[idx2, idx1] = val
    return scoring_matrix, base_idx_map, match_map


def format_and_count_alignment(
    query: str, ref: str, cigar: parasail.Cigar, base_match_map: Dict[str, str]
) -> Tuple[int, int, str, int, int]:
    """
    Return the reference sequence formatted with mismatches as lowercase and gaps as dashes.
    Also count mismatches and indels.
    """

    def degen_match(qbase: str, rbase: str) -> bool:
        # Assumes qbase != rbase
        # Does not count degen ref bases
        for equiv_qbase in base_match_map[qbase]:
            if equiv_qbase == rbase:
                return True
        return False

    formatted_sequence = []
    substitutions = 0
    indels = 0
    query_pos = prev_query_pos = cigar.beg_query
    ref_pos = prev_ref_pos = cigar.beg_ref
    query_length = len(query)
    ref_start = ref_end = None
    for cigar_code in cigar.seq:
        op = cigar.decode_op(cigar_code)
        length = cigar.decode_len(cigar_code)
        if op == b"=":  # Match
            ref_pos += length
            formatted_sequence.append(ref[prev_ref_pos:ref_pos].upper())
            query_pos += length
        elif op == b"X":  # Mismatch (substitution)
            query_pos += length
            ref_pos += length
            query_sub_seq = query[prev_query_pos:query_pos].upper()
            ref_sub_seq = ref[prev_ref_pos:ref_pos].upper()
            # Check whether "mismatch" is to a matching degenerate base
            for query_base, ref_base in zip(query_sub_seq, ref_sub_seq):
                if degen_match(query_base, ref_base):
                    formatted_sequence.append(ref_base)
                else:
                    formatted_sequence.append(ref_base.lower())
                    substitutions += 1
        elif op == b"D":
            ref_pos += length
            if 0 < query_pos < query_length:
                indels += length
                formatted_sequence.append(ref[prev_ref_pos:ref_pos].lower())
        elif op == b"I":  # Insertion in the query sequence (gap in reference)
            query_pos += length
            formatted_sequence.append("-")
            indels += length
        else:
            raise Exception(f'Unrecognized CIGAR operation "{op}"')

        if ref_start is None and query_pos > 0:
            ref_start = prev_ref_pos

        if ref_end is None and query_pos == query_length:
            ref_end = ref_pos

        prev_query_pos = query_pos
        prev_ref_pos = ref_pos

    return ref_start, ref_end, "".join(formatted_sequence), substitutions, indels


DNA_BASE_COMPLEMENT = {
    "A": "T",
    "C": "G",
    "T": "A",
    "G": "C",
    "N": "N",
    "a": "t",
    "c": "g",
    "t": "a",
    "g": "c",
    "n": "n",
    "r": "y",
    "y": "r",
    "R": "Y",
    "Y": "R",
    "s": "s",
    "w": "w",
    "S": "S",
    "W": "W",
    "k": "m",
    "m": "k",
    "K": "M",
    "M": "K",
    "v": "b",
    "b": "v",
    "V": "B",
    "B": "V",
    "d": "h",
    "h": "d",
    "D": "H",
    "H": "D",
    "-": "-",
}


def reverse_complement(string: str) -> str:
    return "".join(DNA_BASE_COMPLEMENT[b] for b in string[::-1])
