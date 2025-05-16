#########################################################################
#########################################################################
##         Helper methods for formatting values for display            ##
#########################################################################
#########################################################################

from typing import List, Union


def format_dna_sequence(motif: str) -> str:
    html_str = "<span style='font-family:monospace'>"
    for letter in motif:
        if letter.upper() == "A":
            html_str += '<span style="color:#008000">{}</span>'.format(letter)
        elif letter.upper() == "C":
            html_str += '<span style="color:#0000FF">{}</span>'.format(letter)
        elif letter.upper() == "G":
            html_str += '<span style="color:#FFC300">{}</span>'.format(letter)
        elif letter.upper() == "T":
            html_str += '<span style="color:#FF0000">{}</span>'.format(letter)
        else:
            html_str += '<span style="color:#000000">{}</span>'.format(letter)
    html_str += "</span>"
    return html_str


def format_coordinate(coord: str) -> str:
    chrom, loc_str = coord.split(":")
    strand = loc_str.split(",")[1]
    locs = loc_str.split(",")[0].split("-")
    formatted_loc_str = "-".join("{:,}".format(int(loc)) for loc in locs)
    return f"{chrom}:{formatted_loc_str},{strand}"


def format_number_list(number_list: List[Union[int, float]]) -> str:
    try:
        numbers = eval(number_list)
    except Exception:
        return number_list
    if isinstance(number_list[0], int):
        return "; ".join("{:,}".format(n) for n in numbers)
    else:
        return "; ".join("{:,.2f}".format(n) for n in numbers)


def format_percentage_list(percentage_list: List[float]) -> str:
    try:
        percentages = eval(percentage_list)
    except Exception:
        return percentage_list
    return "; ".join(format_percentage(p) for p in percentages)


def format_integer(number: int) -> str:
    return "{:,}".format(number)


def format_decimal(number: float) -> str:
    return "{:,.2f}".format(number) if number is not None else None


def format_pvalue(pvalue: float) -> str:
    return "{:.3e}".format(pvalue) if pvalue is not None else None


def format_percentage(percentage: float) -> str:
    return "{:.2%}".format(percentage / 100) if percentage is not None else None
