#!/usr/bin/env python
"""Process all variants from Vlasschaert's script."""


import argparse
import csv
import logging
import re
from typing import Any, Dict, List, Set, Tuple

from pysam import FastaFile  # pylint: disable=no-name-in-module


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


SV_POSSIBLE_ALT = {"<DEL>", "<DUP>", "<INS>", "<INV>", "<CNV>"}


def main():
    """Process all variants from Vlasschaert's script."""
    args = parse_args()

    # Reading the FASTA reference file
    reference = FastaFile(args.reference)

    # Reading and writing the files
    with open(args.in_csv, newline="") as in_csv, \
         open(args.out_tsv, "w") as out_tsv:
        # The CSV reader
        reader = csv.DictReader(in_csv)

        # The new header
        new_header = [
            "post.modif", "SAMPLE", "CHROM", "POS", "REF", "ALT", "TYPE",
            "SYMBOL", "GT", "DP", "VD", "AD", "AF", "VAF", "TopmedClonal",
            "filter_Topmed", "Busque_possible_artifact", "filter_artifacts",
            "aa_pos", "filter_aa_10pct",
        ]
        assert not set(new_header) & set(reader.fieldnames)
        new_header = reader.fieldnames + new_header

        # The header for the TSV
        print(*new_header, sep="\t", file=out_tsv)

        # Parsing the file
        min_vaf = 1
        for row in reader:
            # Generating the new row
            row.update(
                generate_new_row(
                    row=row, vep_fields=args.vep_header.split("|"),
                )
            )

            # Excluding structural variations (if required)
            if args.exclude_sv:
                if row["ALT"] in SV_POSSIBLE_ALT:
                    logger.info(
                        "Removing chr%s:%s:%s:%s for %s: SV",
                        row["CHROM"],
                        row["POS"],
                        row["REF"],
                        row["ALT"],
                        row["SAMPLE"],
                    )
                    continue

            # Checking for the minimal VAF found
            if row["VAF"] < min_vaf:
                min_vaf = row["VAF"]

            # Keeping only VD ≥ threshold and VAF > threshold
            if (row["VAF"] <= args.vaf) or (row["VD"] < args.vd):
                logger.info(
                    "Removing chr%s:%s:%s:%s for %s: VAF=%s, VD=%s",
                    row["CHROM"],
                    row["POS"],
                    row["REF"],
                    row["ALT"],
                    row["SAMPLE"],
                    row["VAF"],
                    row["VD"],
                )
                continue

            # Removing the artifacts
            if row["filter_artifacts"]:
                logger.info(
                    "Removing chr%s:%s:%s:%s for %s: ARTIFACT=%s, VAF=%s",
                    row["CHROM"],
                    row["POS"],
                    row["REF"],
                    row["ALT"],
                    row["SAMPLE"],
                    row["Busque_possible_artifact"],
                    row["VAF"],
                )
                continue

            # Post modifications of whitelist (and other) status
            row = perform_post_modifications(row, set(args.jak2_ok_manual))

            if is_homopolymer(row, reference):
                logger.info(
                    "Removing chr%s:%s:%s:%s for %s: homopolymer frameshift "
                    "(VAF=%s, VD=%s)",
                    row["CHROM"],
                    row["POS"],
                    row["REF"],
                    row["ALT"],
                    row["SAMPLE"],
                    row["VAF"],
                    row["VD"],
                )
                continue

            # Checking position in the peptide
            row.update(annotate_aa_pos(row))

            # Harmonisation of True/False
            for true_false_col in ("filter_artifacts", "filter_Topmed"):
                row[true_false_col] = str(row[true_false_col]).upper()

            # Printing the new file
            print(*[row[field] for field in new_header], sep="\t",
                  file=out_tsv)

        logger.info("Min VAF found is %s", min_vaf)


HOMOPOLYMER_RE = re.compile(r"([ACGT])\1{5}")


def is_homopolymer(row: Dict[str, Any], reference: FastaFile) -> bool:
    """Check for homopolymer in INDELs"""
    if not row["ExonicFunc.refGene"].startswith("frameshift"):
        logger.debug(
            "chr%s:%s:%s:%s for %s: not a frameshift",
            row["CHROM"],
            row["POS"],
            row["REF"],
            row["ALT"],
            row["SAMPLE"],
        )
        return False

    if row["VAF"] >= 0.08 and row["VD"] >= 10:
        logger.debug(
            "chr%s:%s:%s:%s for %s: frameshift, but VAF=%s, VD=%s. DP=%s",
            row["CHROM"],
            row["POS"],
            row["REF"],
            row["ALT"],
            row["SAMPLE"],
            row["VAF"],
            row["VD"],
            row["DP"],
        )
        return False

    # Making sure we have an INDEL
    assert len(row["REF"]) != len(row["ALT"])

    # The position, and starting and ending position
    position = int(row["POS"])
    start = position        # base 0 for FastaFile
    end = position + 1      # end not inclusive for FastaFile (base 0)

    # Checking for deletion
    if len(row["REF"]) > len(row["ALT"]):
        # The start is the position of the deleted nucleotide (ref, base 0) and
        # the end is the position of the deleted nucleotide, plus the length of
        # the deletion.
        start = position
        end = position + len(row["REF"]) - len(row["ALT"])

    # Checking for insertion
    elif len(row["REF"]) < len(row["ALT"]):
        # The start is the position of the insertion (ref, base 0) and the end
        # is the position of the ref, +1 base
        start = position - 1
        end = position + 1

    # The sequence to check
    sequence = reference.fetch(row["CHROM"], start - 4, end + 4)

    # Is this a homopolymer
    if HOMOPOLYMER_RE.search(sequence):
        logger.info(
            "chr%s:%s:%s:%s for %s: homopolymer frameshift (%s)",
            row["CHROM"],
            row["POS"],
            row["REF"],
            row["ALT"],
            row["SAMPLE"],
            sequence,
        )
        return True

    return False


AA_LEN = {
    "NM_015338": 1541,
    "NM_005188": 906,
    "NM_022552": 912,
    "NM_000516": 394,
    "NM_002074": 340,
    "NM_004972": 1132,
    "NM_003620": 605,
    "NM_012433": 1304,
    "NM_003016": 221,
    "NM_001127208": 2002,
    "NM_001126112": 393,
}

AA_POS_RE = re.compile(r"\d+")


def annotate_aa_pos(row: Dict[str, Any]) -> Dict[str, Any]:
    """Check for homopolymer in INDELs"""
    aa_annotation = {
        "aa_pos": "",
        "filter_aa_10pct": "FALSE",
    }

    if row["NonsynOI"] == "" or row["NonsynOI"] == "nan":
        return aa_annotation

    aa_info = row["NonsynOI"].split("_")
    positions = [int(AA_POS_RE.search(_).group()) for _ in aa_info]
    min_pos = min(positions)
    max_pos = max(positions)

    # The 10% of the length of the peptide
    pct_10 = AA_LEN[row["Accession"]] * 0.1

    # Checking at the start of the protein
    if min_pos <= pct_10:
        aa_annotation["aa_pos"] = min_pos
        aa_annotation["filter_aa_10pct"] = "TRUE"

    # Checking at the end of the protein
    elif max_pos >= AA_LEN[row["Accession"]] - pct_10:
        aa_annotation["aa_pos"] = max_pos
        aa_annotation["filter_aa_10pct"] = "TRUE"

    else:
        aa_annotation["aa_pos"] = ",".join(map(str, positions))
        aa_annotation["filter_aa_10pct"] = "FALSE"

    return aa_annotation


# Some variants were missing in the white list, so we need to modify them
POST_MODIFICATIONS = {
    "CREBBP": (
        "NM_004380",
        {"S1680del"},
    ),
    "CSF3R": (
        "NM_000760",
        {"truncatingc.741-791"},
    ),
    "DNMT3A": (
        "NM_022552",
        {"F732del", "F752del"},
    ),
    "EP300": (
        "NM_001429",
        {"VF1148_1149del"},
    ),
    "FLT3": (
        "NM_004119",
        {"FY590-591GD", "del835"},
    ),
    "IDH1": (
        "NM_005896",
        {"V178I"},
    ),
    "KIT": (
        "NM_000222",
        {"del551-559", "del560", "del579", "ins503"},
    ),
    "MPL": (
        "NM_005373",
        {"W515-518KT", "del513"},
    ),
}


def perform_post_modifications(
    row: Dict[str, Any],
    jak2_ok_manual: Set[str],
) -> Dict[str, Any]:
    """Perform post modifications to fix Vlasschaert's script.

    Those annotations were missing, but should have been considered white
    listed.

    +--------+-----------------------------------------------------+
    | Symbol | Variant                                             |
    +========+=====================================================+
    | CREBBP | S1680del                                            |
    +--------+-----------------------------------------------------+
    | CSF3R  | truncatingc.741-791                                 |
    +--------+-----------------------------------------------------+
    | DNMT3A | F732del, F752del                                    |
    +--------+-----------------------------------------------------+
    | EP300  | VF1148_1149del                                      |
    +--------+-----------------------------------------------------+
    | FLT3   | FY590-591GD, del835                                 |
    +--------+-----------------------------------------------------+
    | IDH1   | V178I                                               |
    +--------+-----------------------------------------------------+
    | JAK2   | del/ins537-539L, del/ins538-539L, del/ins540-543MK, |
    |        | del/ins540-544MK, del/ins541-543K, del542-543       |
    |        | del543-544, ins11546-547                            |
    +--------+-----------------------------------------------------+
    | KIT    | del551-559, del560, del579, ins503                  |
    +--------+-----------------------------------------------------+
    | MPL    | W515-518KT, del513                                  |
    +--------+-----------------------------------------------------+

    We also want to exclude variants from ASXL1 which are not located in exon
    11 and 12. These should have `wl.exception` set to FALSE, but `whitelist`
    set to TRUE.

    Finally, a bug was found where synonymous variants were kept, because of
    the way the scripts looks for stop gained / loss (loss of function). It
    looks with a REGEX for a X in the `NonsynOI` column. But if there is a
    synonymous variant which changes a STOP for a STOP (for example, X123X
    means stop à 123 for a stop).

    """
    symbol = row["SYMBOL"]

    # Missing white list
    if symbol in POST_MODIFICATIONS:
        accession, required_variant = POST_MODIFICATIONS[symbol]
        assert row["Accession"] == accession, row["Accession"]
        if row["NonsynOI"] in required_variant and row["whitelist"] != "TRUE":
            row["whitelist"] = "TRUE"
            row["post.modif"] = "TRUE"

    # Confusion about ASXL1 (all variants should be in exons 11 or 12)
    elif symbol == "ASXL1":
        assert row["Accession"] == "NM_015338", row["Accession"]
        if row["whitelist"] == "TRUE":
            if row["wl.exception"] == "FALSE":
                assert (
                    "exon11" not in row["transcriptOI"]
                    and "exon12" not in row["transcriptOI"]
                ), row["transcriptOI"]
                row["whitelist"] = "FALSE"
                row["post.modif"] = "TRUE"

            elif row["wl.exception"] == "TRUE":
                assert (
                    "exon11" in row["transcriptOI"]
                    or "exon12" in row["transcriptOI"]
                ), row["transcriptOI"]

    # After looking at possible manual review INDELS for JAK2, we want to
    # whitelist those
    elif symbol == "JAK2":
        assert row["Accession"] == "NM_004972", row["Accession"]
        if row["manualreview"] == "TRUE":
            if row["NonsynOI"] in jak2_ok_manual:
                row["whitelist"] = "TRUE"
                row["post.modif"] = "TRUE"

    # Synonymous variants
    elif (
        row["ExonicFunc.refGene"].startswith("synonymous")
        and row["whitelist"] == "TRUE"
    ):
        assert row["NonsynOI"].startswith("X")
        assert row["NonsynOI"].endswith("X")
        assert row["wl.lof"] == "TRUE"

        row["whitelist"] = "FALSE"
        row["wl.lof"] = "FALSE"
        row["post.modif"] = "TRUE"

    return row


def generate_new_row(row: Dict[str, str],
                     vep_fields: List[str]) -> Dict[str, Any]:
    """Generate a new row from a current one.

    This function adds the following values:

    - SAMPLE
    - CHROM
    - POS
    - REF
    - ALT
    - SYMBOL
    - GT
    - DP
    - VD
    - AD
    - AF
    - VAF
    - TopmedClonal
    - filter_Topmed
    - Busque_possible_artifact
    - filter_artifacts

    """
    # The INFO field from the VCF is in the 'Otherinfo11' column. Values are
    # split using ';', then there are key=value pairs
    info_field = dict(_.split("=") for _ in row["Otherinfo11"].split(";"))

    # The FORMAT field from the VCF is in the 'Otherinfo12' column. This field
    # gives the value order of the genotype field.
    format_field = {
        name: i for i, name in enumerate(row["Otherinfo12"].split(":"))
    }

    # The genotype field from the VCF is in the 'Otherinfo13' column. This
    # order of this field comes from the FORMAT field.
    geno_field = row["Otherinfo13"].split(":")

    # Some values need to be integer or generated
    vd_value = int(geno_field[format_field["VD"]])
    dp_value = int(geno_field[format_field["DP"]])
    vaf_value = vd_value / dp_value

    # Extra values comming from VEP annotation, if any
    #   - TOPMed
    #   - Artifact
    topmed = ""
    is_topmed = False
    artifact = ""
    is_artifact = False

    # Extracting the VEP annotation
    #   - each entry is split by ','
    #   - each field for each entry is split by '|'
    if "CSQ" in info_field:
        vep_annotations = generate_vep_annotation(
            info_field["CSQ"], vep_fields,
        )

        # The TOPMed annotation
        topmed, is_topmed = check_is_topmed(
            annotations=[_["TopmedClonal"] for _ in vep_annotations],
        )

        # Is this a Busque artifact
        artifact, is_artifact = check_is_artifact(
            annotations=[
                _["Busque_possible_artifact"] for _ in vep_annotations
            ],
            vaf=vaf_value,
        )

    # Some value are straight from the CSV file
    return {
        "SAMPLE": info_field["SAMPLE"],
        "CHROM": row["Otherinfo4"],
        "POS": row["Otherinfo5"],
        "REF": row["Otherinfo7"],
        "ALT": row["Otherinfo8"],
        "TYPE": info_field["TYPE"],
        "SYMBOL": row["Gene.refGene"],
        "GT": geno_field[format_field["GT"]],
        "DP": dp_value,
        "VD": vd_value,
        "AD": geno_field[format_field["AD"]],
        "AF": geno_field[format_field["AF"]],
        "VAF": vaf_value,
        "TopmedClonal": topmed,
        "filter_Topmed": is_topmed,
        "Busque_possible_artifact": artifact,
        "filter_artifacts": is_artifact,
        "post.modif": "FALSE",
    }


def check_is_topmed(annotations: List[str]) -> Tuple[str, bool]:
    """Check if this is a TOPMed CHIP."""
    # Sanity check
    for _ in annotations[1:]:
        assert _ == annotations[0]

    return annotations[0], annotations[0] != ""


def check_is_artifact(annotations: List[str], vaf: float) -> Tuple[str, bool]:
    """Check if this is a Busque artifact.

    Note that if value is 'NM_015338.6:c.1934dupG', we need to check if the VAF
    is less than 6% before excluding.

    When value is True, it means it's an artifact.

    """
    # Sanity check
    for _ in annotations[1:]:
        assert _ == annotations[0]

    # Is this an artifact (i.e. a value in the column)
    is_artifact = annotations[0] != ""

    # We have a special case for ASXL1, we need to exclude as artifact if
    #   - annotation is 'NM_015338.6:c.1934dupG'
    #   - VAF > 0.06
    is_good_asxl1 = (
        (annotations[0] == "NM_015338.6:c.1934dupG")
        and (vaf >= 0.06)
    )

    # Here is the truth table (to explain the XOR)
    #     is_artifact: F F T T
    #   is_good_asxl1: F T F T
    #         exclude: F T T F
    # Note that if is_artifact is F, is_not_artifact can never be T (since the
    # "Busque_possible_artifact" column cannot be empty (i.e. "") and
    # "NM_015338.6:c.1934dupG" at the same time)
    return annotations[0], is_artifact ^ is_good_asxl1


def generate_vep_annotation(csq_field: str,
                            vep_fields: List[str]) -> List[Dict[str, str]]:
    """Generate VEP annotation."""
    annotations = [_.split("|") for _ in csq_field.split(",")]

    # Sanity check
    for _ in annotations:
        assert len(_) == len(vep_fields)

    annotations = [
        dict(zip(vep_fields, _)) for _ in annotations
    ]

    return annotations


def parse_args() -> argparse.Namespace:
    """Parses the arguments and function."""
    parser = argparse.ArgumentParser(
        description="Process all variants from Vlasschaert's script.",
    )

    parser.add_argument(
        "-r", "--reference", type=str, metavar="FASTA", required=True,
        help="The reference file (FASTA).",
    )

    parser.add_argument(
        "--vd", type=int, metavar="INT", default=2,
        help="The VD threshold to keep. [ ≥%(default)s ]",
    )

    parser.add_argument(
        "--vaf", type=float, metavar="FLOAT", default=0.001,
        help="The VAF threshold. [ >%(default)s ]",
    )

    parser.add_argument(
        "--vep-header", type=str, metavar="STR", required=True,
        help="The VEP header (from original VCF files)."
    )

    parser.add_argument(
        "--jak2-ok-manual", type=str, nargs="+", default=[],
        help="Insertions / deletions flagged as manual review, and should be "
             "set to whitelist.",
    )

    parser.add_argument(
        "--exclude-sv", action="store_true",
        help="Exclude structural variation from the dataset. For SV generated "
             "by VarDict, the alternative allele is coded as either <DEL>, "
             "<DUP>, <INV>, <INS> and <CNV>.",
    )

    parser.add_argument(
        dest="in_csv", type=str, metavar="CSV",
        help="The input file (Vlasschaert's results, CSV).",
    )

    parser.add_argument(
        dest="out_tsv", type=str, metavar="TSV",
        help="The output of this script (TSV).",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
