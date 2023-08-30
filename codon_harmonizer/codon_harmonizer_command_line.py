import argparse

import pandas as pd

from codon_harmonizer import codon_harmonizer_from_file, get_codon_table_for_species


def argInputs():
    parser = argparse.ArgumentParser(
        prog="Codon Harmonizer",
        description="Harmonizes Codons",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    g1 = parser.add_argument_group("Required Inputs")
    g2 = parser.add_argument_group("Optional Inputs")
    g3 = parser.add_argument_group("Advanced and Developer Options")

    g1.add_argument(
        "--sequences",
        "-seqs",
        metavar="Files for the sequences to be harmonized",
        help="Upload a file (.gb/.fasta/.gff) of the sequenced you want harmonized",
        type=argparse.FileType("r"),
        required=True,
    )

    g1.add_argument("--source_species", "-srcS", type=str)
    g1.add_argument("--target_species", "-tgtS", type=str)
    g2.add_argument("--source_codon_usage_table", "-srcT", type=argparse.FileType("r"))
    g2.add_argument("--target_codon_usage_table", "-tgtT", type=argparse.FileType("r"))

    args = parser.parse_args()

    # Check for errors
    if args.source_species is None and args.source_codon_usage_table is None:
        parser.error("Either --source_species or --source_codon_usage_table must be provided.")
    elif args.source_species and args.source_codon_usage_table:
        parser.error(
            "Cannot provide both --source_species and --source_codon_usage_table. Choose one."
        )

    if args.target_species is None and args.target_codon_usage_table is None:
        parser.error("Either --target_species or --target_codon_usage_table must be provided.")
    elif args.target_species and args.target_codon_usage_table:
        parser.error(
            "Cannot provide both --target_species and --target_codon_usage_table. Choose one."
        )

    return args


def main():
    args = argInputs()

    if args.source_species:
        src_DF = get_codon_table_for_species(args.source_species)
    else:
        src_DF = pd.read_csv(args.source_codon_usage_table)
        if src_DF.empty:
            raise ValueError(
                "The source codon usage table CSV file is empty or not formatted correctly."
            )

    if args.target_species:
        tgt_DF = get_codon_table_for_species(args.target_species)
    else:
        tgt_DF = pd.read_csv(args.target_codon_usage_table)
        if tgt_DF.empty:
            raise ValueError(
                "The target codon usage table CSV file is empty or not formatted correctly."
            )

    codon_harmonizer_from_file(args.sequences, src_DF, tgt_DF)


if __name__ == "__main__":
    main()
