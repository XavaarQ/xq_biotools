"""
Codon Harmonizer
David Bauwens, Xavaar
"""
import itertools
import sys
import tempfile

import numpy as np
import pandas as pd
import seaborn as sns

sns.set()
import logging
import re
from pathlib import Path

import Bio
import ete3
import matplotlib.pyplot as plt
import requests
from BCBio import GFF
from beartype import beartype
from beartype.typing import List, Tuple, Union
from bs4 import BeautifulSoup
from gff2genbank import _split_fasta_from_gff, gff2genbank


def codon_harmonizer_from_file(
    sequences: List[str],
    source_codon_usage_table: pd.core.frame.DataFrame,
    target_codon_usage_table: pd.core.frame.DataFrame,
) -> Path:
    """
    Transcribes sequences from .gb files into version harmonized for a specific organism
    Codon tables are derived from https://www.kazusa.or.jp/codon

    Parameters
    ----------
    sequences : (file)
        File which containers the sequences to be harmonized,

    source_codon_usage_table : (pd.core.frame.DataFrame)

    target_codon_usage_table : (pd.core.frame.DataFrame)

    Return
    ----------
    A file with the harmonized sequences

    """

    main_wd = Path().absolute()

    name_prefix = "codon_harmonizer_files"
    named_dir = create_named_temporary_directory(name_prefix)

    recs = turn_file_into_seq_records(named_dir, sequences)

    new_recs = []
    for rec in recs:
        strSeq = str(rec.seq)

        newSeq = codon_harmonizer(strSeq, source_codon_usage_table, target_codon_usage_table)

        rec.seq = Bio.Seq.Seq(newSeq)
        rec.name = "Harmonized_" + rec.name
        rec.name = rec.name.replace(" ", "").strip()
        new_recs.append(rec)

    # Make subdirectories in the output dir
    log_dir = Path(named_dir, "logs")
    harmonized_seqs_dir = Path(named_dir, "harmonized_sequence")
    for _dir in [log_dir, harmonized_seqs_dir]:
        _dir.mkdir(parents=True)

    outfile_name = harmonized_seqs_dir / f"{sequences.name.split('.')[-2]}_harmonized_seqs"

    Path("codon_usage_comparison.png").rename(log_dir / "codon_usage_comparison.png")
    Path("codon_comparison.csv").rename((log_dir / "codon_comparison.csv"))
    write_gb(outfile_name, new_recs, "gb")
    print(f"New sequences written to {outfile_name}")

    return outfile_name


def create_named_directory(name_prefix):
    path = Path(name_prefix)
    path.mkdir()
    return path


def turn_file_into_seq_records(temp_dir, arg_file):
    if isinstance(arg_file, str):
        arg_file = open(arg_file, "r", encoding="utf-8")

    arg_file_path = Path(arg_file.name)
    file_ext = arg_file_path.suffix

    if file_ext in (".gff", ".gff3"):
        new_gff, fasta = _split_fasta_from_gff(arg_file, temp_dir)
        new_gb = gff2genbank(new_gff, fasta_file=fasta, molecule_type="DNA")
        recs = get_gb_seq(new_gb)

    elif file_ext in (".gb", ".gbk"):
        recs = get_gb_seq(arg_file)
    elif file_ext in (".fa", ".fasta"):
        recs = get_fasta_seq(arg_file)
    else:
        raise Exception(f"Unacceptable file type {file_ext}")
        logging.error(f"Unacceptable file type {file_ext}")
    return recs


def get_gb_seq(in_file):
    #     seq = Bio.SeqIO.parse(in_file, "genbank")
    recs = []
    counter = 0
    try:  # Some molecules gbs fail for an unexplained reason
        for rec in Bio.SeqIO.parse(in_file, "genbank"):
            #             print(rec)
            #             print('\n\n\n')

            recs.append(rec)
            counter += 1

    except AttributeError as error:
        logging.warning(
            f"Entry number {counter} in the file has a problem with it. It will NOT be turned into a molecule\r\n"
        )

    #     recs = [rec for rec in Bio.SeqIO.parse(in_file, "genbank")]
    #     sys.exit()
    return recs


def get_fasta_seq(fasta_file):
    sequence_record = list(Bio.SeqIO.parse(fasta_file, "fasta"))

    # assign molecule type
    for seq in sequence_record:
        seq.annotations["molecule_type"] = "DNA"

    return sequence_record


# codon_harmonizer
def codon_harmonizer(
    sequence: str,
    source_codon_usage_table: pd.core.frame.DataFrame = pd.DataFrame(),
    target_codon_usage_table: pd.core.frame.DataFrame = pd.DataFrame(),
    source_species: str = False,
    target_species: str = False,
) -> str:
    """
    Transcribes sequences  into version harmonized for a specific organism
    Codon table are derived from https://www.kazusa.or.jp/codon

    Parameters
    ----------
    sequences : (str)
        Nucleotide sequences to be harmonized,

    source_species: (str)
        Pass the scientific name of a species (eg. Homo sapiens) OR the NCBI taxnomic ID of a species (eg. 9606)
        If this value isn't passed then source_codon_usage_table must be provided

    target_species: (str)
        Pass the scientific name of a species (eg. Homo sapiens) OR the NCBI taxnomic ID of a species (eg. 9606)
        If this value isn't passed then target_codon_usage_table must be provided

    source_codon_usage_table: pd.core.frame.DataFrame = pd.DataFrame(),

    target_codon_usage_table: pd.core.frame.DataFrame = pd.DataFrame(),

    Return
    ----------
    harmonized_sequence : (str)

    """

    if target_species:
        target_codon_usage_table = get_codon_table_for_species(target_species)
    if source_species:
        source_codon_usage_table = get_codon_table_for_species(source_species)

    if source_codon_usage_table.empty or target_codon_usage_table.empty:
        raise Exception(
            f"Missing one of the required codon tables, enter a source and target species name/NCBI id or pass a custom codon table directly"
        )

    sequence = sequence.upper().replace("U", "T")
    source_codon_usage_table["codon"] = source_codon_usage_table["codon"].str.replace("U", "T")
    target_codon_usage_table["codon"] = target_codon_usage_table["codon"].str.replace("U", "T")

    ## Break up your string in codons
    n = 3
    sequence = [sequence[i : i + n] for i in range(0, len(sequence), n)]
    df = pd.DataFrame(sequence, columns=["codon_source"])
    # df.to_csv("df.csv")

    ##Use your codon usage table to give the percentage use and the amino-acid it codes for from your string
    df_final_in = df.merge(
        source_codon_usage_table, how="left", left_on="codon_source", right_on="codon"
    ).drop(columns="codon")

    df_final_in["fraction"].plot()

    ##Without harmonizing, how would your sequnece code in the target organism

    df_final_out = df.merge(
        source_codon_usage_table, how="left", left_on="codon_source", right_on="codon"
    ).drop(columns="codon")
    df_final_out["fraction"].plot()

    ## Now what is the absolute difference between host and target organism

    df_diff = abs(df_final_out["fraction"] - df_final_in["fraction"])
    df_diff.plot()

    ## Here we merge the 2 codon usage tables on amino acid (from the source and from the target) and calculate the difference in codon usage percentage per combination
    df2 = pd.merge(
        source_codon_usage_table,
        target_codon_usage_table,
        how="outer",
        on="amino_acid",
        suffixes=["_source", "_target"],
    )
    df2["diff"] = abs(df2["fraction_source"] - df2["fraction_target"])
    df2

    ## Now we're going to find the minimal difference per source codon
    df3 = df2.groupby(["codon_source"])["diff"].min()
    df3 = df3.reset_index()

    ## If we now merge the minimal selection of the source codon and select on the minimal difference, we'll get the amino-acid they are coding for and the target codons.
    df5 = pd.merge(df3, df2, on=["codon_source", "diff"], how="left").drop_duplicates(
        subset=["codon_source"]
    )

    df = df.reset_index()

    df6 = pd.merge(df, df5, on=["codon_source"], how="left")
    df6 = (
        df6.reset_index().dropna()
    )  # If there are nucleotide that don't fit in a codon (and extra T for example) this will eliminate them from the sequence

    a = df6["codon_target"].to_list()

    harmonized_sequence = "".join(a)

    harmonized_sequence = harmonized_sequence.replace("U", "T")

    if not df_final_in["amino_acid"].to_list() == df6["amino_acid"].to_list():
        raise Error("The Amnino Acid sequences are not identical")

    df_diff2 = abs(df6["fraction_target"] - df_final_in["fraction"])

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle("Codon_usage_comparison")

    sns.lineplot(ax=axes[0, 0], data=df6, x="index", y="fraction_source")
    sns.lineplot(ax=axes[0, 1], data=df6, x="index", y="fraction_target")
    sns.lineplot(ax=axes[0, 2], data=df_final_out.reset_index(), x="index", y="fraction")
    sns.lineplot(ax=axes[1, 0], data=df_diff.reset_index(), x="index", y="fraction").set_ylim(0, 1)
    sns.lineplot(ax=axes[1, 1], data=df_diff2.reset_index(), x="index", y=0).set_ylim(0, 1)

    plt.savefig("codon_usage_comparison.png")
    df6.to_csv("codon_comparison.csv")

    return harmonized_sequence


def check_taxonomic_nomenclature(in_string):
    numeric_pattern = re.compile("[0-9]")
    numeric_match = numeric_pattern.match(in_string)

    taxnomic_nomenclature_pattern = re.compile("[A-Z]{1}[a-z]+\s[a-z]+")
    taxonomy_match = taxnomic_nomenclature_pattern.match(in_string)

    if taxonomy_match is None and numeric_match is None:
        raise Exception(
            f"{in_string} does not follow standard taxonmic nomenclature. Genus and species must be present. Genus must start with a capital and all other characters must be lowercase. \
        Example: 'Homo sapiens'"
        )

    if not numeric_match is None:
        return False
    else:
        return True


def write_gb(out_filename: str, output_records: list, filetype: str):
    with open(f"{out_filename}.{filetype}", "w+") as handle:
        for record in output_records:
            try:
                if filetype == "gff":
                    GFF.write([record], handle, include_fasta=True)
                else:
                    Bio.SeqIO.write(record, handle, filetype)

            except NameError as error:
                logging.warning(
                    f"{error}: aka. Molecules with whitespace in the name can't be written to genbank:\n\r {record.name} \n\r"
                )

    return out_filename


# for string in ['9606','Homo sapiens', 'Homo Sapiens', 'homo sapiens','HoMo sapiens','sapiens','Homosapiens']:
#               check_taxonomic_nomenclature(string)


def codon_database_scraper(species_ncbi_id: str):
    url = f"http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species={species_ncbi_id}&aa=1&style=N"

    # Create object page
    page = requests.get(url)

    # Obtain page's information
    soup = BeautifulSoup(page.text, "lxml")  # parser-lxml = Change html to Python friendly format

    try:
        table1 = soup.find("pre").find(text=True)
    except:
        raise Exception(
            f"{url} for {ete3.NCBITaxa().get_taxid_translator([int(species_ncbi_id)])} is empty. If you are using a common organism such as E. coli\
 try passing a numeric NCBI id for a specific substrain instead"
        )

    table2 = (
        table1.replace("(", "").replace(")", "").replace("\n", " ").replace("  ", " ").split(" ")
    )
    table2 = [i for i in table2 if i]

    headers = [
        "codon",
        "amino_acid",
        "fraction",
        "frequency",
        "per_thousand",
    ]
    table_dict = dict([(x, []) for x in headers])

    i = 0
    for x in table2:
        table_dict[list(table_dict.keys())[i % 5]].append(x)
        i += 1

    frequency_table = pd.DataFrame(table_dict)
    frequency_table[
        [
            "fraction",
            "frequency",
            "per_thousand",
        ]
    ] = frequency_table[
        [
            "fraction",
            "frequency",
            "per_thousand",
        ]
    ].astype(float)

    return frequency_table


def get_codon_table_for_species(species: Union[str, int]):
    ncbi = ete3.NCBITaxa()
    species = str(species)
    if check_taxonomic_nomenclature(species):
        species_id = str(ncbi.get_name_translator([species])[species][0])
    else:
        species_id = species

    return codon_database_scraper(species_id)
