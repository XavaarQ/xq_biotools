##Original
##https://github.com/chapmanb/bcbb/blob/master/gff/Scripts/gff/gff_to_genbank.py
import io
import logging
import warnings
from pathlib import Path

import Bio
from BCBio import GFF
from beartype import beartype
from beartype.typing import Iterable, Tuple, Union


@beartype
def gff2genbank(
    gff_file: Union[str, io.IOBase, Path],
    fasta_file: Union[str, None] = None,
    molecule_type: str = "DNA",
) -> Path:
    gff_file_name = Path(gff_file if isinstance(gff_file, str) else gff_file.name)
    out_file = gff_file_name.with_suffix(".gb")
    fasta_input = Bio.SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta")) if fasta_file else {}
    gff_iter = GFF.parse(gff_file, fasta_input)
    Bio.SeqIO.write(_check_gff(_fix_ncbi_id(gff_iter), molecule_type), out_file, "genbank")
    return out_file


@beartype
def _fix_ncbi_id(fasta_iter: Iterable) -> Iterable[Bio.SeqRecord.SeqRecord]:
    """
    GenBank identifiers can only be 16 characters; try to shorten NCBI.
    """
    for rec in fasta_iter:
        if len(rec.name) > 16 and "|" in rec.name:
            new_id = rec.name.split("|")[-1]
            warnings.warn(f"Warning: shortening NCBI name {rec.id} to {new_id}", stacklevel=3)
            rec.id = new_id
            rec.name = new_id

        yield rec


@beartype
def _check_gff(gff_iterator: Iterable, molecule_type: str) -> Iterable[Bio.SeqRecord.SeqRecord]:
    """
    Check GFF files before feeding to Bio.SeqIO to be sure they have sequences.
    Adds molecule type if not present.
    """
    for record in gff_iterator:
        if "molecule_type" not in record.annotations:
            record.annotations["molecule_type"] = molecule_type
        yield _flatten_features(record)


@beartype
def _flatten_features(record: Bio.SeqRecord.SeqRecord) -> Bio.SeqRecord.SeqRecord:
    """
    Make sub_features in an input record flat for output.
    GenBank does not handle nested features, so we make everything top level.
    """
    out = []
    stack = record.features[:]
    while stack:
        current_feature = stack.pop()
        out.append(current_feature)
        stack.extend(current_feature.sub_features)
    record.features = out
    return record


@beartype
def _split_fasta_from_gff(
    in_file: Union[str, Path], temp_dir: Union[str, Path]
) -> Tuple[Path, Path]:
    base_name = Path(in_file).stem
    new_gff = Path(temp_dir, f"{base_name}_2.gff")
    fasta_file = Path(temp_dir, f"{base_name}.fasta")

    out_files = [new_gff, fasta_file]
    out_file_idx = 0

    with Path(in_file).open() as f:
        for line in f:
            line = line.rstrip()
            if "##FASTA" in line:
                out_file_idx += 1
                if out_file_idx >= len(out_files):
                    logging.error("Unexpected '##FASTA' line encountered")
                    raise ValueError("Unexpected '##FASTA' line encountered")
                continue

            with Path(out_files[out_file_idx]).open("a") as out_file:
                out_file.write(f"{line}\\n")

    if out_file_idx != len(out_files) - 1:
        logging.error("Did not encounter expected number of '##FASTA' lines")
        raise ValueError("Did not encounter expected number of '##FASTA' lines")

    return new_gff, fasta_file
