"""
Bacterial Genome Analyser: tests/test_genome_io.py

"""
import pytest
from pathlib import Path

from Bio.SeqFeature import SeqFeature, SimpleLocation

from genome_io import (
    load_sequence,
    extract_gene_info
)

#========================= load_sequence() tests ========================================#

def test_load_sequence_format():

    test_format = "fasta"
    test_path = Path()

    with pytest.raises(ValueError):
        load_sequence(test_path, test_format)


def test_load_sequence_path():

    test_format = "genbank"
    test_path = Path("/not/real/")

    with pytest.raises(FileNotFoundError):
        load_sequence(test_path, test_format)

# TODO: test that load_sequence raises RuntimeError on unparseable file
def test_load_sequence_invalid_file():
    pass

# TODO: test that load_sequence raises ValueError on empty sequence
def test_load_sequence_empty_sequence():
    pass

#========================= extract_gene_info() tests ====================================#

def test_extract_gene_info_returns_expected_keys():
     
    seq_feat = SeqFeature(SimpleLocation(5, 10, strand=+1), type="CDS")
    seq_feat.qualifiers["gene"] = ["test"]
    seq_feat.qualifiers["product"] = ["test"]

    test_dict = extract_gene_info(seq_feat)
    print(test_dict["strand"])

    assert test_dict["gene"] == "test"
    assert test_dict["product"] == "test"
    assert test_dict["start"] == 5
    assert test_dict["end"] == 10
    assert test_dict["length"] == 5
    assert test_dict["strand"] == "+"

def test_extract_gene_info_unknown_gene():
    
    seq_feat = SeqFeature(SimpleLocation(5, 10))
    test_dict = extract_gene_info(seq_feat)

    assert test_dict["gene"] == "Unknown"

def test_extract_gene_info_missing_location():
    
    seq_feat = SeqFeature(location=None)

    with pytest.raises(ValueError):
        test_dict = extract_gene_info(seq_feat)
