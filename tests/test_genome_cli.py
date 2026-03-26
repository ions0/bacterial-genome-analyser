"""
Bacterial Genome Analyser: tests/test_genome_cli.py

"""
import pytest
import argparse

from genome_cli import validate_arguments

def test_validate_arguments(tmp_path):

    test_genome_file = tmp_path / "test_genome.gb"
    test_genome_file.touch()

    test_parser = argparse.Namespace(genome=test_genome_file, output=None, no_display=False)

    arg_dict = validate_arguments(test_parser)

    assert arg_dict["genome_path"].exists()
    assert arg_dict["file_format"] == "genbank"
    assert arg_dict["output_dir"] == None
    assert arg_dict["display_plots"] == True

def test_invalid_extension(tmp_path):

    test_genome_file = tmp_path / "test_genome.gx"
    test_genome_file.touch()

    test_parser = argparse.Namespace(genome=test_genome_file, output=None, no_display=False)

    with pytest.raises(ValueError):
        arg_dict = validate_arguments(test_parser)
    
def test_filenotfound(tmp_path):

    test_genome_file = tmp_path / "test_genome.gx"

    test_parser = argparse.Namespace(genome=test_genome_file, output=None, no_display=False)

    with pytest.raises(ValueError):
        arg_dict = validate_arguments(test_parser)