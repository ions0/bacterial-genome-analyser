"""
Bacterial Genome Analyser: tests/test_genome_analysis.py

"""
import pytest
import pandas as pd

from genome_analysis import (
    get_gc, 
    calculate_strand_switches, 
    calculate_gene_stats, 
    calculate_gc_window,
    calculate_gene_density,
    split_genes_by_strand
)

#========================= get_gc() tests ===============================================#

def test_get_gc_mixed_bases():

    s = "ATGC"
    assert get_gc(s) == 2

def test_get_gc_all_gc():

    s = "GCGCGCGC"
    assert get_gc(s) == 8

def test_get_gc_empty_string():

    s = ""
    assert get_gc(s) == 0

def test_get_gc_lowercase():

    s = "atgc"
    assert get_gc(s) == 2

#========================= calculate_strand_switches() tests ============================#

def test_strand_switch_one_switch():

    test_gene_data = [
        {"start": 100, "strand": "+"},
        {"start": 200, "strand": "-"},
    ]

    assert calculate_strand_switches(test_gene_data) == 1

def test_strand_switch_no_switch():

    test_gene_data = [
        {"start": 100, "strand": "+"},
        {"start": 200, "strand": "+"},
    ]

    assert calculate_strand_switches(test_gene_data) == 0

def test_strand_switch_empty_list():
    
    test_gene_data = []

    with pytest.raises(ValueError):
        calculate_strand_switches(test_gene_data)

def test_strand_switch_order():

    test_gene_data = [
        {"start": 300, "strand": "+"},
        {"start": 200, "strand": "-"},
        {"start": 100, "strand": "+"},
    ]

    assert calculate_strand_switches(test_gene_data) == 2

#========================= calculate_gene_stats() tests =================================#

def test_gene_stats():
    
    test_gene_data = [
        {"length": 100},
        {"length": 200},
        {"length": 300}
    ]

    s = "AAAAAAAAAA"
    test_dict = calculate_gene_stats(test_gene_data, s)

    assert test_dict["average_length"] == 200
    assert test_dict["coding_percentage"] == 6000

def test_gene_stats_empty():

    test_gene_data = []
    s = "AAAAAAAAAA"

    with pytest.raises(ValueError):
        calculate_gene_stats(test_gene_data, s)

#========================= calculate_gc_window() tests ==================================#

def test_gc_window_size():

    seq = "TTAGCTAGC"
    ws = 10

    with pytest.raises(ValueError):
        calculate_gc_window(seq, ws)

def test_gc_window_no_gc():

    seq = "TATATATATA"
    _, gc_values = calculate_gc_window(seq, 5)

    assert gc_values == [0.0]

#========================= calculate_gene_density() tests ===============================#

def test_calc_gene_density_empty():

    test_gene_data = pd.DataFrame()

    with pytest.raises(ValueError):
        pos, density = calculate_gene_density(test_gene_data, 1000, 10, 5)
    
def test_calc_gene_density_expected():

    test_gene_data = [
        {"start": 50},
    ]

    test_gene_df = pd.DataFrame(test_gene_data)
    _, density = calculate_gene_density(test_gene_df, 100, 10, 10)

    assert density == [0, 0, 0, 0, 0, 1, 0, 0, 0]

#========================= split_genes_by_strand() tests ================================#

def test_split_genes_by_strand():

    test_gene_data = [
        {"strand":"+"},
        {"strand":"+"},
        {"strand":"-"},
        {"strand":"-"}
    ]

    test_gene_df = pd.DataFrame(test_gene_data)
    f, r = split_genes_by_strand(test_gene_df)

    assert len(f) == 2
    assert len(r) == 2