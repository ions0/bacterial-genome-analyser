"""
Bacterial Genome Analyser: tests/test_genome_utils.py

"""
from genome_utils import find_gene_by_name, categorise_genes_by_size

def test_find_gene_by_name():

    test_gene_data = [
        {"gene": "lacZ"}
    ]

    test_dict = find_gene_by_name(test_gene_data, "lacZ")

    assert test_dict["gene"] == "lacZ"
    assert find_gene_by_name(test_gene_data, "lacX") == None

def test_catagorise_gene_by_size():

    test_gene_data = [
        {"length": 200},
        {"length": 700},
        {"length": 1900},
        {"length": 3000},
        {"length": 5000},
    ]

    test_cat_dict = categorise_genes_by_size(test_gene_data)

    assert len(test_cat_dict["tiny_genes"]) == 1
    assert len(test_cat_dict["small_genes"]) == 1
    assert len(test_cat_dict["med_genes"]) == 1
    assert len(test_cat_dict["large_genes"]) == 1
    assert len(test_cat_dict["huge_genes"]) == 1