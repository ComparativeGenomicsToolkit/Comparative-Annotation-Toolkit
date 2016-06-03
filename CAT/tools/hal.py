"""
Functions for working with HAL files.
"""
import ete3
from procOps import call_proc_lines


def build_genome_order(hal, ref_genome):
    """
    Constructs an ordered vector of genomes based on genetic distance from a source genome.
    :param hal: HAL file.
    :param ref_genome: Genome we are rooting this order on.
    :return: list of genomes.
    """
    cmd = ['halStats', '--tree', hal]
    newick = call_proc_lines(cmd)[0]
    t = ete3.Tree(newick, format=1)
    distances = [[t.get_distance(ref_genome, x), x.name] for x in t if x.name != ref_genome]
    ordered = sorted(distances, key=lambda (dist, name): dist)
    distances, ordered_names = zip(*ordered)
    return ordered_names


def extract_genomes(hal):
    """
    Constructs a set of genomes present in this alignment
    :param hal: HAL file.
    :return: list of genomes.
    """
    cmd = ['halStats', '--tree', hal]
    newick = call_proc_lines(cmd)[0]
    t = ete3.Tree(newick, format=1)
    return t.get_leaf_names()
