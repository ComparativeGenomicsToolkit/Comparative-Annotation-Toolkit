"""
Functions for working with HAL files.
"""
import ete3
from procOps import call_proc_lines


def get_tree(hal):
    """
    Extracts a Tree object from a HAL
    :param hal: HAL file.
    :return: Tree object
    """
    cmd = ['halStats', '--tree', hal]
    newick = call_proc_lines(cmd)[0]
    return ete3.Tree(newick, format=1)


def build_genome_order(hal, ref_genome, genome_subset=None, include_ancestors=False):
    """
    Constructs an ordered vector of genomes based on genetic distance from a source genome.
    :param hal: HAL file.
    :param ref_genome: Genome we are rooting this order on.
    :param genome_subset: List of genomes to have in final set
    :param include_ancestors: Include ancestral genomes?
    :return: list of genomes.
    """
    t = get_tree(hal)
    if len(t) == 1:  # this is a hacked HAL, we cannot determine distances
        return [t.name]
    if include_ancestors is True:
        distances = [[t.get_distance(ref_genome, x), x.name] for x in t.get_descendants() if x.name != ref_genome]
    else:
        distances = [[t.get_distance(ref_genome, x), x.name] for x in t.get_leaves() if x.name != ref_genome]
    if genome_subset is not None:
        distances = [[dist, name] for dist, name in distances if name in genome_subset]
    ordered = sorted(distances, key=lambda (dist, name): dist)
    distances, ordered_names = zip(*ordered)
    return ordered_names


def extract_genomes(hal, include_ancestors=False, target_genomes=None):
    """
    Constructs a set of genomes present in this alignment
    :param hal: HAL file.
    :param include_ancestors: Should we also extract ancestral genomes?
    :param target_genomes: Subset of genomes to target.
    :return: tuple of genomes.
    """
    t = get_tree(hal)
    if include_ancestors is True:
        if target_genomes is None:
            return tuple(x.name for x in t.get_descendants())
        else:
            anc = t.get_common_ancestor(target_genomes)
            return tuple(x.name for x in anc.get_descendants() if len(x.get_descendants()) != 0 or x.name in target_genomes)
    else:
        if target_genomes is None:
            return tuple(t.get_leaf_names())
        else:
            return tuple(x for x in t.get_leaf_names() if x in target_genomes)
