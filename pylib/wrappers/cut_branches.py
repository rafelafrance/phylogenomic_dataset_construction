"""Cut long internal branches."""

from Bio import Phylo
from pylib import util


def cut_branches(tree_file, output_dir, output_ext, branch_cutoff, min_taxa):
    """Cut long internal branches."""
    tree = Phylo.read(tree_file, 'newick')

    subtrees = cut_deep(tree, branch_cutoff, min_taxa)

    with util.cd(output_dir):
        for i, subtree in enumerate(subtrees):
            output = '{}_{}'.format(tree_file, i)
            output = util.file_name(output, output_ext)
            Phylo.write(subtree, output, 'newick')

    return output


def cut_deep(tree, branch_cutoff, min_taxa):
    """Remove interior nodes if brances are too deep."""
    subtrees = []

    leaves = [n for n in tree.get_terminals()]
    nodes = [n for n in tree.get_nonterminals()]

    for node in nodes:
        children = [n for n in nodes if n in node]

    return subtrees
