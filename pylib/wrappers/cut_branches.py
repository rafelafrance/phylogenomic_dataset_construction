"""Cut long internal branches."""

from Bio import Phylo
from pylib import util


def cut_branches(tree_file, output_dir, output_ext, branch_cutoff, min_taxa):
    """Cut long internal branches."""
    tree = Phylo.read(tree_file, 'newick')

    output = util.file_name(tree_file, output_ext)
    with util.cd(output_dir):
        Phylo.write(tree, output, 'newick')

    return output
