"""Mask both mono- and paraphyletic tips that belong to the same taxon."""

import re
from itertools import groupby
from Bio import Phylo
from pylib.util import taxon_id
from pylib import bio
from pylib import util


IGNORE = re.compile(r'[x*?\-]', re.IGNORECASE)

MIN_TREE = 4


def mask_tips(tree_file, output_dir, output_ext):
    """Wrap tree tip removal."""
    tree = Phylo.read(tree_file, 'newick')

    mask_monophyletic_tips(tree)

    output = util.file_name(tree_file, output_ext)
    with util.cd(output_dir):
        Phylo.write(tree, output, 'newick')

    return output


def mask_monophyletic_tips(tree):
    """Mask monophyletic tips."""
    again = True
    while again and tree.count_terminals() >= MIN_TREE:
        again = False
        parents = [n for n in tree.get_nonterminals() if n.is_preterminal()]
        for parent in parents:
            sibs = [s for s in tree.get_terminals() if parent.is_parent_of(s)]
            sibs = sorted(sibs, key=lambda s: s.branch_length)
            for node in sibs[1:]:
                tree.prune(node)
                again = True
