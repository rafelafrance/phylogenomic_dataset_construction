"""Mask both mono- and paraphyletic tips that belong to the same taxon."""

import re
from Bio import Phylo
import pylib.bio as bio


WILDCARDS = re.compile(r'[x*?\-]', re.IGNORECASE)

MIN_TREE = 4


def mask_tips(fasta_file, tree_file, args):
    """Wrap tree tip removal."""
    print(tree_file)
    tree = Phylo.read(tree_file, 'newick')
    print(dir(tree))
    chars = remove_wildcards(fasta_file)
    mask_monophyletic_tips(tree, chars)
    print(tree)


def remove_wildcards(fasta_file):
    """Get all non-wildcard characters in a sequence."""
    return {s[0]: WILDCARDS.sub('', s[1]) for s in bio.read_fasta(fasta_file)}


def mask_monophyletic_tips(tree, chars):
    """Mask monophyletic tips."""
    going = True
    while going and tree.count_terminals() >= MIN_TREE:
        going = False
        parents = [n for n in tree.get_nonterminals() if n.is_preterminal()]
        for parent in parents:
            print(parent)
            sibs = [s for s in tree.get_terminals() if parent.is_parent_of(s)]
            print(sibs)
            print()


def mask_paraphyletic_tips(tree, chars):
    """Mask paraphyletic tips."""
