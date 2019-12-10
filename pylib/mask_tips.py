"""Mask both mono- and paraphyletic tips that belong to the same taxon."""

import re
from itertools import groupby
from Bio import Phylo
from pylib.util import taxon_id
from . import bio
from . import util


IGNORE = re.compile(r'[x*?\-]', re.IGNORECASE)

MIN_TREE = 4


def mask_tips(fasta_file, tree_file, output_dir):
    """Wrap tree tip removal."""
    tree = Phylo.read(tree_file, 'newick')
    print(tree)
    char_count = {s[0]: len(IGNORE.sub('', s[1]))
                  for s in bio.read_fasta(fasta_file).values()}
    mask_monophyletic_tips(tree, char_count)
    # if args.mask_paraphyletic:
    #     mask_paraphyletic_tips(tree, char_count)

    output = util.file_name(output_dir, tree_file, '.mm')
    Phylo.write(tree, output, 'newick')


def mask_monophyletic_tips(tree, char_count):
    """Mask monophyletic tips."""
    while tree.count_terminals() >= MIN_TREE:
        parents = [n for n in tree.get_nonterminals() if n.is_preterminal()]
        for parent in parents:
            sibs = [s for s in tree.get_terminals() if parent.is_parent_of(s)]
            sibs = sorted(sibs, key=lambda s: s.name)
            for _, group in groupby(sibs, key=lambda s: taxon_id(s.name)):
                group = sorted(
                    group, key=lambda s: char_count[s.name], reverse=True)
                for node in group[1:]:
                    tree.prune(node)


# def mask_paraphyletic_tips(tree, char_count):
#     """Mask paraphyletic tips."""
