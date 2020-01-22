"""Build homology trees."""

import logging
import pylib.util as util
from pylib.wrappers.tree_to_fasta import tree_to_fasta


def tree2fa(args):
    """Mask mono- and paraphyletic-tips that belong to the same taxon."""
    trees = util.get_input_files(args.input_dir, args.tree_filter)
    masks = util.get_input_files(args.input_dir, args.masked_filter)

    for fasta, tree in zip(masks, trees):
        logging.info('tree2fa input: {}'.format(tree))
        new_fasta = tree_to_fasta(fasta, tree, args.output_dir)
        logging.info('tree_to_fasta output: {}'.format(new_fasta))
