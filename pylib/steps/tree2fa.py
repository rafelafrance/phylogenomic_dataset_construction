"""Build homology trees."""

import logging
from pylib.wrappers.tree_to_fasta import tree_to_fasta


def tree2fa(args):
    """Mask mono- and paraphyletic-tips that belong to the same taxon."""
    trees = args.tree_files
    masks = args.mask_files

    for fasta, tree in zip(masks, trees):
        logging.info('tree2fa input: {}'.format(tree))
        new_fasta = tree_to_fasta(
            fasta, tree, args.output_dir, args.output_ext)
        logging.info('tree_to_fasta output: {}'.format(new_fasta))
