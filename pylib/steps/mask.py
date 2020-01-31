"""Build homology trees."""

import logging
from pylib.wrappers.mask_tips import mask_tips


def mask(args):
    """Mask mono- and paraphyletic-tips that belong to the same taxon."""
    trees = args.cleaned_files
    cleaned = args.tree_files

    for fasta, tree in zip(cleaned, trees):
        logging.info('mask input: {}'.format(fasta))
        logging.info('mask input: {}'.format(tree))
        masked = mask_tips(
            fasta, tree, args.output_dir, args.output_ext)
        logging.info('mask_tips output: {}'.format(masked))
