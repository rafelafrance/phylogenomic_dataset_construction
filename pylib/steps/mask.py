"""Build homology trees."""

import logging
import pylib.util as util
from pylib.wrappers.mask_tips import mask_tips


def mask(args):
    """Mask mono- and paraphyletic-tips that belong to the same taxon."""
    trees = util.get_input_files(args.input_dir, args.tree_filter)
    cleaned = util.get_input_files(args.input_dir, args.clean_filter)

    for fasta, tree in zip(cleaned, trees):
        logging.info('mask input: {}'.format(fasta))
        logging.info('mask input: {}'.format(tree))
        masked = mask_tips(
            fasta, tree, args.output_dir, args.output_extension)
        logging.info('mask_tips output: {}'.format(masked))
