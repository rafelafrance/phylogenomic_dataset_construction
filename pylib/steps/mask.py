"""Build homology trees."""

import logging
from pylib.wrappers.mask_tips import mask_tips


def mask(args):
    """Mask monophyletic tree tips that belong to the same taxon."""
    for tree in args.input_files:
        logging.info('mask_tips input: {}'.format(tree))
        masked = mask_tips(tree, args.output_dir, args.output_ext)
        logging.info('mask_tips output: {}'.format(masked))
