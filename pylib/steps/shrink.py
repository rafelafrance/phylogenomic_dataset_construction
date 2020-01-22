"""Build homology trees."""

import logging
import pylib.util as util
from pylib.wrappers.phyx import pxrr
from pylib.wrappers.treeshrink import treeshrink


def shrink(args):
    """Remove long branches from trees."""
    for tree in util.get_input_files(args):
        logging.info('Step shrink: {}'.format(tree))

        trimmed = treeshrink(tree, args.output_dir, args.quantiles)
        logging.info('treeshrink: {}'.format(trimmed))

        unrooted = pxrr(trimmed, args.output_dir)
        logging.info('pxrr: {}'.format(unrooted))
