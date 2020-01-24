"""Build homology trees."""

import logging
import pylib.util as util
from pylib.wrappers.phyx import pxrr
from pylib.wrappers.treeshrink import treeshrink


def shrink(args):
    """Remove long branches from trees."""
    for tree in util.get_input_files(args.input_dir, args.input_filter):
        logging.info('shrink input: {}'.format(tree))

        trimmed = treeshrink(
            tree, args.output_dir, args.quantiles, args.output_extension)
        logging.info('treeshrink output: {}'.format(trimmed))

        unrooted = pxrr(trimmed, args.output_dir, args.output_extension)
        logging.info('pxrr output: {}'.format(unrooted))
