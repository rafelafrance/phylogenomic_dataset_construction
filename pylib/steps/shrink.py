"""Build homology trees."""

from os.path import abspath
import logging
from pylib.wrappers.phyx import pxrr
from pylib.wrappers.treeshrink import treeshrink


def shrink(args):
    """Remove long branches from trees."""
    for tree in args.input_files:
        logging.info('shrink input: {}'.format(tree))

        tree = abspath(tree)

        logging.info('treeshrink started')
        trimmed = treeshrink(tree, args.output_dir, args.output_ext,
                             args.quantiles)
        logging.info('treeshrink output: {}'.format(trimmed))

        logging.info('pxrr started')
        unrooted = pxrr(trimmed, args.output_dir)
        logging.info('pxrr output: {}'.format(unrooted))
