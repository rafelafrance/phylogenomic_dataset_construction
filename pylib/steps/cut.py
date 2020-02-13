"""Cut long internal branches."""

import logging
from pylib.wrappers.cut_branches import cut_branches


def cut(args):
    """Cut long internal branches."""
    for tree in args.input_files:
        logging.info('cut_branches input: {}'.format(tree))
        cut = cut_branches(
            tree, args.output_dir, args.output_ext,
            args.branch_cutoff, args.min_taxa)
        logging.info('cut_branches output: {}'.format(cut))
