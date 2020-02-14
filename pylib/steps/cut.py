"""Cut long internal branches."""

import logging
from pylib import util
from pylib.wrappers.cut_branches import cut_branches


def cut(args):
    """Cut long internal branches."""
    for tree in args.input_files:
        logging.info('cut_branches input: {}'.format(tree))

        count = util.count_tree_taxa(tree)
        if count < args.min_taxa:
            logging.info(
                '"{}" skipped, it has only {} of {} taxa.'.format(
                    tree, count, args.min_taxa))
            continue

        cut_tree, enough_taxa = cut_branches(
            tree, args.output_dir, args.output_ext,
            args.branch_cutoff, args.min_taxa)

        logging.info('cut_branches output: {}'.format(cut_tree))
