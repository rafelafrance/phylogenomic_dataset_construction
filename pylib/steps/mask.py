"""Build homology trees."""

import logging
import pylib.util as util
import pylib.bio as bio
from pylib.wrappers.mask_tips import mask_tips


def mask_tree(args):
    """Mask mono- and paraphyletic-tips that belong to the same taxon."""
    logging.info('Masking tree tips for {}'.format(args.log_name))

    mask_tips(cleaned, unrooted, args.output_dir, args.mask_paraphyletic)
