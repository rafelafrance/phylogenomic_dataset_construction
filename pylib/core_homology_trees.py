#!/usr/bin/env python3

"""Build homology trees."""

from os.path import abspath, basename, join, splitext
from glob import glob
import pylib.util as util
import pylib.bio as bio
import pylib.log as log
from pylib.mafft import mafft
from pylib.raxml import raxml, raxml_bs
from pylib.phyx import pxclsq, pxrr, COLUMN_OCCUPANCY_LG, COLUMN_OCCUPANCY_SM
from pylib.pasta import pasta
from pylib.fasttree import fasttree
from pylib.treeshrink import treeshrink
from pylib.mask_tips import mask_tips


def pipeline(args):
    """Build the homology trees."""
    with util.make_temp_dir(
            where=args.temp_dir,
            prefix='{}_'.format(splitext(basename(__file__))[0]),
            keep=args.keep_temp_dir) as temp_dir:
        setup(args, temp_dir)

        for files in get_fasta_files(args):
            print(files)
            try:
                too_few_records(files)
                seqs_too_long(files, args)
                fasta_to_tree(files, args)
                tree_shrink(files, args)
                mask_tree(files, args)
                fasta_from_tree(files, args)

            except util.StopProcessing:
                pass


def setup(args, temp_dir):
    """Perform pipeline setup."""
    args.temp_dir = temp_dir  # Put the real temp dir into the args
    log_file = join(args.output_dir, 'build_homology_trees.log')
    log.setup(log_file)


def get_fasta_files(args):
    """Get the fasta files to process."""
    pattern = join(args.assemblies_dir, args.file_filter)
    fasta_files = sorted([abspath(p) for p in glob(pattern)])
    if len(fasta_files) == 0:
        log.fatal('No files were found with this mask: "{}".'.format(pattern))
    return [{'fasta': f} for f in fasta_files]


def too_few_records(files):
    """Check if the fasta file is too small to make a good tree."""
    if bio.fasta_record_count(files['fasta']) < bio.MIN_SEQ:
        log.warn('"{}" has fewer than {} records, skipping.'.format(
            files['fasta'], bio.MIN_SEQ))
        raise util.StopProcessing()


def seqs_too_long(files, args):
    """Warn about really long sequences."""
    longest = bio.longest_fasta_seq(files['fasta'])
    if bio.seqs_too_long(longest, args.seq_type):
        seq_count = bio.fasta_record_count(files['fasta'])
        log.warn(util.shorten("""{} has {} sequences.
            The longest is {} characters.
            This is too long and may crash the alignment process.
            """.format(files['fasta'], seq_count, longest)))


def fasta_to_tree(files, args):
    """Build trees from the fasta files."""
    if args.bootstrap:
        files['alignment'] = mafft(files['fasta'], args)
        files['pxclsq'] = pxclsq(files['alignment'], args, COLUMN_OCCUPANCY_LG)
        files['tree'] = raxml_bs(files['pxclsq'], args)
    elif bio.fasta_record_count(files['fasta']) >= bio.SEQ_COUNT_CUTOFF:
        files['alignment'] = pasta(files['fasta'], args)
        files['pxclsq'] = pxclsq(files['alignment'], args, COLUMN_OCCUPANCY_SM)
        files['tree'] = fasttree(files['pxclsq'], args)
    else:
        files['alignment'] = mafft(files['fasta'], args)
        files['pxclsq'] = pxclsq(files['alignment'], args, COLUMN_OCCUPANCY_SM)
        files['tree'] = raxml(files['pxclsq'], args)


def tree_shrink(files, args):
    """Remove long branches from trees."""
    files['treeshrink'] = treeshrink(files['tree'], args)
    files['pxrr'] = pxrr(files['treeshrink'], args)


def mask_tree(files, args):
    """Mask mono- and paraphyletic-tips that belong to the same taxon."""
    mask_tips(files['pxclsq'], files['pxrr'], args)


def fasta_from_tree(files, args):
    """Mask mono- and paraphyletic-tips that belong to the same taxon."""
