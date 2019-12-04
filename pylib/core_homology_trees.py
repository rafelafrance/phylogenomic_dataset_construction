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


def build_trees(args):
    """Build the homology trees."""
    with util.make_temp_dir(
            where=args.temp_dir,
            prefix='{}_'.format(splitext(basename(__file__))[0]),
            keep=args.keep_temp_dir) as temp_dir:
        log_file = join(args.output_prefix, 'log_file.txt')
        log.setup(log_file)
        fasta_paths = get_fasta_paths(args)
        long_seq_check(fasta_paths, args.seq_type)
        trees = fasta_to_trees(args, fasta_paths, temp_dir)
        trees = tree_shrink(args, trees, temp_dir)


def mask_tips(args, trees, temp_dir):
    """Mask mono and paraphyletic tips that belong to the same taxon."""


def tree_shrink(args, old_trees, temp_dir):
    """Remove long branches from trees."""
    trees = []
    for tree in old_trees:
        tree = treeshrink(args, tree, temp_dir)
        trees.append(pxrr(args, tree, temp_dir))
    return trees


def fasta_to_trees(args, fasta_paths, temp_dir):
    """Build trees from the fasta files."""
    trees = []
    for fasta in fasta_paths:
        seq_count = bio.fasta_record_count(fasta)

        if seq_count < bio.MIN_SEQ:
            log.warn('"{}" has fewer than {} records, skipping.'.format(
                fasta, bio.MIN_SEQ))
            continue

        if args.bootstrap:
            alignment = mafft(args, fasta, temp_dir)
            cleaned = pxclsq(args, alignment, temp_dir, COLUMN_OCCUPANCY_LG)
            trees.append(raxml_bs(args, cleaned, temp_dir))
        elif seq_count >= bio.SEQ_COUNT_CUTOFF:
            alignment = pasta(args, fasta, temp_dir)
            cleaned = pxclsq(args, alignment, temp_dir, COLUMN_OCCUPANCY_SM)
            trees.append(fasttree(args, cleaned, temp_dir))
        else:
            alignment = mafft(args, fasta, temp_dir)
            cleaned = pxclsq(args, alignment, temp_dir, COLUMN_OCCUPANCY_SM)
            trees.append(raxml(args, cleaned, temp_dir))

    return trees


def get_fasta_paths(args):
    """Get the fasta files to process."""
    pattern = join(args.assemblies_dir, args.file_filter)
    fasta_paths = sorted([abspath(p) for p in glob(pattern)])
    if len(fasta_paths) == 0:
        log.fatal('No files were found with this mask: "{}".'.format(pattern))
    return fasta_paths


def long_seq_check(fasta_paths, seq_type):
    """Warn about really long sequences."""
    for fasta_path in fasta_paths:
        longest = bio.longest_fasta_seq(fasta_path)
        if bio.too_long(longest, seq_type):
            seq_count = bio.fasta_record_count(fasta_path)
            log.warn(util.shorten("""{} has {} sequences.
                The longest is {} characters.
                This is too long and may crash the alignment process.
                """.format(fasta_path, seq_count, longest)))
