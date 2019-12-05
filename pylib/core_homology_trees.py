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
from oldlib.mask_tips import mask_tips


def pipeline(args):
    """Build the homology trees."""
    with util.make_temp_dir(
            where=args.temp_dir,
            prefix='{}_'.format(splitext(basename(__file__))[0]),
            keep=args.keep_temp_dir) as temp_dir:
        setup(args, temp_dir)

        for data in get_fasta_files(args):
            try:
                too_few_records(data)
                seq_too_long(data, args)
                fasta_to_tree(data, args)
                tree_shrink(data, args)
                mask_tree(data, args)
                fasta_from_tree(data, args)

            except util.StopProcessing:
                pass


def setup(args, temp_dir):
    """Perform pipeline setup."""
    args.temp_dir = temp_dir  # Put the real temp dir into the args
    log_file = join(args.output_dir, 'build_homology_trees.log')
    log.setup(log_file)


def get_fasta_files(args):
    """Get the fasta data to process."""
    log.info('Gathering fasta files.')
    pattern = join(args.assemblies_dir, args.file_filter)
    fasta_files = sorted([abspath(p) for p in glob(pattern)])
    if len(fasta_files) == 0:
        log.fatal('No data were found with this mask: "{}".'.format(pattern))
    return [{'fasta': f} for f in fasta_files]


def too_few_records(data):
    """Check if the fasta file is too small to make a good tree."""
    log.info('Checking fasta for {}'.format(basename(data['fasta'])))

    if bio.fasta_record_count(data['fasta']) < bio.MIN_SEQ:
        log.warn('"{}" has fewer than {} records, skipping.'.format(
            data['fasta'], bio.MIN_SEQ))
        raise util.StopProcessing()


def seq_too_long(data, args):
    """Warn about really long sequences."""
    longest = bio.longest_fasta_seq(data['fasta'])

    if bio.seqs_too_long(longest, args.seq_type):
        seq_count = bio.fasta_record_count(data['fasta'])
        log.warn(util.shorten("""{} has {} sequences.
            The longest is {} characters.
            This is too long and may crash the alignment process.
            """.format(data['fasta'], seq_count, longest)))


def fasta_to_tree(data, args):
    """Build trees from the fasta data."""
    log.info('Converting fasta to tree for {}'.format(
        basename(data['fasta'])))

    if args.bootstrap:
        data['aligned'] = mafft(data['fasta'], args)
        data['cleaned'] = pxclsq(data['aligned'], args, COLUMN_OCCUPANCY_LG)
        data['tree'] = raxml_bs(data['cleaned'], args)
    elif bio.fasta_record_count(data['fasta']) >= bio.SEQ_COUNT_CUTOFF:
        data['aligned'] = pasta(data['fasta'], args)
        data['cleaned'] = pxclsq(data['aligned'], args, COLUMN_OCCUPANCY_SM)
        data['tree'] = fasttree(data['cleaned'], args)
    else:
        data['aligned'] = mafft(data['fasta'], args)
        data['cleaned'] = pxclsq(data['aligned'], args, COLUMN_OCCUPANCY_SM)
        data['tree'] = raxml(data['cleaned'], args)


def tree_shrink(data, args):
    """Remove long branches from trees."""
    log.info('Shrinking tree for {}'.format(basename(data['fasta'])))

    data['trimmed'] = treeshrink(data['tree'], args)
    data['unrooted'] = pxrr(data['trimmed'], args)


def mask_tree(data, args):
    """Mask mono- and paraphyletic-tips that belong to the same taxon."""
    log.info('Masking tree tips for {}'.format(basename(data['fasta'])))

    data['masked'] = mask_tips(data['cleaned'], data['unrooted'], args)


def fasta_from_tree(data, args):
    """Mask mono- and paraphyletic-tips that belong to the same taxon."""
    log.info('Converting masked tree to fasta for {}'.format(
        basename(data['fasta'])))
