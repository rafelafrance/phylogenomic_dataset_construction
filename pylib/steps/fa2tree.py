"""Build homology trees."""

from os.path import abspath
import logging
import pylib.util as util
import pylib.bio as bio
from pylib.wrappers.mafft import mafft
from pylib.wrappers.raxml import raxml, raxml_bs
from pylib.wrappers.phyx import pxclsq
from pylib.wrappers.pasta import pasta
from pylib.wrappers.fasttree import fasttree


def fa2tree(args):
    """Build trees from the fasta data."""
    for fasta in util.get_input_files(args.input_dir, args.input_filter):
        logging.info('fa2tree input: {}'.format(fasta))
        fasta = abspath(fasta)
        if args.bootstrap:
            fa2tree_bs(args, fasta)
        elif bio.fasta_record_count(fasta) >= bio.SEQ_COUNT_CUTOFF:
            fa2tree_big(args, fasta)
        else:
            fa2tree_default(args, fasta)


def fa2tree_bs(args, fasta):
    """Build trees from the fasta data, bootstrap version."""
    aligned = mafft(
        fasta, args.output_dir, args.seq_type, args.cpus, args.anysymbol,
        args.output_extension)
    logging.info('mafft output: {}'.format(aligned))

    cleaned = pxclsq(
        aligned, args.output_dir, args.seq_type, args.min_occupancy,
        args.min_seq_len, args.output_extension)
    logging.info('pxclsq output: {}'.format(cleaned))

    tree = raxml_bs(
        cleaned, args.output_dir, args.seq_type, args.cpus, args.seed,
        args.output_extension)
    logging.info('raxml_bs output: {}'.format(tree))


def fa2tree_big(args, fasta):
    """Build trees from the fasta data, large file version."""
    aligned = pasta(
        fasta, args.output_dir, args.seq_type, args.output_extension,
        args.cpus)
    logging.info('pasta output: {}'.format(aligned))

    cleaned = pxclsq(
        aligned, args.output_dir, args.seq_type, args.min_occupancy,
        args.min_seq_len, args.output_extension)
    logging.info('pxclsq output: {}'.format(cleaned))

    tree = fasttree(
        cleaned, args.output_dir, args.seq_type, args.output_extension)
    logging.info('fasttree output: {}'.format(tree))


def fa2tree_default(args, fasta):
    """Build trees from the fasta data, normal version."""
    aligned = mafft(
        fasta, args.output_dir, args.seq_type, args.cpus, args.anysymbol,
        args.output_extension)
    logging.info('mafft output: {}'.format(aligned))

    cleaned = pxclsq(
        aligned, args.output_dir, args.seq_type, args.min_occupancy,
        args.min_seq_len, args.output_extension)
    logging.info('pxclsq output: {}'.format(cleaned))

    tree = raxml(
        cleaned, args.output_dir, args.seq_type, args.cpus, args.seed,
        args.output_extension)
    logging.info('raxml output: {}'.format(tree))
