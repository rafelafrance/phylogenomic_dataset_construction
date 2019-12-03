#!/usr/bin/env python3

"""Build homology trees."""

import os
from os.path import abspath, basename, expanduser, join, splitext
from glob import glob
from datetime import date
import argparse
import textwrap
import pylib.util as util
import pylib.bio as bio
import pylib.log as log
from pylib.mafft import mafft
import pylib.raxml as raxml


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
        fasta_to_tree(args, fasta_paths, temp_dir)


def fasta_to_tree(args, fasta_paths, temp_dir):
    """Build trees from the fasta files."""
    for fasta in fasta_paths:
        seq_count = bio.fasta_record_count(fasta)

        if seq_count < bio.MIN_SEQ:
            log.warn('"{}" has fewer than {} records, skipping.'.format(
                fasta, bio.MIN_SEQ))
            continue

        if args.bootstrap:
            # alignment = mafft(path, fasta, num_cores, seq_type)
            # cleaned = pxclsq(path, alignment, COLUMN_OCCUPANCY_LG, seq_type)
            # raxml_bs(path, cleaned, num_cores, seq_type)
            raxml.raxml_bs(args, fasta, temp_dir)
        elif seq_count >= bio.SEQ_COUNT_CUTOFF:
            pass
            # alignment = pasta(path, fasta, num_cores, seq_type)
            # cleaned = pxclsq(path, alignment, COLUMN_OCCUPANCY_SM seq_type)
            # fasttree(path, cleaned, seq_type)
        else:
            alignment = mafft(args, fasta, temp_dir)
            # cleaned = pxclsq(path, alignment, COLUMN_OCCUPANCY_SM seq_type)
            # raxml(path, cleaned, num_cores, seq_type)
            raxml.raxml(args, fasta, temp_dir)


def get_fasta_paths(args):
    """Get the fasta files to process."""
    pattern = join(args.assemblies_dir, args.file_filter)
    fasta_paths = sorted([abspath(p) for p in glob(pattern)])
    for path in fasta_paths:
        print(path)
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


def parse_args():
    """Process command-line arguments."""
    description = """Align each cluster, trim alignment, and infer a tree."""
    parser = argparse.ArgumentParser(
        allow_abbrev=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(description),
        fromfile_prefix_chars='@')

    parser.add_argument(
        '--version', '-V', action='version',
        version='%(prog)s v{}'.format(util.__VERSION__))

    parser.add_argument(
        '-a', '--assemblies-dir', metavar='PATH', required=True,
        help="""The path to the DNA contigs.""")

    parser.add_argument(
        '-f', '--file-filter', default='*.fasta', metavar='FILTER',
        help="""Use this to filter files in the assemblies directory. For
            example '*filtered*.fasta' will select all fasta files in the
            assemblies directory with the word filtered in them. The default
            is to select all fasta files in the assemblies directory
            '*.fasta'.""")

    cpus = min(10, os.cpu_count() - 4 if os.cpu_count() > 4 else 1)
    parser.add_argument(
        '-c', '--cpus', '--processes', type=int, default=cpus,
        help="""Number of CPU processors to use. This is passed to wrapped
            programs that use this option like: mafft.
            The default will use {} out of {} CPUs.
            """.format(cpus, os.cpu_count()))

    parser.add_argument(
        '-t', '--seq-type', default='dna',
        choices=['dna', 'aa'],
        help="""Are building trees from DNA or amino acid sequences. The
            default is "dna".""")

    parser.add_argument(
        '-b', '--bootstrap', action='store_true',
        help="""Turn on rapid bootstrapping.""")

    parser.add_argument(
        '-s', '--seed', type=int, default=12345,
        help="""A random number seed. This allows you to reproduce your
            results and helps with debugging the program.""")

    parser.add_argument(
        '-o', '--output-prefix',
        help="""This is the prefix of all of the output files. So you can
            identify different output file sets. You may include a directory
            as part of the prefix. This program will add suffixes to
            differentiate output files.""")

    parser.add_argument(
        '--temp-dir', metavar='DIR',
        help="""Place temporary files in this directory. All files will be
            deleted after aTRAM completes. The directory must exist.""")

    parser.add_argument(
        '--keep-temp-dir', action='store_true',
        help="""This flag will keep the temporary files in the --temp-dir
            around for debugging.""")

    parser.add_argument(
        '--anysymbol', action='store_true',
        help="""A mafft only option to handle when there are "U"s in aa
            sequences.""")

    args = parser.parse_args()

    util.temp_dir_exists(args.temp_dir)

    if not args.output_prefix:
        prefix = '{}_{}'.format(splitext(
            basename(__file__)), date.today().isoformat())
        args.output_prefix = join('.', prefix)
    args.output_prefix = abspath(expanduser(args.output_prefix))

    return args


if __name__ == "__main__":
    ARGS = parse_args()
    build_trees(ARGS)
