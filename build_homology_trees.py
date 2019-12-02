#!/usr/bin/env python3

"""Build homology trees."""

import os
from os.path import basename, join, splitext
from glob import glob
from datetime import date
import argparse
import textwrap
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pylib.util as util
import pylib.bio as bio
import pylib.log as log
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
    for path in fasta_paths:
        if args.bootstrap:
            raxml.raxml_bs(args, path, temp_dir)
        else:
            raxml.raxml(args, path, temp_dir)


def get_fasta_paths(args):
    """Get the fasta files to process."""
    pattern = join(args.assemblies_dir, args.file_filter)
    fasta_paths = sorted(glob(pattern))
    if not len(fasta_paths):
        log.fatal('No files were found with this mask: "{}".'.format(pattern))
    return fasta_paths


def long_seq_check(fasta_paths, seq_type):
    """Warn about really long sequences."""
    for path in fasta_paths:
        longest, i = 0, 0

        with open(path) as fasta_file:
            for i, (_, seq) in enumerate(SimpleFastaParser(fasta_file), 1):
                longest = max(longest, len(seq.replace('-', '')))

        if bio.too_long(seq_type, longest):
            log.warn(util.shorten("""{} has {} sequences. 
                The longest is {} characters.
                This is too long and may crash the alignment process.
                """.format(path, i, longest)))


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
        '-b', '--bootstrap',  action='store_true',
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

    args = parser.parse_args()

    util.temp_dir_exists(args.temp_dir)

    if not args.output_prefix:
        prefix = '{}_{}'.format(splitext(
            basename(__file__)), date.today().isoformat())
        args.output_prefix = join('.', prefix)

    return args


if __name__ == "__main__":
    ARGS = parse_args()
    build_trees(ARGS)
