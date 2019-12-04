#!/usr/bin/env python3

"""Build homology trees."""

import os
from os.path import abspath, basename, expanduser, join, splitext
from datetime import date
import argparse
import textwrap
import pylib.util as util
from pylib.core_homology_trees import build_trees


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

    parser.add_argument(
        '-q', '--quantiles', type=float, default=0.05,
        help="""A TreeShrink only option for tree trimming quantiles. The
            default is 0.05.""")

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
