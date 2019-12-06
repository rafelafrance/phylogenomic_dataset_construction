#!/usr/bin/env python3

"""Build homology trees."""

import os
from os.path import abspath, expanduser
import argparse
import textwrap
import pylib.util as util
from pylib.core_homology_trees import pipeline


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
        '-i', '--input-dir', metavar='PATH', required=True,
        help="""The directory containing the input fasta files.""")

    parser.add_argument(
        '-f', '--input-filter', default='*.fasta', metavar='FILTER',
        help="""Use this to filter files in the input directory. For
            example '*filtered*.fasta' will select all fasta files in the
            input directory with the word filtered in them. The default
            is to select all fasta files in the input directory
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
        '-o', '--output-dir', default='.',
        help="""Place output files in this directory. The default is the
            current directory.""")

    parser.add_argument(
        '--anysymbol', action='store_true',
        help="""A mafft only option to handle when there are "U"s in aa
            sequences.""")

    parser.add_argument(
        '-q', '--quantiles', type=float, default=0.05,
        help="""A TreeShrink only option for tree trimming quantiles. The
            default is 0.05.""")

    parser.add_argument(
        '--mask-paraphyletic', action='store_true',
        help="""When masking tree tips, do you want to also mask paraphyletic
            tips.""")

    parser.add_argument(
        '--temp-dir', metavar='DIR',
        help="""Place temporary files in this directory. All files will be
            deleted after aTRAM completes. The directory must exist.""")

    parser.add_argument(
        '--keep-temp-dir', action='store_true',
        help="""This flag will keep the temporary files in the --temp-dir
            around for debugging.""")

    args = parser.parse_args()

    args.output_dir = abspath(expanduser(args.output_dir))
    util.temp_dir_exists(args.temp_dir)

    return args


if __name__ == "__main__":
    ARGS = parse_args()
    pipeline(ARGS)
