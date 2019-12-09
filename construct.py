#!/usr/bin/env python3

"""Build homology trees."""

import os
from os.path import abspath, expanduser
import argparse
import textwrap
import pylib.util as util
from pylib.core_construct import pipeline


def parse_args():
    """Process command-line arguments."""
    # python scripts/prune_paralogs_MI.py data/locus_001/fasta/ tt 0.02
    # 0.02 2 data/locus_001/
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=True,
        fromfile_prefix_chars='@',
        description=textwrap.dedent("""Phylogenomic dataset construction."""),
        epilog=textwrap.dedent("""
            Pruning options (--prune)
            -------------------------
              1to1: Only look at homologs that are strictly one-to-one.
                    No cutting is carried out.
              mi:   Prune by using homologs with monophyletic, non-repeating
                    out-groups, reroot and cut paralog from root to tip.
              mo:   Prune by using homologs with monophyletic, non-repeating
                    out-groups reroot and cut paralog from root to tip.
              rt:   rune by extracting in-group clades and then cut paralogs
                    from root to tip. If no out-group, only use those that do
                    not have duplicated taxa.
            """))

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

    parser.add_argument(
        '-o', '--output-dir', default='.',
        help="""Place output files in this directory. The default is the
            current directory.""")

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
        help="""Are we building trees from DNA or amino acid sequences. The
            default is "dna".""")

    parser.add_argument(
        '-b', '--bootstrap', action='store_true',
        help="""Turn on rapid bootstrapping.""")

    parser.add_argument(
        '-s', '--seed', type=int, default=12345,
        help="""A random number seed. This allows you to reproduce your
            results and helps with debugging the program.""")

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
        '-p', '--prune', choices=['1to1', 'mi', 'mo', 'rt'], required=True,
        help="""How will you prune the input trees. See below for what each
            option means.""")

    parser.add_argument(
        '-m', '--min-taxa', type=int, default=2,
        help="""""")

    parser.add_argument(
        '-a', '--absolute-tip-cutoff', type=float, default=0.02,
        help="""""")

    parser.add_argument(
        '-r', '--relative-tip-cutoff', type=float, default=0.02,
        help="""""")

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
