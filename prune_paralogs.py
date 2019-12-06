#!/usr/bin/env python3

"""Prune paralogs to infer homologs."""

from os.path import abspath, expanduser
import argparse
import textwrap
import pylib.util as util
from pylib.core_prune_paralogs import pipeline


def parse_args():
    """Process command-line arguments."""
    # python scripts/prune_paralogs_MI.py data/locus_001/fasta/ tt 0.02
    # 0.02 2 data/locus_001/
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=True,
        fromfile_prefix_chars='@',
        description=textwrap.dedent("""Prune paralogs to infer homologs."""),
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
        help="""The directory containing the input Newick tree files.""")

    parser.add_argument(
        '-f', '--input-filter', default='*.tt', metavar='FILTER',
        help="""Use this to filter tree files in the input directory. For
            example '*filtered*.tt' will select all 'tt'' files in the
            input directory with the word filtered in them. The default
            is to select all '*.tt' files in the input directory.""")

    parser.add_argument(
        '-o', '--output-dir', default='.',
        help="""Place output files in this directory. The default is the
            current directory.""")

    parser.add_argument(
        '-p', '--prune', choices=['1to1', 'mi', 'mo', 'rt'], required=True,
        help="""How will you prune the input trees. See below for what each
            option means.""")

    parser.add_argument(
        '-m', '--minimum-taxa', type=int, default=2,
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
