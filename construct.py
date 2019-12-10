#!/usr/bin/env python3

"""Build homology trees."""

import re
import sys
import os
from os.path import abspath, exists, expanduser
import argparse
import textwrap
from pylib import util
from pylib.core_construct import pipeline


def parse_args():
    """Process command-line arguments."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=True,
        fromfile_prefix_chars='@',
        description=textwrap.dedent(util.__TITLE__))

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

    parser.add_argument(
        '-p', '--prune', choices=['1to1', 'mi', 'mo', 'rt'], required=True,
        help="""How will you prune the input trees.
              "1to1" = Only look at homologs that are strictly one-to-one.
                       No cutting is carried out.
              "mi" = Prune by using homologs with monophyletic, non-repeating
                     out-groups, reroot and cut paralog from root to tip.
              "mo" = Prune by using homologs with monophyletic, non-repeating
                     out-groups reroot and cut paralog from root to tip.
              "rt" = Prune by extracting in-group clades and then cut paralogs
                     from root to tip. If no out-group, only use those that do
                     not have duplicated taxa.""")

    parser.add_argument(
        '-t', '--seq-type', default='dna',
        choices=['dna', 'aa'],
        help="""Are we building trees from DNA or amino acid sequences. The
            default is "dna".""")

    cpus = min(10, os.cpu_count() - 4 if os.cpu_count() > 4 else 1)
    parser.add_argument(
        '--cpus', '--processes', type=int, default=cpus,
        help="""Number of CPU processors to use. This is passed to wrapped
            programs that use this option like: mafft.
            The default will use {} out of {} CPUs.
            """.format(cpus, os.cpu_count()))

    parser.add_argument(
        '--bootstrap', action='store_true',
        help="""Turn on rapid bootstrapping.""")

    parser.add_argument(
        '--min-taxa', type=int, default=2,
        help="""""")

    parser.add_argument(
        '--min-occupancy', type=float, default=0.3,
        help="""""")

    parser.add_argument(
        '--min-seq-len', type=int, default=10,
        help="""""")

    parser.add_argument(
        '--seed', type=int, default=12345,
        help="""A random number seed. This allows you to reproduce your
            results and helps with debugging the program.""")

    parser.add_argument(
        '--anysymbol', action='store_true',
        help="""A mafft only option to handle when there are "U"s in aa
            sequences.""")

    parser.add_argument(
        '--quantiles', type=float, default=0.05,
        help="""A TreeShrink only option for tree trimming quantiles. The
            default is 0.05.""")

    parser.add_argument(
        '--mask-paraphyletic', action='store_true',
        help="""When masking tree tips, do you want to also mask paraphyletic
            tips.""")

    parser.add_argument(
        '--min-bootstrap', type=float, default=0.0,
        help="""""")

    parser.add_argument(
        '--absolute-tip-cutoff', type=float, default=0.02,
        help="""""")

    parser.add_argument(
        '--relative-tip-cutoff', type=float, default=0.02,
        help="""""")

    parser.add_argument(
        '--out-groups',
        help="""This is a comma separated list of out-groups used while
            pruning trees. You may need to quote this argument.""")

    parser.add_argument(
        '--taxon-code-file', metavar='CODE-FILE',
        help="""Path to the taxon code file.""")

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

    return args


def check_args(args):
    """Check arguments are consistent."""
    if args.temp_dir and not exists(args.temp_dir):
        sys.exit('The temporary directory must exist.')

    if args.prune == 'mo' and not args.in_groups:
        sys.exit(util.shorten("""You must specify out-groups
            when --prune=mo."""))

    if args.prune == 'rt' and not args.taxon_code_file:
        sys.exit(util.shorten("""You must specify a taxon code file
            when --prune=rt."""))


def parse_out_groups(args):
    """Check all sequences are a member of an in- or out-group."""
    if args.prune != 'mo':
        return

    args.out_groups = args.out_groups.strip(r'\'"')
    args.out_groups = [g.strip() for g
                       in re.split(r'\s*[\'",]+\s*', args.out_groups)]
    args.out_groups = [g for g in args.out_groups if g]


if __name__ == "__main__":
    ARGS = parse_args()
    check_args(ARGS)
    parse_out_groups(ARGS)
    pipeline(ARGS)
