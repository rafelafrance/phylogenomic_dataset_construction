#!/usr/bin/env python3

"""Build homology trees."""

import re
import sys
import os
from os.path import abspath, expanduser, isfile
import logging
from glob import glob
import argparse
from pylib import util
from pylib import bio
from pylib.steps.check import check
from pylib.steps.fa2tree import fa2tree
from pylib.steps.shrink import shrink
from pylib.steps.mask import mask
from pylib.steps.tree2fa import tree2fa
from pylib.steps.prune import prune_paralogs
from pylib.steps.orth2fa import orthologs_to_fasta


STEP = 0
INPUT_ATTRS = set()


def construct():
    """The main function."""
    formatter = '%(asctime)s %(levelname)s: %(message)s'
    logging.basicConfig(
        level=logging.INFO, format=formatter,  datefmt='%Y-%m-%d %H:%M:%S')

    args = parse_args()

    step_name = args.func.__name__  # Get entered step via its function name
    if step_name == 'prune':
        check_args(args)
        parse_out_groups(args)

    args.func(args)


def parse_args():
    """Process command-line arguments and run the chose function."""

    description = util.shorten("""
        Infer orthologies for phylogenomic analyses.""")

    epilog = util.shorten("""
        To get help on a particular step type the step name and then -h.
        For example to get help on step 1 you can type:
        "./construct.py fa2tree -h".""")

    usage = """{} [-h] STEP""".format(__file__)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=True,
        fromfile_prefix_chars='@',
        usage=usage,
        description=description,
        epilog=epilog)

    subparsers = parser.add_subparsers()
    check_step(subparsers)
    fasta2tree_step(subparsers)
    shrink_step(subparsers)
    mask_step(subparsers)
    tree2fa_step(subparsers)
    prune_step(subparsers)
    orth2fa_step(subparsers)

    args = parser.parse_args()

    expand_files(args)

    return args


def expand_files(args):
    """Handle unexpanded globs."""
    global INPUT_ATTRS

    for attr in INPUT_ATTRS:
        if not hasattr(args, attr):
            continue

        names = getattr(args, attr)
        names = names if isinstance(names, list) else [names]
        setattr(args, attr, names)

        arg_files = []
        try:
            for name in getattr(args, attr):
                if isfile(name):
                    arg_files += [name]
                elif any(name.find(x) > -1 for x in ('*', '?', '[')):
                    files = glob(name)
                    if not files:
                        raise ValueError(name)
                    arg_files += files
                else:
                    raise ValueError(name)
        except ValueError as err:
            logging.critical('"{}" Did not match files'.format(err))
            sys.exit(1)
        setattr(args, attr, arg_files)


def check_step(subparsers):
    """Add check step."""
    check_parser = subparsers.add_parser(
        'check', help=helper("""Check that this utility can handle the input
            fasta files. Are there too few fasta records in the file to make
            a tree? Or are there really long sequences that may crash the
            alignment process?"""))
    input_files(check_parser, './*.fasta')
    seq_type_arg(check_parser)
    check_parser.set_defaults(func=check)


def fasta2tree_step(subparsers):
    """Add fa2tree step."""
    fa2tree_parser = subparsers.add_parser(
        'fa2tree', help=helper("""Build homolog trees from fasta files. This
            performs three steps: 1) It will align the sequences. 2) It will
            clean the aligned sequences. 3) Create a tree from the cleaned and
            aligned sequences. Normally, this program will use "mafft" to
            align the sequences, "pxclsq" to clean them, and "raxml" to build
            the tree. If you select the "--bootstrap" option it still uses
            "mafft" and "pxclsq" for aligning and cleaning but it uses
            "raxml_bs" for tree building. And it the fasta file has more than
            {} sequences then it will use "pasta" for alignment, "pxclsq" for
            cleaning and "fasttree" for tree building.""".format(
            bio.SEQ_COUNT_CUTOFF)))
    io_args(fa2tree_parser, '*.fasta', '.tre')
    seq_type_arg(fa2tree_parser)
    cpus_arg(fa2tree_parser)
    fa2tree_parser.add_argument(
        '--bootstrap', action='store_true',
        help="""Turn on rapid bootstrapping.""")
    fa2tree_parser.add_argument(
        '--anysymbol', action='store_true',
        help="""A mafft only option to handle when there are "U"s in aa
            sequences.""")
    fa2tree_parser.add_argument(
        '--min-occupancy', type=float, default=0.3,
        help="""""")
    fa2tree_parser.add_argument(
        '--min-seq-len', type=int, default=10,
        help="""""")
    fa2tree_parser.add_argument(
        '--seed', type=int, default=12345,
        help="""A random number seed. This allows you to reproduce your
            results and helps with debugging the program.""")
    fa2tree_parser.set_defaults(func=fa2tree)


def shrink_step(subparsers):
    """Add shrink step."""
    shrink_parser = subparsers.add_parser(
        'shrink', help=helper("""Trim spurious tips with TreeShrink."""))
    io_args(shrink_parser, '*.tre', '.tt')
    shrink_parser.add_argument(
        '--quantiles', type=float, default=0.05,
        help="""A TreeShrink only option for tree trimming quantiles. The
            default is 0.05.""")
    shrink_parser.set_defaults(func=shrink)


def mask_step(subparsers):
    """Add mask step."""
    mask_parser = subparsers.add_parser(
        'mask', help=helper("""Mask both mono- and (optional) paraphyletic
            tips that belong to the same taxon."""))
    input_files(mask_parser, '*.cln', long='--cleaned-files', short='-c')
    input_files(mask_parser, '*.ts', long='--tree-files', short='-t')
    mask_parser.add_argument(
        '--mask-paraphyletic', action='store_true',
        help="""When masking tree tips, do you want to also mask paraphyletic
            tips.""")
    mask_parser.set_defaults(func=mask)


def tree2fa_step(subparsers):
    """Add tree2fa step."""
    tree2fa_parser = subparsers.add_parser(
        'tree2fa', help=helper(""""""))
    input_files(tree2fa_parser, '*.t', long='--tree-files', short='-t')
    input_files(tree2fa_parser, '*.m', long='--mask-files', short='-m')
    tree2fa_parser.set_defaults(func=tree2fa)


def prune_step(subparsers):
    """Add prune step."""
    prune_parser = subparsers.add_parser(
        'prune', help=helper(""""""))
    prune_parser.set_defaults(func=prune_paralogs)


def orth2fa_step(subparsers):
    """Add orth2fa step."""
    orth2fa_parser = subparsers.add_parser(
        'orth2fa', help=helper(""""""))
    orth2fa_parser.set_defaults(func=orthologs_to_fasta)


def helper(msg):
    """Build a help message."""
    global STEP
    STEP += 1
    return 'STEP {}: {}'.format(STEP, msg)


def io_args(parser, input_filter, output_ext):
    """Add input and output file args to the parser."""
    input_files(parser, input_filter)
    output_args(parser, output_ext)


def input_files(parser, input_filter, long='--input-files', short='-i'):
    """Add input filter arg."""
    global INPUT_ATTRS
    arg_name = long[2:].replace('-', '_')
    INPUT_ATTRS.add(arg_name)

    parser.add_argument(
        short, long, default=input_filter, metavar='FILTER', nargs='+',
        help="""Use this to filter files in an input directory. For example
            'my_project/*filtered*{0}' will select all "{0}" files in the 
            local directory "my_project" with the word "filtered" in them.
            The default is '{0}'.""".format(input_filter))


def output_args(parser, output_ext):
    """Add output file args to the parser."""
    parser.add_argument(
        '-o', '--output-dir', metavar='PATH', default='.',
        help="""Place output files in this directory. The default is the
            current directory.""")

    parser.add_argument(
        '-e', '--output-ext', metavar='EXT', default=output_ext,
        help="""The file extention to use for the output files. 
            The default is '{}'.""".format(output_ext))


def seq_type_arg(parser):
    """Set flag for DNA or amino acid sequences."""
    parser.add_argument(
        '-t', '--seq-type', default='dna', choices=['dna', 'aa'],
        help="""Are we building trees from DNA or amino acid sequences. The
            default is "dna".""")


def cpus_arg(parser):
    """How many CPUS to use."""
    parser.add_argument(
        '--cpus', type=int, default=1,
        help="""Number of CPU processors to use. This is passed to wrapped
            programs that use this option like: mafft or pasta.
            The default will use 1 out of {} CPUs.
            """.format(os.cpu_count()))


def other_args(parser):
    """Placeholder."""
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
        '--min-taxa', type=int, default=2,
        help="""""")

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

    args.out_groups = args.out_groups.strip().strip(r'\'"')
    args.out_groups = [g.strip() for g
                       in re.split(r'\s*[\'",]+\s*', args.out_groups)]
    args.out_groups = [g for g in args.out_groups if g]


if __name__ == "__main__":
    construct()
