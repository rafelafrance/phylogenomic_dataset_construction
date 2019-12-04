"""Wrap fasttree functions."""

from os.path import basename, join, splitext
import pylib.util as util
import pylib.log as log


def fasttree(fasta_file, args):
    """Build a tree with fasttree."""
    cmd = ['fasttree', '-quiet']
    cmd += ['-wag'] if args.seq_type == 'aa' else ['-nt', '-gtr']
    cmd.append(fasta_file)
    cmd = ' '.join(cmd)

    tree_file = join(args.output_dir, splitext(basename(fasta_file))[0])
    tree_file += '.aligned'

    with util.cd(args.temp_dir):
        log.subcommand(cmd, out_path=tree_file)

    return tree_file
