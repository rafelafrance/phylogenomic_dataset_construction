"""Wrap fasttree functions."""

from os.path import basename, join, splitext
import pylib.util as util
import pylib.log as log


def fasttree(args, fasta_path, temp_dir):
    """Build a tree with fasttree."""
    cmd = ['fasttree', '-quiet']
    cmd += ['-wag'] if args.seq_type == 'aa' else ['-nt', '-gtr']
    cmd.append(fasta_path)
    cmd = ' '.join(cmd)

    tree = join(args.output_prefix, splitext(basename(fasta_path))[0])
    tree += '.aligned'

    with util.cd(temp_dir):
        log.subcommand(cmd, out_path=tree)

    return tree
