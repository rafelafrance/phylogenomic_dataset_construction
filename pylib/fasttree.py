"""Wrap fasttree functions."""

from os.path import basename, join, splitext
from . import util
from . import log


def fasttree(fasta_file, output_dir, temp_dir, seq_type):
    """Build a tree with fasttree."""
    cmd = ['fasttree', '-quiet']
    cmd += ['-wag'] if seq_type == 'aa' else ['-nt', '-gtr']
    cmd.append(fasta_file)
    cmd = ' '.join(cmd)

    tree_file = util.file_name(output_dir, fasta_file, '.aligned')

    with util.cd(temp_dir):
        log.subcommand(cmd, out_path=tree_file)

    return tree_file
