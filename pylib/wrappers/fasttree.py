"""Wrap fasttree functions."""

import subprocess
from pylib import util


def fasttree(fasta_file, output_dir, output_ext, seq_type):
    """Build a tree with fasttree."""
    cmd = ['fasttree', '-quiet']
    cmd += ['-wag'] if seq_type == 'aa' else ['-nt', '-gtr']
    cmd.append(fasta_file)
    cmd = ' '.join(cmd)

    tree_file = util.file_name(fasta_file, output_ext)

    with util.cd(output_dir):
        subprocess.check_call(cmd, shell=True)

    return tree_file
