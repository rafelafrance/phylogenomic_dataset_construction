"""Wrap raxml functions."""

# pylint: disable=too-many-arguments

from os.path import basename, join, splitext
from shutil import move
import subprocess
from pylib import util


def raxml(fasta_file, output_dir, output_ext, seq_type, cpus, seed):
    """Build a tree with raxml."""
    model = "PROTCATWAG" if seq_type == "aa" else "GTRCAT"
    tree = util.file_name(fasta_file, output_ext, output_dir)
    cmd = ' '.join([
        'raxml',
        '-T {}'.format(cpus),
        '-p {}'.format(seed),
        '-m {}'.format(model),
        '-s {}'.format(fasta_file),
        '-n {}'.format(tree)])

    with util.cd(output_dir):
        subprocess.check_call(cmd, shell=True)

        tree_src = join('RAxML_bestTree.' + tree)
        tree_dst = join(output_dir, tree)
        move(tree_src, tree_dst)

    return tree_dst


def raxml_bs(fasta_file, output_dir, output_ext, seq_type, cpus, seed,
             replicates=100):
    """Build a bootstrapped tree with raxml."""
    model = "PROTCATWAG" if seq_type == "aa" else "GTRCAT"
    tree = util.file_name(fasta_file, output_ext)
    cmd = ' '.join([
        'raxml',
        '-T {}'.format(cpus),
        '-f a',
        '-x {}'.format(seed),
        '-p {}'.format(seed),
        '-m {}'.format(model),
        '-# {}'.format(replicates),
        '-s {}'.format(fasta_file),
        '-n {}'.format(tree)])

    with util.cd(output_dir):
        subprocess.check_call(cmd, shell=True)

        tree_src = join('RAxML_bipartitions.' + tree)
        tree_dst = join(output_dir, tree)
        move(tree_src, tree_dst)

    return tree_dst
