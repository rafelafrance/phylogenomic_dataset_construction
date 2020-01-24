"""Wrap raxml-ng functions."""

# pylint: disable=too-many-arguments

import subprocess
from os.path import basename, join, splitext
from shutil import move
from pylib import util


def raxml_ng(
        fasta_file, output_dir, temp_dir,
        seq_type, cpus, seed, output_extension):
    """Build a tree with raxml."""
    model = "Blosum62" if seq_type == "aa" else "GTR"
    tree = util.file_name(output_dir, fasta_file, output_extension)
    cmd = ' '.join([
        'raxml-ng',
        '-T {}'.format(cpus),
        '-p {}'.format(seed),
        '-m {}'.format(model),
        '-s {}'.format(fasta_file),
        '-n {}'.format(tree)])

    with util.cd(temp_dir):
        subprocess.check_call(cmd, shell=True)

        tree_src = join('RAxML_bestTree.' + tree)
        tree_dst = join(output_dir, tree)
        move(tree_src, tree_dst)

    return tree_dst


def raxml_ng_bs(
        fasta_file, output_dir, temp_dir,
        seq_type, cpus, seed, output_extension, replicates=100):
    """Build a bootstrapped tree with raxml."""
    model = "Blosum62" if seq_type == "aa" else "GTR"
    tree = util.file_name(output_dir, fasta_file, output_extension)
    cmd = ' '.join([
        'raxml-ng',
        '-T {}'.format(cpus),
        '-f a',
        '-x {}'.format(seed),
        '-p {}'.format(seed),
        '-m {}'.format(model),
        '-# {}'.format(replicates),
        '-s {}'.format(fasta_file),
        '-n {}'.format(tree)])

    with util.cd(temp_dir):
        subprocess.check_call(cmd, shell=True)

        tree_src = join('RAxML_bipartitions.' + tree)
        tree_dst = join(output_dir, tree)
        move(tree_src, tree_dst)

    return tree_dst
