"""Wrap raxml functions."""

# pylint: disable=too-many-arguments

from shutil import move
import subprocess
from pylib import util


def raxml(fasta_file, output_dir, output_ext, seq_type, cpus, seed):
    """Build a tree with raxml."""
    model = "PROTCATWAG" if seq_type == "aa" else "GTRCAT"
    tree = util.file_name(fasta_file, output_ext)
    cmd = ' '.join([
        'raxml',
        '-T {}'.format(cpus),
        '-p {}'.format(seed),
        '-m {}'.format(model),
        '-s {}'.format(fasta_file),
        '-n {}'.format(tree)])

    with util.cd(output_dir):
        subprocess.check_call(cmd, shell=True)
        tree_src = 'RAxML_bestTree.' + tree
        move(tree_src, tree)
        util.remove_files('RAxML_*')

    return tree


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
        tree_src = 'RAxML_bipartitions.' + tree
        move(tree_src, tree)
        util.remove_files('RAxML_*')

    return tree
