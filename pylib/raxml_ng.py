"""Wrap raxml-ng functions."""

# pylint: disable=too-many-arguments

from os.path import basename, join, splitext
from shutil import move
from . import util
from . import log


def raxml_ng(
        fasta_file, output_dir, temp_dir,
        seq_type, cpus, seed):
    """Build a tree with raxml."""
    model = "Blosum62" if seq_type == "aa" else "GTR"
    tree, _ = splitext(basename(fasta_file))
    tree = tree + '.tre'
    cmd = ' '.join([
        'raxml-ng',
        '-T {}'.format(cpus),
        '-p {}'.format(seed),
        '-m {}'.format(model),
        '-s {}'.format(fasta_file),
        '-n {}'.format(tree)])

    with util.cd(temp_dir):
        log.subcommand(cmd)

        tree_src = join('RAxML_bestTree.' + tree)
        tree_dst = join(output_dir, tree)
        move(tree_src, tree_dst)

    return tree_dst


def raxml_ng_bs(
        fasta_file, output_dir, temp_dir,
        seq_type, cpus, seed, replicates=100):
    """Build a bootstrapped tree with raxml."""
    model = "Blosum62" if seq_type == "aa" else "GTR"
    tree, _ = splitext(basename(fasta_file))
    tree = tree + '.tre'
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
        log.subcommand(cmd)

        tree_src = join('RAxML_bipartitions.' + tree)
        tree_dst = join(output_dir, tree)
        move(tree_src, tree_dst)

    return tree_dst
