"""Wrap raxml functions."""

from os.path import basename, join, splitext
from shutil import move
from . import util
from . import log


def raxml(fasta_file, args):
    """Build a tree with raxml."""
    model = "PROTCATWAG" if args.seq_type == "aa" else "GTRCAT"
    tree, _ = splitext(basename(fasta_file))
    tree = tree + '.tre'
    cmd = ' '.join([
        'raxml',
        '-T {}'.format(args.cpus),
        '-p {}'.format(args.seed),
        '-m {}'.format(model),
        '-s {}'.format(fasta_file),
        '-n {}'.format(tree)])

    with util.cd(args.temp_dir):
        log.subcommand(cmd)

        tree_src = join('RAxML_bestTree.' + tree)
        tree_dst = join(args.output_dir, tree)
        move(tree_src, tree_dst)

    return tree_dst


def raxml_bs(fasta_file, args, replicates=100):
    """Build a bootstrapped tree with raxml."""
    model = "PROTCATWAG" if args.seq_type == "aa" else "GTRCAT"
    tree, _ = splitext(basename(fasta_file))
    tree = tree + '.tre'
    cmd = ' '.join([
        'raxml',
        '-T {}'.format(args.cpus),
        '-f a',
        '-x {}'.format(args.seed),
        '-p {}'.format(args.seed),
        '-m {}'.format(model),
        '-# {}'.format(replicates),
        '-s {}'.format(fasta_file),
        '-n {}'.format(tree)])

    with util.cd(args.temp_dir):
        log.subcommand(cmd)

        tree_src = join('RAxML_bipartitions.' + tree)
        tree_dst = join(args.output_dir, tree)
        move(tree_src, tree_dst)

    return tree_dst
