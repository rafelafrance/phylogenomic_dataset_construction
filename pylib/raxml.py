"""Wrap raxml functions."""

from os.path import basename, join, splitext
import shutil
import pylib.util as util
import pylib.log as log


def raxml(args, fasta_path, temp_dir):
    """Run a raxml tree build."""
    model = "PROTCATWAG" if args.seq_type == "aa" else "GTRCAT"
    tree = basename(fasta_path)
    tree, _ = splitext(tree)
    tree = tree + '.tre'
    cmd = ' '.join([
        'raxml',
        '-T {}'.format(args.cpus),
        '-p {}'.format(args.seed),
        '-m {}'.format(model),
        '-s {}'.format(fasta_path),
        '-n {}'.format(tree)])
    log.subcommand(cmd, temp_dir)

    tree_src = join('RAxML_bestTree.' + tree)
    tree_dst = join(args.output_prefix, tree)
    shutil.move(tree_src, tree_dst)

    util.remove_all('RAxML_log.*')
    util.remove_all('RAxML_info.*')
    util.remove_all('RAxML_result.*')
    util.remove_all('RAxML_parsimonyTree.*')

    return tree_dst


def raxml_bs(args, fasta_path, temp_dir, replicates=100):
    """Run a boot strapped raxml tree build."""
    model = "PROTCATWAG" if args.seq_type == "aa" else "GTRCAT"
    tree = basename(fasta_path)
    tree, _ = splitext(tree)
    tree = tree + '.tre'
    cmd = ' '.join([
        'raxml',
        '-T {}'.format(args.cpus),
        '-f a',
        '-x {}'.format(args.seed),
        '-p {}'.format(args.seed),
        '-m {}'.format(model),
        '-# {}'.format(replicates),
        '-s {}'.format(fasta_path),
        '-n {}'.format(tree)])
    log.subcommand(cmd, temp_dir)

    tree_src = join('RAxML_bipartitions.' + tree)
    tree_dst = join(args.output_prefix, tree)
    shutil.move(tree_src, tree_dst)

    util.remove_all('RAxML_log.*')
    util.remove_all('RAxML_info.*')
    util.remove_all('RAxML_result.*')
    util.remove_all('RAxML_bestTree.*')
    util.remove_all('RAxML_bootstrap.*')
    util.remove_all('RAxML_parsimonyTree.*')
    util.remove_all('RAxML_bipartitionsBranchLabels.*')

    return tree_dst
