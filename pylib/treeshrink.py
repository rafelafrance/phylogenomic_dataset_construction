"""Wrapper for the treeshrink program."""

from os.path import basename, join, splitext
import pylib.util as util
import pylib.log as log


def treeshrink(args, tree, temp_dir):
    """Remove long branches from a tree."""
    subdir = join(temp_dir, splitext(basename(tree))[0])

    cmd = ' '.join([
        'run_treeshrink.py',
        '--tree {}'.format(tree),
        '--centroid',
        '--mode per-gene',
        '--quantiles {}'.format(args.quantiles),
        '--outdir {}'.format(basename(subdir)),
        '--tempdir {}'.format(basename(subdir))])

    with util.cd(temp_dir):
        log.subcommand(cmd)

    tree_src = join(subdir, tree)
    tree_dst = join(args.output_prefix, splitext(tree)[0] + '.ts')

    with open(tree_src) as in_file, open(tree_dst, 'w') as out_file:
        content = in_file.read()
        out_file.write(content.replace("'", ''))

    return tree_dst
