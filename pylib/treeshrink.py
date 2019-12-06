"""Wrapper for the treeshrink program."""

from os.path import basename, join, splitext
from . import util
from . import log


def treeshrink(tree_file, output_dir, temp_dir, quantiles):
    """Remove long branches from a tree."""
    subdir = join(temp_dir, splitext(basename(tree_file))[0])

    cmd = ' '.join([
        'run_treeshrink.py',
        '--tree {}'.format(tree_file),
        '--centroid',
        '--mode per-gene',
        '--quantiles {}'.format(quantiles),
        '--outdir {}'.format(basename(subdir)),
        '--tempdir {}'.format(basename(subdir))])

    with util.cd(temp_dir):
        log.subcommand(cmd)

    tree_src = join(subdir, tree_file)
    tree_dst = join(output_dir, splitext(tree_file)[0] + '.ts')

    with open(tree_src) as in_file, open(tree_dst, 'w') as out_file:
        content = in_file.read()
        out_file.write(content.replace("'", ''))

    return tree_dst
