"""Wrapper for the treeshrink program."""

from glob import glob
from shutil import rmtree
import subprocess
from pylib import util

EXT_IN = '.tre'
EXT_OUT = '.ts'


def treeshrink(tree_file, output_dir, output_ext, quantiles):
    """Remove long branches from a tree."""
    subdir = util.file_name(tree_file)

    cmd = ' '.join([
        'run_treeshrink.py',
        '--tree {}'.format(tree_file),
        '--centroid',
        '--mode per-gene',
        '--quantiles {}'.format(quantiles),
        '--outdir {}'.format(subdir),
        '--tempdir {}'.format(subdir)])

    with util.cd(output_dir):
        subprocess.check_call(cmd, shell=True)

        mask = util.file_name(subdir + '_*', ext=EXT_IN, dir_=subdir)
        tree_src = glob(mask)[0]
        tree_dst = util.file_name(tree_file, output_ext + EXT_OUT)

        with open(tree_src) as in_file, open(tree_dst, 'w') as out_file:
            content = in_file.read()
            out_file.write(content.replace("'", ''))

        rmtree(subdir)

    return tree_dst
