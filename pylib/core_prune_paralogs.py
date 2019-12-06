"""Paralogy pruning to infer homologs."""

from os.path import abspath, basename, join, splitext
from glob import glob
import pylib.util as util
import pylib.log as log


def pipeline(args):
    """Build the homology trees."""
    with util.make_temp_dir(
            where=args.temp_dir,
            prefix='{}_'.format(splitext(basename(__file__))[0]),
            keep=args.keep_temp_dir) as temp_dir:
        setup(args, temp_dir)

        for data in get_tree_files(args):
            try:
                pass
            except util.StopProcessing:
                pass


def setup(args, temp_dir):
    """Perform pipeline setup."""
    args.temp_dir = temp_dir  # Put the real temp dir into the args
    log_file = join(args.output_dir, 'build_homology_trees.log')
    log.setup(log_file)


def get_tree_files(args):
    """Get the fasta data to process."""
    log.info('Gathering fasta files')
    pattern = join(args.input_dir, args.input_filter)
    fasta_files = sorted([abspath(p) for p in glob(pattern)])
    if len(fasta_files) == 0:
        log.fatal('No data were found with this mask: "{}".'.format(pattern))
    return [{'fasta': f} for f in fasta_files]


def log_name(path):
    """Update the file path to something more readable for the logs."""
    return '"{}"'.format(splitext(basename(path))[0])
