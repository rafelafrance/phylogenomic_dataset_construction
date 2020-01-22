"""Misc. utilities and constants."""

import sys
import os
from os.path import basename, expanduser, join, splitext
import logging
from glob import glob
from shutil import rmtree
from tempfile import mkdtemp
from contextlib import contextmanager


__VERSION__ = '0.0.1'
__TITLE__ = 'Phylogenomic Dataset Construction'


def shorten(text):
    """Collapse whitespace in a string."""
    return ' '.join(text.split())


@contextmanager
def make_temp_dir(where=None, prefix=None, keep=False):
    """Handle creation and deletion of temporary directory."""
    temp_dir = mkdtemp(prefix=prefix, dir=where)
    try:
        yield temp_dir
    finally:
        if not keep or not where:
            rmtree(temp_dir)


@contextmanager
def cd(new_dir):
    """Change to a new working directory."""
    prev_dir = os.getcwd()
    os.chdir(expanduser(new_dir))
    try:
        yield
    finally:
        os.chdir(prev_dir)


def remove_files(pattern):
    """Remove all files matching the given pattern."""
    for path in glob(pattern):
        os.remove(path)


def taxon_id(header):
    """Split the name and return the taxon ID part."""
    return header.split('@')[0]


def file_name(output_dir, path, suffix=None):
    """Build the output file name."""
    path = join(output_dir, splitext(basename(path))[0])
    if suffix:
        path += suffix
    return path


def get_input_files(args):
    """Get a list of the input files."""
    pattern = join(args.input_dir, args.input_filter)
    in_files = sorted([p for p in glob(pattern)])
    if len(in_files) == 0:
        logging.critical(
            'No files were found with this mask: "{}".'.format(pattern))
        sys.exit(1)
    return in_files
