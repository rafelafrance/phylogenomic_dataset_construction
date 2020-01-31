"""Misc. utilities and constants."""

import os
from os.path import basename, expanduser, join, splitext
from glob import glob
from shutil import rmtree
from tempfile import mkdtemp
from contextlib import contextmanager


__VERSION__ = '0.0.1'
__TITLE__ = 'Phylogenomic Dataset Construction'


class Break(Exception):
    """Used to break out of nested loops."""


def shorten(text):
    """Collapse whitespace in a string."""
    return ' '.join(text.split())


@contextmanager
def make_temp_dir(where=None, prefix=None, keep=False):
    """Handle creation and deletion of a temporary directory."""
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


def file_name(base_name, ext=None, dir_=None):
    """Build the output file name."""
    path = splitext(basename(base_name))[0]
    if dir_:
        path = join(dir_, path)
    if ext:
        path += ext
    return path
