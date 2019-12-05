"""Misc. utilities and constants."""

import sys
import os
from os.path import expanduser, exists, join, split
from glob import glob
from shutil import rmtree
from tempfile import mkdtemp
from contextlib import contextmanager


__VERSION__ = '2.0.0'
__TITLE__ = 'Phylogenomic Dataset Construction'


class StopProcessing(Exception):
    """Stop processing the item in the pipeline."""


def shorten(text):
    """Collapse whitespace in a string."""
    return ' '.join(text.split())


def temp_dir_exists(temp_dir):
    """Make sure the temporary directory exits."""
    if temp_dir and not exists(temp_dir):
        sys.exit('The temporary directory must exist.')


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
