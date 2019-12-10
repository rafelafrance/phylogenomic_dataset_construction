"""Misc. utilities and constants."""

import os
from os.path import expanduser
from glob import glob
from shutil import rmtree
from tempfile import mkdtemp
from contextlib import contextmanager


__VERSION__ = '0.0.1'
__TITLE__ = 'Phylogenomic Dataset Construction'


class StopProcessing(Exception):
    """Stop processing the item in the pipeline."""


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
