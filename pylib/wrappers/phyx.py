"""Wrapper for phyx programs."""

# pylint: disable=too-many-arguments

from os.path import basename
import subprocess
from Bio.SeqIO.FastaIO import SimpleFastaParser
from pylib import util
from pylib import bio


MIN_LEN = 10
EXT_PXCLSQ = '.cln'


def pxclsq(fasta_file, output_dir, output_ext, seq_type, min_occupancy,
           min_len):
    """Filter aligned sequences for occupancy and length."""
    ext = output_ext + EXT_PXCLSQ
    temp_cleaned = util.file_name(fasta_file, ext)

    cmd = ' '.join([
        'pxclsq',
        '--aminoacid' if seq_type == 'aa' else '',
        '--prop {}'.format(min_occupancy),
        '--seqf {}'.format(fasta_file),
        '--outf {}'.format(basename(temp_cleaned))])

    cleaned = util.file_name(fasta_file, output_ext)

    with util.cd(output_dir):
        subprocess.check_call(cmd, shell=True)
        with open(temp_cleaned) as in_file, open(cleaned, 'w') as out_file:
            for header, seq in SimpleFastaParser(in_file):
                if len(seq.replace('-', '')) >= min_len:
                    bio.write_fasta_record(out_file, header, seq)

        util.remove_files('phyx.logfile')

    return cleaned


def pxrr(tree_file, output_dir):
    """Unroot the tree returned by treeshrink."""
    unrooted = util.file_name(tree_file)
    cmd = ' '.join([
        'pxrr',
        '--unroot',
        '--treef {}'.format(tree_file),
        '--outf {}'.format(unrooted)])

    with util.cd(output_dir):
        subprocess.check_call(cmd, shell=True)
        util.remove_files('phyx.logfile')

    return unrooted
