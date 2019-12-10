"""Wrapper for phyx programs."""

# pylint: disable=too-many-arguments

from os.path import basename
from Bio.SeqIO.FastaIO import SimpleFastaParser
from . import util
from . import log
from . import bio


MIN_LEN = 10


def pxclsq(
        fasta_file, output_dir, temp_dir, seq_type, min_occupancy, min_len):
    """Filter aligned sequences for occupancy and length."""
    temp_cleaned = util.file_name(output_dir, fasta_file, '_cleaned.fasta')

    cmd = ' '.join([
        'pxclsq',
        '--aminoacid' if seq_type == 'aa' else '',
        '--prop {}'.format(min_occupancy),
        '--seqf {}'.format(fasta_file),
        '--outf {}'.format(basename(temp_cleaned))])

    with util.cd(temp_dir):
        log.subcommand(cmd)

    cleaned = util.file_name(output_dir, fasta_file, '.cln')

    with open(temp_cleaned) as in_file, open(cleaned, 'w') as out_file:
        for header, seq in SimpleFastaParser(in_file):
            if len(seq.replace('-', '')) >= min_len:
                bio.write_fasta_record(out_file, header, seq)

    return cleaned


def pxrr(tree_file, output_dir, temp_dir):
    """Unroot the tree returned by treeshrink."""
    unrooted = util.file_name(output_dir, tree_file, '.tt')
    cmd = ' '.join([
        'pxrr',
        '--unroot',
        '--treef {}'.format(tree_file),
        '--outf {}'.format(unrooted)])

    with util.cd(temp_dir):
        log.subcommand(cmd)

    return unrooted
