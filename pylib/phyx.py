"""Wrapper for phyx programs."""

from os.path import basename, join, splitext
from Bio.SeqIO.FastaIO import SimpleFastaParser
from . import util
from . import log
from . import bio


COLUMN_OCCUPANCY_SM = 0.1
COLUMN_OCCUPANCY_LG = 0.2
MIN_LEN = 10


def pxclsq(fasta_file, args, min_occupancy, min_len=MIN_LEN):
    """Filter aligned sequences for occupancy and length."""
    temp_cleaned = join(args.temp_dir, splitext(basename(fasta_file))[0])
    temp_cleaned += '_cleaned.fasta'

    cmd = ' '.join([
        'pxclsq',
        '--aminoacid' if args.seq_type == 'aa' else '',
        '--prop {}'.format(min_occupancy),
        '--seqf {}'.format(fasta_file),
        '--outf {}'.format(basename(temp_cleaned))])

    with util.cd(args.temp_dir):
        log.subcommand(cmd)

    cleaned = join(args.output_dir, splitext(basename(fasta_file))[0])
    cleaned += '.cln'

    with open(temp_cleaned) as in_file, open(cleaned, 'w') as out_file:
        for header, seq in SimpleFastaParser(in_file):
            if len(seq.replace('-', '')) >= min_len:
                bio.write_fasta_record(out_file, header, seq)

    return cleaned


def pxrr(tree_file, args):
    """Unroot the tree returned by treeshrink."""
    unrooted = join(args.output_dir, splitext(basename(tree_file))[0] + '.tt')
    cmd = ' '.join([
        'pxrr',
        '--unroot',
        '--treef {}'.format(tree_file),
        '--outf {}'.format(unrooted)])

    with util.cd(args.temp_dir):
        log.subcommand(cmd)

    return unrooted
