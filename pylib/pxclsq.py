"""Wrapper for pxclsq from phyx."""

from os.path import basename, join, splitext
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pylib.util as util
import pylib.log as log
import pylib.bio as bio


COLUMN_OCCUPANCY_SM = 0.1
COLUMN_OCCUPANCY_LG = 0.2
MIN_LEN = 10


def pxclsq(args, fasta_path, temp_dir, min_occupancy, min_len=MIN_LEN):
    """Filter aligned sequences for occupancy and length."""
    temp_cleaned = join(temp_dir, splitext(basename(fasta_path))[0])
    temp_cleaned += '_cleaned.fasta'

    cmd = ' '.join([
        'pxclsq',
        '--aminoacid' if args.seq_type == 'aa' else '',
        '--prop {}'.format(min_occupancy),
        '--seqf {}'.format(fasta_path),
        '--outf {}'.format(basename(temp_cleaned))])

    with util.cd(temp_dir):
        log.subcommand(cmd)

    cleaned = join(args.output_prefix, splitext(basename(fasta_path))[0])
    cleaned += '_cleaned.fasta'

    with open(temp_cleaned) as in_file, open(cleaned, 'w') as out_file:
        for header, seq in SimpleFastaParser(in_file):
            if len(seq.replace('-', '')) >= min_len:
                bio.write_fasta_record(out_file, header, seq)

    return cleaned
