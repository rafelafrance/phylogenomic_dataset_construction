"""Build homology trees."""

import logging
import pylib.util as util
import pylib.bio as bio


def check(args):
    """Check the input files for good data."""
    for fasta in util.get_input_files(args.input_dir, args.input_filter):
        logging.info('check input: {}'.format(fasta))
        too_few_records(fasta)
        seq_too_long(fasta, args.seq_type)


def too_few_records(fasta):
    """Check if the fasta file is too small to make a good tree."""
    if bio.fasta_record_count(fasta) < bio.MIN_SEQ:
        logging.warning('{} has fewer than {} records, skipping.'.format(
            fasta, bio.MIN_SEQ))


def seq_too_long(fasta, seq_type):
    """Warn about really long sequences."""
    longest = bio.longest_fasta_seq(fasta)

    if bio.seqs_too_long(longest, seq_type):
        seq_count = bio.fasta_record_count(fasta)
        logging.warning(util.shorten("""{} has {} sequences.
            The longest is {} characters.
            This is too long and may crash the alignment process.
            """.format(fasta, seq_count, longest)))
