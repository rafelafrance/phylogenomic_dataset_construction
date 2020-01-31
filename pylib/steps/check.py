"""Build homology trees."""

import logging
import pylib.util as util
import pylib.bio as bio


def check(args):
    """Check the input files for good data."""
    for fasta in args.input_files:
        logging.info('check input: {}'.format(fasta))
        duplicate_names(fasta)
        too_few_records(fasta)
        seq_too_long(fasta, args.seq_type)


def duplicate_names(fasta):
    """Duplicate names are not allowed."""
    seen = set()
    mentioned = set()
    for seq_name, _ in bio.read_fasta_records(fasta):
        if seq_name in seen and seq_name not in mentioned:
            logging.error(
                'The sequence name {} is in {} more than once'.format(
                    seq_name, fasta))
            mentioned.add(seq_name)
        seen.add(seq_name)


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
