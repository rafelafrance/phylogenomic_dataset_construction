"""Utilities for working with sequences."""

CODON_LEN = 3

COMPLEMENT = str.maketrans('ACGTUWSMKRYBDHVNXacgtuwsmkrybdhvnx-',
                           'TGCAAWSKMYRVHDBNXtgcaawskmyrvhdbnx-')

TOO_LONG = 10_000


def reverse_complement(seq):
    """Reverse complement a nucleotide sequence. We added some wildcards."""
    return seq.translate(COMPLEMENT)[::-1]


def too_long(seq_type, seq_len):
    """Check if the sequence is too long for the alignment process."""
    max_len = TOO_LONG if seq_type == 'aa' else 3 * TOO_LONG
    return seq_len >= max_len


def fasta_record_count(fasta_path):
    """Count the number of records in a fasta file."""
