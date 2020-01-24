"""Utilities for working with sequences."""

import re
from functools import reduce
from Bio.SeqIO.FastaIO import SimpleFastaParser
from . import util

CODON_LEN = 3

COMPLEMENT = str.maketrans('ACGTUWSMKRYBDHVNXacgtuwsmkrybdhvnx-',
                           'TGCAAWSKMYRVHDBNXtgcaawskmyrvhdbnx-')
AA_REPLACE = str.maketrans('Uu', 'Xx')

SEQ_LEN_CUTOFF = 10_000
SEQ_COUNT_CUTOFF = 1_000
MIN_SEQ = 4


def reverse_complement(seq):
    """Reverse complement a nucleotide sequence. We added some wildcards."""
    return seq.translate(COMPLEMENT)[::-1]


def seqs_too_long(seq_len, seq_type):
    """Check if the sequence is too long for the alignment process."""
    if seq_type == 'aa':
        max_len = SEQ_LEN_CUTOFF
    else:
        max_len = CODON_LEN * SEQ_LEN_CUTOFF
    return seq_len >= max_len


def read_fasta(fasta_file):
    """Read in a fasta file for further processing."""
    with open(fasta_file) as fasta_file:
        return {s[0]: s[1] for s in SimpleFastaParser(fasta_file)}


def read_fasta_records(fasta_file):
    """Read in a fasta file for further processing."""
    with open(fasta_file) as fasta_file:
        for seq_name, seq in SimpleFastaParser(fasta_file):
            yield seq_name, seq


def fasta_record_count(fasta):
    """Count the number of records in a fasta file."""
    with open(fasta) as fasta_file:
        return sum(1 for _ in SimpleFastaParser(fasta_file))


def longest_fasta_seq(fasta):
    """Get the longest sequence length in a fasta file."""
    with open(fasta) as fasta_file:
        return reduce(
            max, [len(s[1]) for s in SimpleFastaParser(fasta_file)], 0)


def adjust_aa_seq(seq):
    """Replace "U"s with "X"s and remove everything after a stop codon."""
    seq = seq.translate(AA_REPLACE)
    return re.sub(r'\*.*', '', seq)


def adjust_aa_seqs(fasta_file, output_dir):
    """Fix up amino acid sequences."""
    file_name = util.file_name(fasta_file, '_aa.fasta')

    with util.cd(output_dir):
        with open(fasta_file) as in_file, open(file_name, 'w') as out_file:
            for header, seq in SimpleFastaParser(in_file):
                seq = adjust_aa_seq(seq)
                write_fasta_record(out_file, header, seq)

    return file_name


def write_fasta_record(out_file, header, seq):
    """Write a fasta record to the file."""
    out_file.write('>')
    out_file.write(header)
    out_file.write('\n')

    out_file.write(seq)
    out_file.write('\n')
