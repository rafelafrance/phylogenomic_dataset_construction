"""Wrap mafft alignment tool."""

import re
from os.path import basename, join, splitext
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pylib.util as util
import pylib.log as log
import pylib.bio as bio


MAX_ITERATE = 10_000
AA_REPLACE = str.maketrans('Uu', 'Xx')


def mafft(args, fasta_path, temp_dir):
    """Align sequences."""
    cmd = [
        'mafft',
        '--thread {}'.format(args.cpus)]

    if args.seq_type == 'dna':
        cmd.append('--nuc')
        in_path = fasta_path
    else:
        cmd.append('--amino')
        in_path = adjust_aa_seqs(fasta_path, temp_dir)

    if bio.fasta_record_count(in_path) >= bio.SEQ_COUNT_CUTOFF \
            or bio.longest_fasta_seq(in_path) >= bio.SEQ_LEN_CUTOFF:
        cmd.append('--auto')
    else:
        cmd += [
            '--genafpair',
            '--maxiterate {}'.format(MAX_ITERATE)]

    cmd.append(in_path)

    cmd = ' '.join(cmd)

    aligned = join(args.output_prefix, splitext(basename(fasta_path))[0])
    aligned += '.aligned'

    with util.cd(temp_dir):
        log.capture(cmd, aligned, '.')

    return aligned


def adjust_aa_seqs(fasta_path, temp_dir):
    """Fix up amino acid sequences."""
    out_path = join(temp_dir, splitext(basename(fasta_path))[0]) + '_aa.fasta'

    with open(fasta_path) as in_file, open(out_path, 'w') as out_file:
        for header, seq in SimpleFastaParser(in_file):
            seq = seq.translate(AA_REPLACE)
            seq = re.sub(r'\*.*', '', seq)

            bio.write_fasta_record(out_file, header, seq)

    return out_path
