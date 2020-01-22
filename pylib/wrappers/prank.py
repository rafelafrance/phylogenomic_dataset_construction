"""Wrap prank alignment tool."""

from pylib import util
from pylib import log
from pylib import bio


def prank(fasta_file, output_dir, temp_dir, seq_type):
    """Align sequences."""
    in_path = fasta_file
    if seq_type == 'aa':
        in_path = bio.adjust_aa_seqs(fasta_file, temp_dir)

    aligned = util.file_name(output_dir, fasta_file, 'ortho.aln')

    cmd = [
        'prank',
        '-d {}'.format(in_path),
        '-o {}'.format(aligned),
        '-protein' if seq_type == 'aa' else '-DNA',
        ]

    cmd = ' '.join(cmd)

    with util.cd(temp_dir):
        subprocess.check_call(cmd, out_path=aligned)

    return aligned
