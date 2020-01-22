"""Wrap pasta functions."""

from os.path import abspath, basename, join, splitext
from shutil import move
import subprocess
from pylib import util
from pylib import bio


def pasta(fasta_file, output_dir, seq_type):
    """Align sequences."""
    in_path = fasta_file
    if seq_type == 'aa':
        in_path = bio.adjust_aa_seqs(fasta_file, output_dir)

    cmd = ' '.join([
        'run_pasta.py',
        '--datatype {}'.format('Protein' if seq_type == 'aa' else 'DNA'),
        '--input {}'.format(in_path),
        '--output-directory {}'.format(abspath(output_dir))])

    with util.cd(output_dir):
        subprocess.check_call(cmd)

    base_name = splitext(basename(fasta_file))[0]
    temp_aligned = join(
        output_dir, 'pastajob.marker001.' + base_name + '.aln')
    aligned = join(output_dir, base_name + '.aln')
    move(temp_aligned, aligned)

    util.remove_files(join(output_dir, 'pastajob*'))

    return aligned
