"""Wrap pasta functions."""

from os.path import abspath, basename, join, splitext
from shutil import move
from . import util
from . import log
from . import bio


def pasta(fasta_file, output_dir, temp_dir, seq_type):
    """Align sequences."""
    in_path = fasta_file
    if seq_type == 'aa':
        in_path = bio.adjust_aa_seqs(fasta_file, temp_dir)

    cmd = ' '.join([
        'run_pasta.py',
        '--datatype {}'.format('Protein' if seq_type == 'aa' else 'DNA'),
        '--input {}'.format(in_path),
        '--output-directory {}'.format(abspath(temp_dir))])

    with util.cd(temp_dir):
        log.subcommand(cmd)

    base_name = splitext(basename(fasta_file))[0]
    temp_aligned = join(temp_dir, 'pastajob.marker001.' + base_name + '.aln')
    aligned = join(output_dir, base_name + '.aln')
    move(temp_aligned, aligned)

    util.remove_files(join(temp_dir, 'pastajob*'))

    return aligned
