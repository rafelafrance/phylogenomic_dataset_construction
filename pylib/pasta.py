"""Wrap pasta functions."""

from os.path import abspath, basename, join, splitext
from shutil import copy
import pylib.util as util
import pylib.log as log
import pylib.bio as bio


def pasta(args, fasta_path, temp_dir):
    """Align sequences."""
    in_path = fasta_path
    if args.seq_type == 'aa':
        in_path = bio.adjust_aa_seqs(fasta_path, temp_dir)

    cmd = ' '.join([
        'run_pasta.py',
        '--datatype {}'.format('Protein' if args.seq_type == 'aa' else 'DNA'),
        '--input {}'.format(in_path),
        '--output-directory {}'.format(abspath(temp_dir))])

    with util.cd(temp_dir):
        log.subcommand(cmd)

    base_name = splitext(basename(fasta_path))[0]
    temp_aligned = join(temp_dir, 'pastajob.marker001.' + base_name + '.aln')
    aligned = join(args.output_prefix, base_name + '.aligned')
    copy(temp_aligned, aligned)

    util.remove_files(join(temp_dir, 'pastajob*'))

    return aligned
