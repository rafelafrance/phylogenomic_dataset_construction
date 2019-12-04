"""Wrap pasta functions."""

from os.path import abspath, basename, join, splitext
from shutil import move
import pylib.util as util
import pylib.log as log
import pylib.bio as bio


def pasta(fasta_file, args):
    """Align sequences."""
    in_path = fasta_file
    if args.seq_type == 'aa':
        in_path = bio.adjust_aa_seqs(fasta_file, args.temp_dir)

    cmd = ' '.join([
        'run_pasta.py',
        '--datatype {}'.format('Protein' if args.seq_type == 'aa' else 'DNA'),
        '--input {}'.format(in_path),
        '--output-directory {}'.format(abspath(args.temp_dir))])

    with util.cd(args.temp_dir):
        log.subcommand(cmd)

    base_name = splitext(basename(fasta_file))[0]
    temp_aligned = join(
        args.temp_dir, 'pastajob.marker001.' + base_name + '.aln')
    aligned = join(args.output_dir, base_name + '.aln')
    move(temp_aligned, aligned)

    util.remove_files(join(args.temp_dir, 'pastajob*'))

    return aligned
