"""Wrap mafft alignment tool."""

import subprocess
from pylib import util
from pylib import bio

MAX_ITERATE = 10_000


def mafft(fasta_file, output_dir, output_ext, seq_type, cpus, anysymbol):
    """Align sequences."""
    in_path = fasta_file
    if seq_type == 'aa':
        in_path = bio.adjust_aa_seqs(fasta_file, output_dir)

    cmd = [
        'mafft',
        '--amino' if seq_type == 'aa' else '--nuc',
        '--thread {}'.format(cpus),
        '--anysymbol' if anysymbol else '']

    if (bio.fasta_record_count(in_path) >= bio.SEQ_COUNT_CUTOFF
            or bio.longest_fasta_seq(in_path) >= bio.SEQ_LEN_CUTOFF):
        cmd.append('--auto')
    else:
        cmd += [
            '--genafpair',
            '--maxiterate {}'.format(MAX_ITERATE),
            '--anysymbol' if anysymbol else '']

    cmd.append(in_path)
    cmd = ' '.join(cmd)

    aligned = util.file_name(fasta_file, output_ext)

    with util.cd(output_dir):
        result = subprocess.check_output(cmd, shell=True)
        with open(aligned, 'wb') as out_file:
            out_file.write(result)

    return aligned
