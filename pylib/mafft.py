"""Wrap mafft alignment tool."""

from os.path import basename, join, splitext
from . import util
from . import log
from . import bio


MAX_ITERATE = 10_000


def mafft(fasta_file, args):
    """Align sequences."""
    in_path = fasta_file
    if args.seq_type == 'aa':
        in_path = bio.adjust_aa_seqs(fasta_file, args.temp_dir)

    cmd = [
        'mafft',
        '--amino' if args.seq_type == 'aa' else '--nuc',
        '--thread {}'.format(args.cpus),
        '--anysymbol' if args.anysymbol else '']

    if bio.fasta_record_count(in_path) >= bio.SEQ_COUNT_CUTOFF \
            or bio.longest_fasta_seq(in_path) >= bio.SEQ_LEN_CUTOFF:
        cmd.append('--auto')
    else:
        cmd += [
            '--genafpair',
            '--maxiterate {}'.format(MAX_ITERATE),
            '--anysymbol' if args.anysymbol else '']

    cmd.append(in_path)
    cmd = ' '.join(cmd)

    aligned = join(args.output_dir, splitext(basename(fasta_file))[0])
    aligned += '.aln'

    with util.cd(args.temp_dir):
        log.subcommand(cmd, out_path=aligned)

    return aligned
