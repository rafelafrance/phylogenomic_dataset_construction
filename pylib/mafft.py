"""Wrap mafft alignment tool."""

from os.path import basename, join, splitext
import pylib.util as util
import pylib.log as log
import pylib.bio as bio


MAX_ITERATE = 10_000


def mafft(args, fasta_path, temp_dir):
    """Align sequences."""
    in_path = fasta_path
    if args.seq_type == 'aa':
        in_path = bio.adjust_aa_seqs(fasta_path, temp_dir)

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

    aligned = join(args.output_prefix, splitext(basename(fasta_path))[0])
    aligned += '.aligned'

    with util.cd(temp_dir):
        log.subcommand(cmd, out_path=aligned)

    return aligned
