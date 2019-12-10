"""Build homology trees."""

from os.path import abspath, basename, join, splitext
from glob import glob
import pylib.util as util
import pylib.bio as bio
import pylib.log as log
from pylib.mafft import mafft
from pylib.raxml import raxml, raxml_bs
from pylib.phyx import pxclsq, pxrr
from pylib.pasta import pasta
from pylib.fasttree import fasttree
from pylib.treeshrink import treeshrink
from pylib.tree_to_fasta import tree_to_fasta, ortholog_to_fasta
from oldlib.mask_tips import mask_tips
from oldlib.prune_paralogs_mi import prune_mi
from oldlib.prune_paralogs_mo import prune_mo
from oldlib.prune_paralogs_rt import prune_rt
from oldlib.prune_orthologs_1to1 import prune_1to1


def pipeline(args):
    """Build the homology trees."""
    with util.make_temp_dir(
            where=args.temp_dir,
            prefix='{}_'.format(splitext(basename(__file__))[0]),
            keep=args.keep_temp_dir) as temp_dir:
        setup(args, temp_dir)

        for data in get_fasta_files(args):
            setattr(
                args, 'log_name',
                '"{}"'.format(splitext(basename(data['fasta']))[0]))
            try:
                too_few_records(data, args)
                seq_too_long(data, args)
                fasta_to_tree(data, args)
                tree_shrink(data, args)
                mask_tree(data, args)
                new_fasta(data, args)
                prune_paralogs(data, args)
                orthologs_to_fasta(data, args)

            except util.StopProcessing:
                pass


def setup(args, temp_dir):
    """Perform pipeline setup."""
    args.temp_dir = temp_dir  # Put the real temp dir into the args
    log_file = join(args.output_dir, 'build_homology_trees.log')
    log.setup(log_file)
    setattr(args, 'groups_processed', False)


def get_fasta_files(args):
    """Get the fasta data to process."""
    log.info('Gathering fasta files')
    pattern = join(args.input_dir, args.input_filter)
    fasta_files = sorted([abspath(p) for p in glob(pattern)])
    if len(fasta_files) == 0:
        log.fatal('No data were found with this mask: "{}".'.format(pattern))
    return [{'fasta': f} for f in fasta_files]


def too_few_records(data, args):
    """Check if the fasta file is too small to make a good tree."""
    log.info('Checking fasta for {}'.format(args.log_name))

    if bio.fasta_record_count(data['fasta']) < bio.MIN_SEQ:
        log.warn('"{}" has fewer than {} records, skipping.'.format(
            data['fasta'], bio.MIN_SEQ))
        raise util.StopProcessing()


def seq_too_long(data, args):
    """Warn about really long sequences."""
    longest = bio.longest_fasta_seq(data['fasta'])

    if bio.seqs_too_long(longest, args.seq_type):
        seq_count = bio.fasta_record_count(data['fasta'])
        log.warn(util.shorten("""{} has {} sequences.
            The longest is {} characters.
            This is too long and may crash the alignment process.
            """.format(data['fasta'], seq_count, longest)))


def fasta_to_tree(data, args):
    """Build trees from the fasta data."""
    log.info('Converting fasta to tree for {}'.format(args.log_name))

    if args.bootstrap:
        data['aligned'] = mafft(
            data['fasta'], args.output_dir, args.temp_dir,
            args.seq_type, args.cpus, args.anysymbol)
        data['cleaned'] = pxclsq(
            data['aligned'], args.output_dir, args.temp_dir,
            args.seq_type, args.min_occupancy, args.min_seq_len)
        data['tree'] = raxml_bs(
            data['cleaned'], args.output_dir, args.temp_dir,
            args.seq_type, args.cpus, args.seed)

    elif bio.fasta_record_count(data['fasta']) >= bio.SEQ_COUNT_CUTOFF:
        data['aligned'] = pasta(
            data['fasta'], args.output_dir, args.temp_dir, args.seq_type)
        data['cleaned'] = pxclsq(
            data['aligned'], args.output_dir, args.temp_dir,
            args.seq_type, args.min_occupancy, args.min_seq_len)
        data['tree'] = fasttree(
            data['cleaned'], args.output_dir, args.temp_dir, args.seq_type)

    else:
        data['aligned'] = mafft(
            data['fasta'], args.output_dir, args.temp_dir,
            args.seq_type, args.cpus, args.anysymbol)
        data['cleaned'] = pxclsq(
            data['aligned'], args.output_dir, args.temp_dir,
            args.seq_type, args.min_occupancy, args.min_seq_len)
        data['tree'] = raxml(
            data['cleaned'], args.output_dir, args.temp_dir,
            args.seq_type, args.cpus, args.seed)


def tree_shrink(data, args):
    """Remove long branches from trees."""
    log.info('Shrinking tree for {}'.format(args.log_name))

    data['trimmed'] = treeshrink(
        data['tree'], args.output_dir, args.temp_dir, args.quantiles)
    data['unrooted'] = pxrr(data['trimmed'], args.output_dir, args.temp_dir)


def mask_tree(data, args):
    """Mask mono- and paraphyletic-tips that belong to the same taxon."""
    log.info('Masking tree tips for {}'.format(args.log_name))

    data['masked'] = mask_tips(
        data['cleaned'], data['unrooted'], args.output_dir,
        args.mask_paraphyletic)


def new_fasta(data, args):
    """Mask mono- and paraphyletic-tips that belong to the same taxon."""
    log.info('Converting masked tree to fasta for {}'.format(args.log_name))

    data['fasta2'] = tree_to_fasta(
        data['fasta'], data['masked'], args.output_dir)


def prune_paralogs(data, args):
    """Prune paralogs."""
    log.info('Pruning paralogs for {}'.format(args.log_name))

    if args.prune == 'mi':
        data['orthos'] = prune_mi(
            data['masked'], args.output_dir, args.min_taxa,
            args.relative_tip_cutoff, args.absolute_tip_cutoff)
    elif args.prune == 'mo':
        data['orthos'] = prune_mo(
            data['masked'], args.output_dir, args.min_taxa, args.out_groups)
    elif args.prune == 'rt':
        data['orthos'] = prune_rt(
            data['masked'], args.output_dir, args.min_taxa,
            args.taxon_code_file)
    else:
        data['orthos'] = prune_1to1(
            data['masked'], args.output_dir, args.min_taxa)


def orthologs_to_fasta(data, args):
    """Write ortholog trees to fasta files."""
    data['ortho_fa'] = []
    for tree_file in data['orthos']:
        fasta = ortholog_to_fasta(
            data['fasta2'], tree_file, args.output_dir, args.min_taxa)
        if fasta:
            data['ortho_fa'].append_(fasta)
