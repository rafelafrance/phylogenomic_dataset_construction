"""Build homology trees."""

import logging
import pylib.util as util
import pylib.bio as bio
# from pylib.mafft import mafft
# from pylib.raxml import raxml, raxml_bs
# from pylib.phyx import pxclsq, pxrr
# from pylib.pasta import pasta
# from pylib.fasttree import fasttree
# from pylib.treeshrink import treeshrink
# from pylib.tree_to_fasta import tree_to_fasta, ortholog_to_fasta
# from oldlib.mask_tips import mask_tips
# from oldlib.prune_paralogs_mi import prune_mi
# from oldlib.prune_paralogs_mo import prune_mo
# from oldlib.prune_paralogs_rt import prune_rt
# from oldlib.prune_orthologs_1to1 import prune_1to1


def check_input(args):
    """Check the input files for good data."""
    for in_file in util.get_input_files(args.input_dir, args.input_filter):
        logging.info('Checking fasta "{}"'.format(in_file))
        too_few_records(in_file)
        seq_too_long(in_file, args.seq_type)


def too_few_records(fasta):
    """Check if the fasta file is too small to make a good tree."""
    if bio.fasta_record_count(fasta) < bio.MIN_SEQ:
        logging.warning('"{}" has fewer than {} records, skipping.'.format(
            fasta, bio.MIN_SEQ))


def seq_too_long(fasta, seq_type):
    """Warn about really long sequences."""
    longest = bio.longest_fasta_seq(fasta)

    if bio.seqs_too_long(longest, seq_type):
        seq_count = bio.fasta_record_count(fasta)
        logging.warning(util.shorten("""{} has {} sequences.
            The longest is {} characters.
            This is too long and may crash the alignment process.
            """.format(fasta, seq_count, longest)))


def fasta_to_tree(args):
    """Build trees from the fasta data."""

    for fasta in util.get_input_files(args.input_dir, args.input_filter):
        if args.bootstrap:
            fasta_to_tree_bs(args, fasta)
        elif bio.fasta_record_count(fasta) >= bio.SEQ_COUNT_CUTOFF:
            pass
        else:
            aligned = mafft(
                fasta, args.output_dir, args.seq_type, args.cpus,
                args.anysymbol)
            cleaned = pxclsq(
                aligned, args.output_dir, args.seq_type, args.min_occupancy,
                args.min_seq_len)
            raxml(
                cleaned, args.output_dir, args.seq_type, args.cpus, args.seed)


def fasta_to_tree_bs(args, fasta):
    """Build trees from the fasta data, bootstrap version."""

    aligned = mafft(
        fasta, args.output_dir, args.seq_type, args.cpus, args.anysymbol)

    cleaned = pxclsq(
        aligned, args.output_dir, args.seq_type, args.min_occupancy,
        args.min_seq_len)

    raxml_bs(cleaned, args.output_dir, args.seq_type, args.cpus, args.seed)


def fasta_to_tree_big(args, fasta):
    """Build trees from the fasta data, large file version."""
    aligned = pasta(fasta, args.output_dir, args.seq_type)

    cleaned = pxclsq(
        aligned, args.output_dir, args.seq_type, args.min_occupancy,
        args.min_seq_len)

    fasttree(cleaned, args.output_dir, args.seq_type)


def fasta_to_tree_normal(args, fasta):
    """Build trees from the fasta data, normal version."""


def tree_shrink(data, args):
    """Remove long branches from trees."""
    logging.info('Shrinking tree for {}'.format(args.log_name))

    trimmed = treeshrink(
        tree, args.output_dir, args.quantiles)
    pxrr(trimmed, args.output_dir)


def mask_tree(data, args):
    """Mask mono- and paraphyletic-tips that belong to the same taxon."""
    logging.info('Masking tree tips for {}'.format(args.log_name))

    mask_tips(cleaned, unrooted, args.output_dir, args.mask_paraphyletic)


def new_fasta(data, args):
    """Mask mono- and paraphyletic-tips that belong to the same taxon."""
    logging.info(
        'Converting masked tree to fasta for {}'.format(args.log_name))

    tree_to_fasta(fasta, masked, args.output_dir)


def prune_paralogs(data, args):
    """Prune paralogs."""
    logging.info('Pruning paralogs for {}'.format(args.log_name))

    if args.prune == 'mi':
        prune_mi(
            masked, args.output_dir, args.min_taxa,
            args.relative_tip_cutoff, args.absolute_tip_cutoff)
    elif args.prune == 'mo':
        prune_mo(masked, args.output_dir, args.min_taxa, args.out_groups)
    elif args.prune == 'rt':
        prune_rt(
            masked, args.output_dir, args.min_taxa,
            args.taxon_code_file)
    else:
        prune_1to1(masked, args.output_dir, args.min_taxa)


def orthologs_to_fasta(data, args):
    """Write ortholog trees to fasta files."""
    ortho_fa = []
    for tree_file in orthos:
        fasta = ortholog_to_fasta(
            fasta2, tree_file, args.output_dir, args.min_taxa)
        if fasta:
            ortho_fa.append_(fasta)
