"""Convert a Newick tree to a fasta file."""

from Bio import Phylo
from pylib import bio
from pylib import util


def tree_to_fasta(old_fasta, tree_file, output_dir, output_extension):
    """Convert a Newick tree to a fasta file."""
    tree = Phylo.read(tree_file, 'newick')
    fasta = bio.read_fasta(old_fasta)

    fasta_path = util.file_name(output_dir, tree_file, output_extension)

    with open(fasta_path, 'w') as out_file:
        for node in tree.get_terminals():
            bio.write_fasta_record(out_file, node.name, fasta[node.name])

    return fasta_path


def ortholog_to_fasta(
        old_fasta, tree_file, output_dir, min_taxa, output_extension):
    """Convert a Newick tree to a fasta file using extra checks."""
    tree = Phylo.read(tree_file, 'newick')
    fasta = bio.read_fasta(old_fasta)

    fasta_path = util.file_name(output_dir, tree_file, output_extension)

    taxa = set(n.name.split('@')[0]
               for n in tree.get_terminals() if '@' in n.name)
    if len(taxa) < min_taxa:
        return None

    with open(fasta_path, 'w') as out_file:
        for node in tree.get_terminals():
            bio.write_fasta_record(out_file, node.name, fasta[node.name])

    return fasta_path
