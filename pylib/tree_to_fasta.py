"""Convert a Newick tree to a fasta file."""
from os.path import basename, join, splitext
from Bio import Phylo
import pylib.bio as bio


def tree_to_fasta(old_fasta, tree_file, args):
    """Convert a Newick tree to a fasta file."""
    tree = Phylo.read(tree_file, 'newick')
    fasta = bio.read_fasta(old_fasta)

    fasta_path = join(args.output_dir, splitext(basename(tree_file))[0])
    fasta_path += 'rr.fa'

    with open(fasta_path, 'w') as out_file:
        for node in tree.get_terminals():
            bio.write_fasta_record(out_file, node.name, fasta[node.name])

    return fasta_path
