"""
make sure that raxml, mafft, run_pasta.py, pxclsq and fastree are in the path
"""

import os
import sys
from .mafft_wrapper import mafft
from .pasta_wrapper import pasta
from .pxclsq_wrapper import pxclsq
from .fasttree_wrapper import fasttree
from .raxml_wrapper import raxml
from .raxml_bs_wrapper import raxml_bs
from .seq import read_fasta_file

NUM_SEQ_CUTOFF = 1000  # Use different alignment and tree inference tools


def get_fasta_size(fasta):
    """
    given a fasta file
    output the number of seqs and the length of the longest seq
    """
    longest = 0
    seq_list = read_fasta_file(fasta)
    for seq in seq_list:
        longest = max(longest, len(seq.seq.replace('-', '')))
    return len(seq_list), longest


def fasta_to_tree(path, fasta, num_cores, seq_type):
    """
    given a fasta file
    align, trim alignment and build a tree
    choose appropriate tools depending on size of the fasta file
    """
    if path[-1] != '/':
        path += '/'
    seq_count, _ = get_fasta_size(path + fasta)
    assert seq_count >= 4, 'Less than four sequences in ' + path + fasta
    print((fasta, seq_count, 'sequences'))
    if seq_count >= NUM_SEQ_CUTOFF:  # large cluster
        print('running pasta')
        alignment = pasta(path, fasta, num_cores, seq_type)
        cleaned = pxclsq(path, alignment, 0.1, seq_type)
        if len(read_fasta_file(path + cleaned)) >= 4:
            _ = fasttree(path, cleaned, seq_type)
        else:
            print(('Less than 4 taxa in', cleaned))
    else:  # small cluster
        alignment = mafft(path, fasta, num_cores, seq_type)
        cleaned = pxclsq(path, alignment, 0.1, seq_type)
        if len(read_fasta_file(path + cleaned)) >= 4:
            _ = raxml(path, cleaned, num_cores, seq_type)
        else:
            print(('Less than 4 taxa in', cleaned))


def fasta_to_bs_tree(path, fasta, num_cores, seq_type):
    """
    given a fasta file for the final homolog
    align, trim alignment and build a tree with bootstrap support
    """
    if path[-1] != '/':
        path += '/'
    seq_count, _ = get_fasta_size(path + fasta)
    assert seq_count >= 4, 'Less than four sequences in ' + path + fasta
    print((fasta, seq_count, 'sequences'))
    alignment = mafft(path, fasta, num_cores, seq_type)
    cleaned = pxclsq(path, alignment, 0.2, seq_type)
    if len(read_fasta_file(path + cleaned)) >= 4:
        raxml_bs(path, cleaned, num_cores, seq_type)
    else:
        print(('Less than 4 taxa in', cleaned))


def main(path, num_cores, seq_type, bs, test=False):
    """If test, only process cluster ID that ends with 0."""
    assert seq_type in ('aa', 'dna'), 'Input data type: dna or aa'
    assert bs in ('y', 'n'), 'bootstrap? (y/n)'
    if path[-1] != '/':
        path += '/'

    # check for really long sequences or large alignments.
    # These crashes the alignment program
    for i in os.listdir(path):
        if i.endswith('.fa'):
            seqcount, max_len = get_fasta_size(path + i)
            if (max_len >= 10000 and seq_type == 'aa') or (
                    max_len >= 30000 and seq_type == 'dna'):
                print((i, 'has', seqcount, 'sequences'))
                print(('longest sequence has', max_len, 'characters'))
                print(
                    'Warning: sequence too long. May crash alignment process')

    filecount = 0
    for i in os.listdir(path):
        if not i.endswith('.fa'):
            continue
        if test and (i.split('.')[0])[-1] != '0':
            continue
        filecount += 1
        if bs == 'n':
            fasta_to_tree(path=path, fasta=i, num_cores=num_cores,
                          seq_type=seq_type)
        else:
            fasta_to_bs_tree(path=path, fasta=i, num_cores=num_cores,
                             seq_type=seq_type)
    assert filecount > 0, 'No file end with .fa found in ' + path


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print(
            'python fasta_to_tree.py path number_cores dna/aa bootstrap(y/n)')
        sys.exit(0)

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
