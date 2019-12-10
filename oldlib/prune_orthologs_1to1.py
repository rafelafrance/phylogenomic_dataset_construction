from os.path import basename, join, splitext
from shutil import copyfile
from . import newick3
from .tree_utils import *


def prune_1to1(tree_file, output_dir, min_taxa, min_bootstrap=0.0):
    output_files = []
    with open(tree_file) as infile:
        intree = newick3.parse(infile.readline())
    names = get_front_names(intree)
    num_tips, num_taxa = len(names), len(set(names))
    print("number of tips:", num_tips, "number of taxa:", num_taxa)
    if num_tips == num_taxa and num_taxa >= min_taxa:
        if min_bootstrap > 0.0 and not pass_boot_filter(intree, min_bootstrap):
            return output_files
        output_file = join(output_dir, splitext(basename(tree_file))[0])
        output_file += '_1to1ortho.tre'
        copyfile(tree_file, output_file)
        output_files.append(output_file)
    return output_files
