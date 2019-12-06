import os
from . import newick3
from .tree_utils import *


def get_clusterID(filename):
    return filename.split(".")[0]


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "python prune_paralogs_from_rooted_trees.py homoTreeDIR "
            "tree_file_ending minimal_taxa outDIR")
        sys.exit(0)

    inDIR = sys.argv[1] + "/"
    tree_file_ending = sys.argv[2]
    MIN_TAXA = int(sys.argv[3])
    outDIR = sys.argv[4] + "/"
    for i in os.listdir(inDIR):
        if not i.endswith(tree_file_ending):
            continue
        print(i)
        outID = outDIR + get_clusterID(i)
        with open(inDIR + i, "r") as infile:
            intree = newick3.parse(infile.readline())
        orthologs = get_ortho_from_rooted_inclade(intree)
        count = 1
        for ortho in orthologs:
            if len(set(get_front_names(ortho))) >= MIN_TAXA:
                with open(outID + ".ortho" + str(count) + ".tre",
                          "w") as outfile:
                    outstring = newick3.tostring(ortho)
                    outfile.write(outstring + ";\n")
                count += 1
