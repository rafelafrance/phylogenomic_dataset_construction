"""
If no out-group, only output in-group clades with no taxon repeats
If out-group present, extract rooted in-group clades and prune paralogs

Prepare a taxon file, with each line look like (separated by tab):
IN  taxonID1
IN  taxonID2
OUT taxonID3
"""

from os.path import basename, join, splitext
import sys
import oldlib.newick3 as newick3
from . import tree_utils


def prune_rt(tree_file, output_dir, min_taxa, taxon_code_file):
    output_files = []
    in_groups = []
    out_groups = []
    with open(taxon_code_file, "r") as infile:
        for line in infile:
            if len(line) < 3:
                continue
            spls = line.strip().split("\t")
            if spls[0] == "IN":
                in_groups.append(spls[1])
            elif spls[0] == "OUT":
                out_groups.append(spls[1])
            else:
                print("Check taxon_code_file file format")
                sys.exit()
    if len(set(in_groups) & set(out_groups)) > 0:
        print("Taxon ID", set(in_groups) & set(out_groups),
              "in both ingroups and outgroups")
        sys.exit(0)
    print(len(in_groups), "ingroup taxa and", len(out_groups),
          "outgroup taxa read")
    print("Ingroups:", in_groups)
    print("Outgroups:", out_groups)

    with open(tree_file) as infile:
        intree = newick3.parse(infile.readline())
    curroot = intree
    all_names = tree_utils.get_front_names(curroot)
    num_taxa = len(set(all_names))

    # check taxonIDs
    ingroup_names = []
    outgroup_names = []
    for name in all_names:
        if name in in_groups:
            ingroup_names.append(name)
        elif name in out_groups:
            outgroup_names.append(name)
        else:
            print(name, "not in ingroups or outgroups")
            sys.exit()
    if len(set(ingroup_names)) < min_taxa:
        print("not enough ingroup taxa in tree")
        return output_files

    if len(outgroup_names) > 0:  # >= one outgroup, root & cut inclades
        inclades = tree_utils.extract_rooted_ingroup_clades(
            curroot, in_groups, out_groups, min_taxa)
        inclade_count = 0
        for inclade in inclades:
            inclade_count += 1
            output_file = join(
                output_dir, splitext(basename(tree_file))[0])
            output_file += '.inclade{}'.format(inclade_count)
            output_files.append(output_file)
            with open(output_file, "w") as outfile:
                outfile.write(newick3.tostring(inclade) + ";\n")
            orthologs = tree_utils.get_ortho_from_rooted_inclade(inclade)
            ortho_count = 0
            for ortho in orthologs:
                if len(tree_utils.get_front_labels(ortho)) >= min_taxa:
                    ortho_count += 1
                    output_file = join(
                        output_dir, splitext(basename(tree_file))[0])
                    output_file += '.ortho{}.tre'.format(inclade_count)
                    output_files.append(output_file)
                    with open(output_file, "w") as outfile:
                        outfile.write(newick3.tostring(ortho) + ";\n")

    elif len(all_names) == num_taxa:
        # only output ortho tree when there is no taxon repeats
        output_file = join(
            output_dir, splitext(basename(tree_file))[0])
        output_file += '.unrooted-ortho.tre'
        output_files.append(output_file)
        with open(output_file, "w") as outfile:
            outfile.write(newick3.tostring(curroot) + ";\n")

    else:  # do not attempt to infer direction of gene duplication
        # without outgroup info
        print("duplicated taxa in unrooted tree")

    return output_files
