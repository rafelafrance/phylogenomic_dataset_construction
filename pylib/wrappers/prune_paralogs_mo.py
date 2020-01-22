"""
Taxon duplication? --No--> output one-to-one orthologs
        |
       Yes
        |
Out-group present? --No--> ignore this homolog
        |
       Yes
        |
Out-group taxon duplication? --Yes--> ignore this homolog
        |
        No
        |
Out-group monophyletic? --No--> ignore this homolog
        |
       Yes
        |
Infer orthologs by using monophyletic, non-repeating out-groups

If not to output 1-to-1 orthologs, for example, already analysed these
set OUTPUT_1TO1_ORTHOLOGS to False
"""

from shutil import copyfile
from pylib import util, phylo3, newick3

OUTPUT_1TO1_ORTHOLOGS = True


# if pattern changes, change it here
# given tip label, return taxon name identifier
def get_name(label):
    return label.split("@")[0]


def get_front_labels(node):
    leaves = node.leaves()
    return [i.label for i in leaves]


def get_back_labels(node, root):
    all_labels = get_front_labels(root)
    front_labels = get_front_labels(node)
    return set(all_labels) - set(front_labels)  # labels do not repeat


def get_front_names(node):  # may include duplicates
    labels = get_front_labels(node)
    return [get_name(i) for i in labels]


def get_front_outgroup_names(node, out_groups):
    names = get_front_names(node)
    return [i for i in names if i in out_groups]


def get_back_names(node, root):  # may include duplicates
    back_labels = get_back_labels(node, root)
    return [get_name(i) for i in back_labels]


def remove_kink(node, cur_root):
    if node == cur_root and cur_root.nchildren == 2:
        # move the root away to an adjacent none-tip
        if cur_root.children[0].istip:  # the other child is not tip
            cur_root = phylo3.reroot(cur_root, cur_root.children[1])
        else:
            cur_root = phylo3.reroot(cur_root, cur_root.children[0])
    # ---node---< all nodes should have one child only now
    length = node.length + (node.children[0]).length
    par = node.parent
    kink = node
    node = node.children[0]
    # parent--kink---node<
    par.remove_child(kink)
    par.add_child(node)
    node.length = length
    return node, cur_root


# check if out-groups are monophyletic and non-repeating and reroot
# otherwise return None
def reroot_with_monophyletic_outgroups(root, out_groups):
    lvs = root.leaves()
    outgroup_matches = {}  # key is label, value is the tip node object
    # Since no taxon repeat in out-groups name and leaf is one-to-one
    outgroup_labels = []
    for leaf in lvs:
        label = leaf.label
        name = get_name(label)
        if name in out_groups:
            outgroup_matches[label] = leaf
            outgroup_labels.append(leaf.label)
    if len(outgroup_labels) == 1:  # one single out-group
        # cannot reroot on a tip so have to go one more node into the in-group
        new_root = outgroup_matches[outgroup_labels[0]].parent
        return phylo3.reroot(root, new_root)

    newroot = None
    for node in root.iternodes():
        if node == root:
            continue  # skip the root
        front_names = get_front_names(node)
        back_names = get_back_names(node, root)
        front_in_names, front_out_names, back_in_names, back_out_names =\
            0, 0, 0, 0
        for i in front_names:
            if i in out_groups:
                front_out_names += 1
            else:
                front_in_names += 1
        for j in back_names:
            if j in out_groups:
                back_out_names += 1
            else:
                back_in_names += 1
        if front_in_names == 0 and front_out_names > 0 \
                and back_in_names > 0 and back_out_names == 0:
            newroot = node  # in-group at back, out-group in front
            break
        if front_in_names > 0 and front_out_names == 0 \
                and back_in_names == 0 and back_out_names > 0:
            newroot = node.parent  # in-group in front, out-group at back
            break
    if newroot is not None:
        return phylo3.reroot(root, newroot)
    return None


def prune_paralogs_from_rerooted_homotree(root, out_groups):
    if len(get_front_names(root)) == len(set(get_front_names(root))):
        return root  # no pruning needed
    # check for duplications at the root first
    # one or two of the trifurcating root clades are in-group clades
    node0, node1, node2 = root.children[0], root.children[1], root.children[2]
    out0 = len(get_front_outgroup_names(node0, out_groups))
    out1 = len(get_front_outgroup_names(node1, out_groups))
    out2 = len(get_front_outgroup_names(node2, out_groups))
    if out0 == 0 and out1 == 0:  # 0 and 1 are the in-group clades
        name_set0 = set(get_front_names(node0))
        name_set1 = set(get_front_names(node1))
        if len(name_set0.intersection(name_set1)) > 0:
            if len(name_set0) > len(name_set1):  # cut the side with less taxa
                root.remove_child(node1)
                node1.prune()
            else:
                root.remove_child(node0)
                node0.prune()
    elif out1 == 0 and out2 == 0:  # 1 and 2 are the in-group clades
        name_set1 = set(get_front_names(node1))
        name_set2 = set(get_front_names(node2))
        if len(name_set1.intersection(name_set2)) > 0:
            if len(name_set1) > len(name_set2):  # cut the side with less taxa
                root.remove_child(node2)
                node2.prune()
            else:
                root.remove_child(node1)
                node1.prune()
    elif out0 == 0 and out2 == 0:  # 0 and 2 are the in-group clades
        name_set0 = set(get_front_names(node0))
        name_set2 = set(get_front_names(node2))
        if len(name_set0.intersection(name_set2)) > 0:
            if len(name_set0) > len(name_set2):  # cut the side with less taxa
                root.remove_child(node2)
                node2.prune()
            else:
                root.remove_child(node0)
                node0.prune()
    while len(get_front_names(root)) > len(set(get_front_names(root))):
        for node in root.iternodes(order=0):  # PREORDER, root to tip
            if node.istip or node == root:
                continue
            child0, child1 = node.children[0], node.children[1]
            name_set0 = set(get_front_names(child0))
            name_set1 = set(get_front_names(child1))
            if len(name_set0.intersection(name_set1)) > 0:
                if len(name_set0) > len(name_set1):  # cut side w/ fewer taxa
                    node.remove_child(child1)
                    child1.prune()
                else:
                    node.remove_child(child0)
                    child0.prune()
                node, root = remove_kink(node, root)  # no rerooting here
                break
    return root


def prune_mo(tree_file, output_dir, min_taxa, out_groups):
    output_files = []

    # read in the tree and check number of taxa
    with open(tree_file) as infile:
        intree = newick3.parse(infile.readline())
    curroot = intree
    names = get_front_names(curroot)
    num_tips, num_taxa = len(names), len(set(names))
    if num_taxa < min_taxa:
        return output_files  # not enough taxa

    # If the homolog has no taxon duplication, no cutting is needed
    if num_tips == num_taxa:
        if OUTPUT_1TO1_ORTHOLOGS:
            output_file = util.file_name(
                output_dir, tree_file, '_1to1ortho.tre')
            copyfile(tree_file, output_file)
            output_files.append(output_file)
    else:
        # now need to deal with taxon duplications
        # check to make sure that the ingroup and outgroup names were
        # set correctly
        outgroup_names = get_front_outgroup_names(curroot, out_groups)

        # if no out-group at all, do not resolve gene duplication
        if len(outgroup_names) == 0:
            print("duplicated taxa in unrooted tree")

        # skip the homolog if there are duplicated out-group taxa
        elif len(outgroup_names) > len(set(outgroup_names)):
            print("outgroup contains taxon repeats")

        else:  # at least one out-group present and there's no out-group
            # duplication
            if curroot.nchildren == 2:  # need to reroot
                _, curroot = remove_kink(curroot, curroot)
            curroot = reroot_with_monophyletic_outgroups(curroot, out_groups)
            # only return one tree after pruning
            if curroot is not None:
                output_file = util.file_name(output_dir, tree_file, '.reroot')
                output_files.append(output_file)
                with open(output_file, "w") as outfile:
                    outfile.write(newick3.tostring(curroot) + ";\n")
                ortho = prune_paralogs_from_rerooted_homotree(
                    curroot, out_groups)
                if len(set(get_front_names(curroot))) >= min_taxa:
                    output_file = util.file_name(
                        output_dir, tree_file, '.ortho.tre')
                    output_file += '.ortho.tre'
                    output_files.append(output_file)
                    with open(output_file, "w") as outfile:
                        outfile.write(newick3.tostring(ortho) + ";\n")
                else:
                    print("not enough taxa after pruning")
            else:
                print("out-group non-monophyletic")

    return output_files
