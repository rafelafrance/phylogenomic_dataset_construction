"""
Mask both mono- and paraphyletic tips that belong to the same taxon
If only mask monophyletic tips, comment out this line:
Keep the tip that has the most un-ambiguous, well-aligned characters in the
trimmed alignment
"""

import re
import pylib.bio as bio
import pylib.newick3 as newick3
from .tree_utils import get_name, remove_kink
from pylib import util


IGNORE = re.compile(r'[x*?\-]', re.IGNORECASE)


def mask_tips(fasta_file, tree_file, output_dir, mask_paraphyletic):
    """Wrap tree tip removal."""
    with open(tree_file) as in_file:
        in_tree = newick3.parse(in_file.readline())

    chars = {k: len(IGNORE.sub('', v))
             for k, v in bio.read_fasta(fasta_file).items()}

    output = util.file_name(tree_file, '.mm')

    root = mask_monophyletic_tips(in_tree, chars)

    if mask_paraphyletic:
        root = mask_paraphyletic_tips(root, chars)

    with open(output, 'w') as outfile:
        outfile.write(newick3.tostring(root) + ";\n")

    return output


def mask_monophyletic_tips(root, chars):
    going = True
    while going and len(root.leaves()) >= 4:
        going = False
        for node in root.iternodes():  # walk through nodes
            if not node.istip:
                continue  # only look at tips
            for sister in node.get_sisters():
                if sister.istip and get_name(node.label) == get_name(
                        sister.label):  # masking
                    if chars[node.label] > chars[sister.label]:
                        node = sister.prune()
                    else:
                        node = node.prune()
                    if len(root.leaves()) >= 4:
                        if (node == root and node.nchildren == 2) or (
                                node != root and node.nchildren == 1):
                            node, root = remove_kink(node, root)
                    going = True
                    break
    return root


def mask_paraphyletic_tips(root, chars):
    going = True
    while going and len(root.leaves()) >= 4:
        going = False
        for node in root.iternodes():  # walk through nodes
            if not node.istip:
                continue  # only look at tips
            parent = node.parent
            if root in (node, parent):
                continue  # no paraphyletic tips for the root
            for para in parent.get_sisters():
                if para.istip and get_name(node.label) == get_name(para.label):
                    if chars[node.label] > chars[para.label]:
                        node = para.prune()
                    else:
                        node = node.prune()
                    if len(root.leaves()) >= 4:
                        if (node == root and node.nchildren == 2) or (
                                node != root and node.nchildren == 1):
                            node, root = remove_kink(node, root)
                    going = True
                    break
    return root
