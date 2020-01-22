import sys
from shlex import shlex
from pylib.phylo3 import Node
from io import StringIO


class Tokenizer(shlex):
    """Provides tokens for parsing Newick-format trees"""
    def __init__(self, infile):
        shlex.__init__(self, infile)
        self.commenters = ''
        self.wordchars = self.wordchars+'-.@'
        self.quotes = "'"

    def parse_comment(self):
        while 1:
            token = self.get_token()
            if token == '':
                sys.stdout.write('EOF encountered mid-comment!\n')
                break
            elif token == ']':
                break
            elif token == '[':
                self.parse_comment()
            else:
                pass


def parse(raw, ttable=None):
    """
    Parse a Newick-formatted tree description
    input is any file-like object that can be coerced into shlex,
    or a string (converted to StringIO)
    """
    if type(raw) is str:
        raw = StringIO(raw)

    start_pos = raw.tell()
    tokens = Tokenizer(raw)

    node = None
    lp = 0
    rp = 0

    prev_tok = None

    while 1:
        token = tokens.get_token()
        if token == ';' or token == '':
            assert lp == rp, \
                   'unbalanced parentheses in tree description'
            break

        # internal node
        elif token == '(':
            lp = lp+1
            new_node = Node()
            new_node.istip = False
            if node:
                node.add_child(new_node)
            node = new_node

        elif token == ')':
            rp = rp + 1
            node = node.parent

        elif token == ',':
            node = node.parent

        # branch length
        elif token == ':':
            token = tokens.get_token()

            if not (token == ''):
                try:
                    brlen = float(token)
                except ValueError:
                    raise Exception(
                        "NewickError invalid literal for branch "
                        "length, '%s'" % token)
            else:
                raise Exception(
                    "NewickError unexpected end-of-file "
                    "(expecting branch length)")

            node.length = brlen
        # comment
        elif token == '[':
            tokens.parse_comment()

        # leaf node or internal node label
        else:
            if prev_tok != ')': # leaf node
                if ttable:
                    ttoken = ttable.get(token) or ttable.get(int(token))
                    if ttoken:
                        token = ttoken
                newnode = Node()
                newnode.label = token
                newnode.istip = True
                node.add_child(newnode)
                node = newnode
            else: # label
                # translation table for internal nodes labels?
                node.label = token

        prev_tok = token

    raw.seek(start_pos)

    return node


def traverse(node):
    if node.istip:
        return node.back
    else:
        return node.next.back


def to_string(node, length_fmt=":%s"):
    if not node.istip:
        node_str = "(%s)%s" % \
                   (",".join([ to_string(child, length_fmt) \
                               for child in node.children ]),
                    node.label or ""
                    )
    else:
        node_str = "%s" % node.label

    if node.length is not None:
        length_str = length_fmt % node.length
    else:
        length_str = ""

    s = "%s%s" % (node_str, length_str)
    return s


tostring = to_string


def parse_from_file(filename):
    if filename == '-':
        file = sys.stdin
    else:
        file = open(filename, 'r')
    content = file.read()
    content = content.strip()
    treedescs = content.split(';')
    tree = parse(treedescs[0])
    file.close()
    return tree
