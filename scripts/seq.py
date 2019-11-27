offset = 33  # 33 for most newer reads, 64 for old Illumina reads

COMPLEMENT = str.maketrans('ACGTUWSMKRYBDHVNXacgtuwsmkrybdhvnx-',
                           'TGCAAWSKMYRVHDBNXtgcaawskmyrvhdbnx-')


class Sequence:
    def __init__(self, name="", seq=""):
        self.name = name
        self.seq = seq
        self.qualstr = ""  # string of quality scores in characters
        self.qualarr = []  # list of quality scores in ASCII numbers

    # input a string of quality characters, set both qualstr and qualarr
    def set_qualstr(self, qual):
        self.qualstr = qual
        if len(self.qualarr) == 0:
            for j in self.qualstr:
                quality_score = ord(j) - offset
                self.qualarr.append(quality_score)
                assert 0 < quality_score <= 41, \
                    "Change the offset in seq.py\nqual"

    def get_fastq(self):
        """Return the record in fastq format."""
        return '@{}\n{}\n+\n{}\n'.format(self.name, self.seq, self.qualstr)

    def get_fasta(self):
        """Return the record in fasta format."""
        return '>{}\n{}\n'.format(self.name, self.seq)

    def rev_comp(self):
        """Reverse complement a nucleotide sequence."""
        return self.seq.translate(COMPLEMENT)[::-1]


def fastq_generator(infile):
    line = infile.readline()
    while len(line) > 0:  # no empty string at end-of-file in readline
        if line[0] == "@":
            name = line[1:].strip()  # name of sequence, minus "@"
            seq = infile.readline().strip()  # the actual sequence
            line = infile.readline().strip()  # "+"
            qual = infile.readline().strip()  # the quality string
            tseq = Sequence(name=name, seq=seq)
            tseq.set_qualstr(qual)
            yield tseq
        line = infile.readline()


def read_fasta_file(file_name):
    """given the path to a fasta file, return a list of seq objects"""
    fl = open(file_name, "r")
    seqlist = []  # list of sequence objects
    templab = ""
    tempseq = ""
    first = True
    for i in fl:
        if i[0] == ">":
            if first:
                first = False
            else:  # if not first store lastseq read
                seqlist.append(Sequence(templab, tempseq))
            templab = i.strip()[1:]
            tempseq = ""
        else:
            tempseq = tempseq + i.strip()
    fl.close()
    seqlist.append(Sequence(templab, tempseq))
    return seqlist
