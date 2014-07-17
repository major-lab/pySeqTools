from script.seqTools import *

__author__ = 'nath'


def scan(mRNA,miR):
    """
    scan a sequence of mRNA 5' -> 3' looking for a MRE. The seed is defined by the nt 2->7 of the miRNA
    return a dictionary
        position (int): type of seed (str)
    """

    seed = miR.seq[1:7]
    seed_rc = revCompl(seed)
    listHits = re.finditer(seed_rc,mRNA.seq)
    A1 = False
    m8 = False
    seedType = ""
    res = {}
    for m in listHits:
        print (mRNA.seq[m.start()-1:m.end()+1])
        if mRNA.seq[m.end()] == 'A':
            A1 = True
        if mRNA.seq[m.start()-1] == revCompl(miR.seq[7]):
            m8=True
        if A1 and m8:
            seedtype = "8mer"
            #in the case of 8mer, position on m8

        elif m8:
            seedtype = "7mer-m8"
            #in the case of 7mer-m8, position on m8

        elif A1:
            seedtype = "7mer-A1"
            #in the case of 7mer-A1, position on m7
        else:
            seedtype = "6mer"
        res[m.start()] = seedtype
    return res




if __name__ == '__main__':
    fc = fastaContainer("../seq_p21.fa")
    mRNA = fc.findSeq("SM_000001")
    fc = fastaContainer("../data/miR_wu.fa")
    mir = fc.findSeq("MIMAT0000063","MIMAT")
    print (mir.seq)
    res = scan(mRNA,mir)
    print (res)

