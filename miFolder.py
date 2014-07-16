from script.foldSeq import foldSeq, getE, foldSeqConstr, foldDupl
from script.ppPrint import getPPrint

__author__ = 'nath'

from script.seqTools import fastaContainer


def doItOn(listPos, mir, mRNA, verbose=False):
    res = []
    for pos in listPos:
        s = mRNA.getMRE(pos, mir.seq)
        offset = 30
        maskprefix = min(offset,pos)
        masksuffix = min(offset+1, len(mRNA.seq) - pos -6 )
        #print(maskprefix)
        #print (len(mRNA.seq) - pos-7)
        #print (masksuffix)
        mask = "x" * maskprefix + "x" * (len(mir.seq) - 8) + "." * 7 + "x" * masksuffix
        #pp = getPPrint(mRNA.seq, mir.seq, pos, 30)
        #print (pp)
        eFold = getE(foldSeq(s))
        eFoldConstr = getE(foldSeqConstr(s, mask))
        eFoldDuppl = foldDupl(s[offset + len(mir.seq) - 8:offset + len(mir.seq) - 1], mir.seq[1:8])
        eTot = eFoldConstr + eFoldDuppl - eFold
        if verbose:
            print("#" * 40 + " pos: " + str(pos) + " " + "#" * 40)
            pp = getPPrint(mRNA.seq, mir.seq, pos, 30)
            print(pp)
            print(s)
            print(foldSeq(s))
            print(s)
            print(mask)
            print(foldSeqConstr(s, mask))
            print(foldDupl(s, mir.seq, mode="-sd", getE=False))
            print("energy seed : " + str(eFoldDuppl))
            print("energy totale : " + str(eTot))
        res.append(eTot)
    return res


if __name__ == '__main__':
    fc = fastaContainer("../seq_p21.fa")
    mRNA = fc.findSeq("SM_000001")

    fc = fastaContainer("../data/miR_wu.fa")

    fPos = open("../data/posMir_Pos_sorted.txt")
    fout = open("res.txt","w")

    for l in fPos:
        l = l.rstrip("\n")
        buff = l.split("\t")
        #print (buff[0],buff[1])
        #mir = fc.findSeq("MIMAT0004614", "MIMAT")
        mir = fc.findSeq(buff[0], "MIMAT")
        listPos = [int(buff[1])]
        #listPos = [26]
        if int(buff[1])+8< len(mir.seq)+30:
            continue
        res = doItOn(listPos, mir,mRNA)
        fout.write (str(buff[0])+"\t"+str(buff[1])+"\t"+str(res[0])+"\n")
        #break
