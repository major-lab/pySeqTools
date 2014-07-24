from script.foldSeq import foldSeq, getE, foldSeqConstr, foldDupl
from script.ppPrint import getPPrint, doPPrint

__author__ = 'nath'

from script.seqTools import fastaContainer


def foldSeedMRNA(listPos, mir, mRNA, verbose=False):
    res = []
    offset = 30
    for pos in listPos:
        if pos < len(mir.seq)+offset :
            continue
        s = mRNA.getMRE(pos, mir.seq)

        ##MASK:
        maskprefix = min(offset,pos)
        masksuffix = min(offset+1, len(mRNA.seq) - pos -6 )
        mask = "x" * maskprefix + "x" * (len(mir.seq) - 8) + "." * 7 + "x" * masksuffix
        ##END MASK

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

def foldSeedTypeMRNA(listPos, mir, mRNA,UTR3 = True):
    res = dict()
    offset = 30
    seedLength = 7
    for pos in listPos.keys():
        if pos < len(mir.seq)+offset :
            continue

        if UTR3 and pos <  mRNA.UTR3  :
            continue
        if listPos.get(pos) == "6mer" or listPos.get(pos) == "7mer-A1" :
            seedLength = 6

        elif listPos.get(pos) == "8mer" or listPos.get(pos) == "7mer-m8" :
            seedLength = 7
        s = mRNA.getMRE(pos, mir.seq,seedLength=seedLength)
        ##MASK:
        maskprefix = min(offset,pos)
        masksuffix = min(offset+1, len(mRNA.seq) - pos - (seedLength-2) )
        mask = "x" * maskprefix + "x" * (len(mir.seq) - (seedLength+1)) + "." * seedLength + "x" * masksuffix
        ##END MASK

        eFold = getE(foldSeq(s))


        eFoldConstr = getE(foldSeqConstr(s, mask))
        tempMRNA =  mRNA.getMRE(pos, mir.seq[1:seedLength+1],seedLength=seedLength,offset=0)# s[offset + len(mir.seq) - seedLength -1:offset + len(mir.seq) -  1]


        eFoldDuppl = foldDupl(tempMRNA, mir.seq[1:seedLength+1])
        eTot = eFoldConstr + eFoldDuppl - eFold

        res[pos] = eTot
    return res


def foldSeed(listPos, mir, mRNA,UTR3 = True):
    res = dict()
    offset = 30
    seedLength = 7
    for pos in listPos.keys():
        if pos < len(mir.seq)+offset :
            continue

        if UTR3 and pos <  mRNA.UTR3  :
            continue
        if listPos.get(pos) == "6mer" or listPos.get(pos) == "7mer-A1" :
            seedLength = 6

        elif listPos.get(pos) == "8mer" or listPos.get(pos) == "7mer-m8" :
            seedLength = 7
        s = mRNA.getMRE(pos, mir.seq,seedLength=seedLength)
        ##MASK:
        maskprefix = min(offset,pos)
        masksuffix = min(offset+1, len(mRNA.seq) - pos - (seedLength-2) )
        mask = "x" * maskprefix + "x" * (len(mir.seq) - (seedLength+1)) + "." * seedLength + "x" * masksuffix
        ##END MASK

        #eFold = getE(foldSeq(s))


        #eFoldConstr = getE(foldSeqConstr(s, mask))
        tempMRNA =  mRNA.getMRE(pos, mir.seq[1:seedLength+1],seedLength=seedLength,offset=0)# s[offset + len(mir.seq) - seedLength -1:offset + len(mir.seq) -  1]


        eFoldDuppl = foldDupl(tempMRNA, mir.seq[1:seedLength+1])

        #eTot = eFoldConstr + eFoldDuppl - eFold

        res[pos] = eFoldDuppl
    return res
def foldMir(listPos, mir, mRNA,UTR3 = True):
    res = dict()
    offset = 30
    seedLength = 7
    for pos in listPos.keys():
        if pos < len(mir.seq)+offset :
            continue

        if UTR3 and pos <  mRNA.UTR3  :
            continue
        if listPos.get(pos) == "6mer" or listPos.get(pos) == "7mer-A1" :
            seedLength = 6

        elif listPos.get(pos) == "8mer" or listPos.get(pos) == "7mer-m8" :
            seedLength = 7
        s = mRNA.getMRE(pos, mir.seq,seedLength=seedLength)
        ##MASK:
        maskprefix = min(offset,pos)
        masksuffix = min(offset+1, len(mRNA.seq) - pos - (seedLength-2) )
        mask = "x" * maskprefix + "x" * (len(mir.seq) - (seedLength+1)) + "." * seedLength + "x" * masksuffix
        ##END MASK

        #eFold = getE(foldSeq(s))


        #eFoldConstr = getE(foldSeqConstr(s, mask))
        tempMRNA = s[offset -30:offset + len(mir.seq) +1]
        #print()
        #print (tempMRNA)
        #print (mir.seq[::-1])
        #print(foldDupl(tempMRNA, mir.seq,getE=False))

        eFoldDuppl = foldDupl(tempMRNA, mir.seq)

        #eTot = eFoldConstr + eFoldDuppl - eFold

        res[pos] = eFoldDuppl
    return res


if __name__ == '__main__':
    fc = fastaContainer("../seq_p21.fa")
    mRNA = fc.findSeq("SM_000001")

    fc = fastaContainer("../data/miR_wu.fa")

    fPos = open("../data/posMir_Pos_sorted.txt")
    #fout = open("res.txt","w")

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
        res = foldSeedMRNA(listPos, mir,mRNA)
        print (str(buff[0])+"\t"+str(buff[1])+"\t"+str(res[0]))
        #fout.write (str(buff[0])+"\t"+str(buff[1])+"\t"+str(res[0])+"\n")
        #break
