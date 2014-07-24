__author__ = 'weilln'
import sys, re


def readmiRNA(file):
    dico = {}
    f = open(file, 'r')
    key = ""

    for l in f:
        l = l.rstrip("\n")
        if (l[0] == '>'):
            #print (l)
            tab = re.split('\s+', l)
            #print (tab)
            key = tab[1]
        else:
            seq = l
            dico[key] = seq
    return dico


def getMatches(seq1, seq2, offsetSeq1=0, offsetSeq2=0):
    res = ""
    mySeq1 = seq1.replace("T","U")
    mySeq2 = seq2.replace("T","U")
    for i in range(0, min(len(seq1),len(seq2))):

        if mySeq1[i + offsetSeq1] == 'A' and mySeq2[i + offsetSeq2] == 'U':
            res += "|"
        elif mySeq1[i + offsetSeq1] == 'U' and mySeq2[i + offsetSeq2] == 'A':
            res += "|"
        elif mySeq1[i + offsetSeq1] == 'G' and mySeq2[i + offsetSeq2] == 'C':
            res += "|"
        elif mySeq1[i + offsetSeq1] == 'C' and mySeq2[i + offsetSeq2] == 'G':
            res += "|"
        else:
            res += "."
    return res

def getPPrint(mRNA,miRNA,mRNAStart,offsetmRNA = 30,seedLength = 7):
    if offsetmRNA <0 :
        offsetmRNA = 0
    res = ""

    seq1 = mRNA[mRNAStart + seedLength - len(miRNA) - offsetmRNA : mRNAStart + seedLength + offsetmRNA]
    seq2 = miRNA[::-1]
    res += ("5' - " + seq1 + " mRNA "+str(mRNAStart)+'\n')
    res += (offsetmRNA * ' ' + "     " + getMatches(seq1.upper(), seq2.upper(), offsetmRNA, 0)+'\n')
    res += (offsetmRNA * ' ' + "3' - " + seq2+'\n')
    return res




def doPPrint(mRNA,miRNA,mRNAStart,offsetmRNA = 30,seedLength = 7):
    res = getPPrint(mRNA, miRNA, mRNAStart, offsetmRNA, seedLength)
    print (res)



if __name__ == '__main__':

    #miRNAFile:
    miRNAFile = "../data/mir_base_mature_hsa.fa"
    dico = readmiRNA(miRNAFile)
    seqFile = sys.argv[1]
    miRNA = dico[sys.argv[2]]
    startSeq = int(sys.argv[3])

    #test:
    offset = 30
    f = open(seqFile, 'r')
    seq = f.readline()
    seq = seq.rstrip("\n")
    print(getPPrint(seq,miRNA,startSeq,30))
    doPPrint(seq,miRNA,startSeq,30)

    seq = "sdfsdfafda\n"
    seq.rstrip("\n")
    print (seq)




