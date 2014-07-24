from script.ScanSeed import scan
from script.miFolder import *

__author__ = 'nath'
from script.seqTools import fastaContainer



if __name__ == '__main__':
    fc = fastaContainer("../seq_p21.fa")
    mRNA = fc.findSeq("SM_000001")
    mRNA.CDS = 10
    mRNA.UTR3 = 947


    f = open( "res_foldScanSeed_miR_foldmRNA.txt","w")
    fc = fastaContainer("../data/miR_wu.fa")

    for mir in fc:
        print (mir.name,end="")
        f.write(mir.name)
        res = scan(mRNA,mir)
        folds = foldSeedTypeMRNA(res,mir,mRNA,)
        E = 0
        for k in folds.keys():

            E+=folds[k]
            #print (mir.name+"\t"+res[k]+"\t"+str(folds[k])+"\t"+str(k))

        print("\t"+str(E))
        f.write("\t"+str(E)+"\n")
    f.close()
