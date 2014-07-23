from pip.util import dist_in_site_packages
from script.ppPrint import doPPrint
from script.seqTools import fastaContainer
from script.ScanSeed import *

__author__ = 'nath'



if __name__ == '__main__':
    fc = fastaContainer("../seq_p21.fa")
    mRNA = fc.findSeq("SM_000001")
    mRNA.CDS = 10
    mRNA.UTR3 = 947



    fc = fastaContainer("../data/miR_wu.fa")
    fout = open("res_scanSeed.txt","w")
    print ("mimat\tmer6\tmer7-m8\tmer7-A1\tmer8\tdist6MerMin")
    fout.write("mimat\tmer6\tmer7-m8\tmer7-A1\tmer8\tdist6MerMin\n")
    for mir in fc:

        res = scan(mRNA,mir)
        print (mir.name+"\t", end ='')
        fout.write(mir.name+"\t")
        fp = [0,0,0,0]
        dist6Mer = 100000
        lastPos6mer=0
        tabPos6mer = []
        for pos in res.keys():
            if pos >= mRNA.UTR3:
                #print ('',pos,res.get(pos))
                #doPPrint(mRNA.seq,mir.seq,pos)
                #if res.get(pos)== "6mer":
                tabPos6mer.append(pos)

                newFp = getFPScan(res.get(pos))

                fp = sumFPScan(fp,newFp)
        if len(tabPos6mer) > 1:

            for i in range(0,len(tabPos6mer)):
                for j in range(0,len(tabPos6mer)):
                    if i==j :
                        continue
                    else :
                        dist6Mer = min(dist6Mer,abs(tabPos6mer[i]-tabPos6mer[j]))

        print ("\t".join(map(str,fp))+"\t"+str(dist6Mer))
        fout.write("\t".join(map(str,fp))+"\t"+str(dist6Mer)+"\n")
