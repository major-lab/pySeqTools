import re

__author__ = 'nath'
from numpy import sum
from Sequence import Sequence


class fastaContainer:

    def __init__(self,file):
        self.end=False
        self.file = file
        self.f = open(self.file,'r')
        self.curName = self.f.readline().rstrip("\n").lstrip(">")

    def __del__(self):
        self.f.close()


    def nextSeq(self):
        """return the next sequence in the fasta file"""
        buffName = self.curName
        buffSeq = ""
        while ( 1 ):
            l = self.f.readline()
            if not l:
                self.end=True
                break
            if l[0:1] == '>':
                self.curName = l.rstrip("\n").lstrip(">")
                break
            else:
                buffSeq += l.rstrip("\n")

        res = Sequence(buffSeq,buffName)
        return res

    def resetFile(self):
        """reset the position of the file handler to the top of the file
        """
        self.f.close()
        self = self.__init__(self.file)



    def findAll(self,listIDs,mode="NM"):
        """return a list of sequences where the Ids are provided in the list
            mode can be set to "NM" (refseq IDs) or "MIMAT" (mirBase IDs)
        """
        self.resetFile()
        res = []

        tests = [0]*len(listIDs)

        while ( 1 ):
            s = self.nextSeq()

            if self.end or sum(tests) == len(listIDs):
                break
            id = ""
            if mode == "NM":
                id = s.parseNM()
            if mode == "MIMAT":
                id = s.parseMIMAT()
            if id in listIDs:
                tests[listIDs.index(id)] = 1
                res.append(s)


        return res

    def findSeq(self,ID,mode="NM"):
        self.resetFile()
        listTemp = []
        listTemp.append(ID)
        res = self.findAll(listTemp,mode)
        if (len(res) < 1 ):
            return None
        return res[0]

def revCompl(seq):
    res = seq[::-1].upper()
    transtab = "".maketrans("AUTCG","UAAGC")
    res= res.translate(transtab)
    return res.upper()





if __name__ == '__main__':
    fc = fastaContainer("/Users/nath//IRIC/paul/human.rna-NOTPREDICTED-mRNA.fasta")
    print(fc.nextSeq())
    print(fc.nextSeq())
    fc.resetFile()
    s = fc.nextSeq()
    print (s)
    print (s.parseNM())
    #print (fc.findAll(["NM_000314","NM_989898","NM_001145712"]))
    print("1")
    s  = fc.findSeq("NM_000314")
    print (s)
    tag = "*5(1-1031)c(1032-2243)3(2244-5572)"
    #tag = "*c(1-486)"

    print (tag)
    tab=[]

    """m = re.search(r"5\((\d+)-(\d+)\)c\((\d+)-(\d+)\)3\((\d+)-(\d+)\)",tag)

    for i in range(1,7):
        tab.append(int(m.group(i)))

    print (tab)
    tab = []"""

    """exemple = s.seq
    m = re.search(r"5\((\d+)-(\d+)\)",tag)
    if m:
        tab.append(int(m.group(1)))
        tab.append(int(m.group(2)))
    m = re.search(r"c\((\d+)-(\d+)\)",tag)
    if m:
        tab.append(int(m.group(1)))
        tab.append(int(m.group(2)))
    m = re.search(r"3\((\d+)-(\d+)\)",tag)
    if m:
        tab.append(int(m.group(1)))
        tab.append(int(m.group(2)))
    print (tab)
    seqs = []
    for i in range (0,int(len(tab)/2)):
        seqs.append(exemple[tab[i]:tab[i+1]])
    """
    print (s.seq)
    print (s.seq.upper())
    print (revCompl(s.seq))

    print ()
    print ('A'+revCompl('A'))




