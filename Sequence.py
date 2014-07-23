__author__ = 'nath'
import re

class Sequence:


    def __init__(self,seq,name):
        self.name = name
        self.seq = seq
        self.ID=""
        self.CDS=-1
        self.UTR3=-1


    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return ">"+self.name+"\n"+self.seq

    def parseNM(self):
        """return the refSeq ID contained in the name"""
        m = re.search(r'.._\d+',self.name)
        if m != None:
            self.ID = m.group(0)
            return self.ID
        return ""


    def parseMIMAT(self):
        """return the mirBase ID contained in the name of the sequence"""
        m = re.search(r"MIMAT\d+",self.name)
        if m != None:
            self.ID = m.group(0)
            return self.ID
        return ""

    def getMRE(self,start,miRNA,offset=30,seedLength = 7 ):
        return self.seq[start + seedLength - len(miRNA) - offset:start + seedLength + offset]

if __name__ == '__main__':
    s = Sequence("toto","tata")
    print (s)

