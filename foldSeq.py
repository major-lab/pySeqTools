__author__ = 'weilln'

import subprocess,re
dir = "/Users/nath/IRIC/flashFold/FlashFold"
bin = dir + "/macOSXbin/mcff"
table = dir + "/tables"




def foldSeq(seq):
    res = subprocess.check_output([bin,"-s",seq,"-tables",table],universal_newlines=True)
    res = res.rstrip("\n")
    return res

def foldSeqConstr(seq,constr):
    res = subprocess.check_output(
        [bin,"-s",seq,"-tables",table,"-m",constr],
        universal_newlines=True)
    res = res.rstrip("\n")
    return res

def foldDupl(seq1,seq2,mode="-sd",getE=True):
    res = subprocess.check_output([bin,"-s",seq1,mode,seq2,"-tables",table, "-v"],universal_newlines=True,stderr=subprocess.STDOUT)
    if not getE:
        return res
    m =re.search(r'mfe\((-\d+\.\d+)',res)
    if m != None:
        return float(m.group(1))
    return ""


def getE(fold):
    buff = re.split('\s+', fold)
    return float(buff[1])


if __name__ == '__main__':
    seq1 = "AAAAAAGGGTTTTTTTT"
    seq2 = "GGGTTTAAATAAAAAAA"
    res = foldSeq(seq1)
    print (getE(res))
    res = foldSeqConstr("AAAAAAGGGTTTTTTTT",
                        "xxx..xxxxxxxxxxxx")
    print (getE(res))

    res = foldDupl(seq1,seq2)
    print (res)
