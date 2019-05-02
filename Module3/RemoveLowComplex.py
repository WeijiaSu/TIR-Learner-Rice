from Bio import SeqIO
import multiprocessing
from multiprocessing import Pool

def TArepeats(s):
    t=s.upper().count("T")
    a=s.upper().count("A")
    ta=t+a
    if (ta>len(s)*0.7):
        return True
    else:
        return False


def checkN(s):
    n=s.upper().count("N")
    if n>0:
        return True
    else:
        return False


def getSeqID(file):
    remove=[]
    records=list(SeqIO.parse(file,"fasta"))
    for rec in records:
        s=str(rec.seq)[0:10]
        if (TArepeats(s)==True or checkN(s)==True):
            remove.append(rec.id)
    records=list(SeqIO.parse(file,"fasta"))
    SeqIO.write((rec for rec in records if rec.id not in remove),file+"_p","fasta")

if __name__ == '__main__':
    l=[]
    for i in range(1,13):
        l.append("IRF_chr%s.fa"%(i))
    pool = multiprocessing.Pool(16)
    pool.map(getSeqID,l)
    pool.close()
    pool.join()


