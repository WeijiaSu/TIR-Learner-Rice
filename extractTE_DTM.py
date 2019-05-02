import pandas as pd
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO
import os
import multiprocessing
from multiprocessing import Pool
import time





def GetFastaFromFile(line):

    p1 = int(line.split("\t")[8])
    p2 = int(line.split("\t")[9])
    entry=int(line.split("\t")[1])
    o=open("DTM_fullcover_chr%s.fa"%(entry), 'a+')
    if (p1 < p2):
        p_start = p1-20
        p_end = p2+20
    else:
        p_start = p2-20
        p_end = p1+20
    seq1 = subprocess.check_output("blastdbcmd -db '%s' -entry '%s' -range '%s'-'%s'" % ("/work/LAS/thomasp-lab/weijia/research/Rice/Genome/ncbi/Rice.fa", entry, int(p_start), int(p_end)), shell=True)
    seq1 = seq1.decode("utf-8")
    out_seq1 = ''
    split = seq1.split('\n')
    for sp in split:
        if any([i.isdigit() for i in sp]):
            continue
        out_seq1 += sp
    o.write(">%s_%s_%s_%s_%s_%s" % (line.split("\t")[0],entry, p1, p2, p_start, p_end) + "\n" + str(out_seq1) + "\n")
    f.close()
    o.close()

# #
# #
if __name__ == '__main__':
     for i in range (1,13):
         f=open("%s_DTM_select.csv"%(i),"r+")
         lines=f.readlines()
         pool = multiprocessing.Pool(16)
         pool.map(GetFastaFromFile,lines)
         pool.close()
         pool.join()
