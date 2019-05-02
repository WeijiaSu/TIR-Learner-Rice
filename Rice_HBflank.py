import pandas as pd
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO
import os
import multiprocessing
from multiprocessing import Pool
import time




TIR={}
TIR["DTA"]=11
TIR["DTT"]=13
TIR["DTM"]=40
TIR["DTH"]=14
TIR["DTC"]=13
#
TSD={}
TSD["DTA"]=8
TSD["DTT"]=2
TSD["DTM"]=9
TSD["DTH"]=3
TSD["DTC"]=3


def GetFastaFromFile(line):

    p1 = int(line.split("\t")[3])
    p2 = int(line.split("\t")[4])
    entry=int(line.split("\t")[0])
    family=str(line.split("\t")[2])
    o=open("Rice_HB200F.fa", 'a+')
    TSD_l=TSD[family]
    if (p1 < p2):
        p_start = p1-200
        p_end = p2+200
    else:
        p_start = p2-200
        p_end = p1+200
    seq1 = subprocess.check_output("blastdbcmd -db '%s' -entry '%s' -range '%s'-'%s'" % ("/work/LAS/thomasp-lab/weijia/research/Rice/Genome/ncbi/Rice.fa", entry, int(p_start), int(p_end)), shell=True)

    seq1 = seq1.decode("utf-8")
    out_seq1 = ''
    split = seq1.split('\n')
    for sp in split:
        if any([i.isdigit() for i in sp]):
            continue
        out_seq1 += sp
    o.write(">%s_%s_%s_%s_%s_%s_%s" % ("Rice", entry, p1, p2, p_start, p_end,family) + "\n" + str(out_seq1) + "\n")
    f.close()
    o.close()

# #
# #
if __name__ == '__main__':
     f=open("Rice_RefHB.gff3","r+")
     lines=f.readlines()
     pool = multiprocessing.Pool(16)
     pool.map(GetFastaFromFile,lines)
     pool.close()
     pool.join()
