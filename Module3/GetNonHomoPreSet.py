
import pandas as pd
from Bio import SeqIO
import time
import random
import os


def getNonHomoTopre(chr):
    allData = pd.read_table("chr%s_toPre.csv"%(chr), header=0, sep=",")
    records = list(SeqIO.parse("chr%s_NonHomoSeq.fa" % (chr),"fasta"))
    NonhomoList = [str(rec.id) for rec in records]
    nonHomo = allData.loc[allData["ID"].isin(NonhomoList)]
    nonHomo.to_csv("chr%s_NonHomoToPre.csv"%(chr), index=None,sep=",")

for i in range(1,13):
    getNonHomoTopre(i)
