import pandas as pd
from Bio import SeqIO
from itertools import permutations
from itertools import product
import time
import regex as re
import random
import os
import multiprocessing
from multiprocessing import Pool

def getK_mers(k):
    N="ATCG"
    L=product(N,repeat=k)
    t=[i for i in list(L)]
    k_mers=[]
    for i in t:
        s=''
        for j in i:
            s+=j
        k_mers.append(s)
    return k_mers


def getFeatureList(n_kmer):
    featureList=["ID"]
    for i in range(1,n_kmer):
        l=getK_mers(i)
        featureList+=l
    return featureList


def getTrainingset(rec):
    featureList=getFeatureList(5)
    Dic=[]
    Dic.append(str(rec.id))
    for i in featureList[1:]:
        seq=str(rec.seq)[20:-20].upper()
        Dic.append(len(re.findall(i, seq, overlapped=True)))
    rec_df=pd.DataFrame([Dic],columns=featureList)
    v=""
    for k in Dic[:-1]:
        v+=str(k)+","
    v+=str(Dic[-1])
    return Dic

#change
out=open("chr2_toPre.csv","a+")
s=""
featureList=getFeatureList(5)
for j in featureList[:-1]:
    s+=j+","
s+=featureList[-1]
out.write(s+"\n")
#change
file="IRF_chr2.fa"
if __name__ == '__main__':
    records = list(SeqIO.parse(file,"fasta"))
    print(len(records))
    pool = multiprocessing.Pool(16)
    d=pool.map(getTrainingset,records)
    pool.close()
    pool.join()
    print(len(d))
    for single in d:
        data=""
        for init in single[:-1]:
            data+=str(init)+","
        data+=str(single[-1])
        out.write(data+"\n")
    out.close()
