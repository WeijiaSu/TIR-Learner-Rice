
import pandas as pd
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO
import os
import multiprocessing
from multiprocessing import Pool
import time


TIR={}
TIR["DTA"]=9
TIR["DTT"]=9
TIR["DTM"]=9
TIR["DTH"]=9
TIR["DTC"]=9
#
TSD={}
TSD["DTA"]=8
TSD["DTT"]=2
TSD["DTM"]=9
TSD["DTH"]=3
TSD["DTC"]=3



def getLengthDic(DATfile):
    TIRlength = {}
    f=pd.read_table(DATfile,header=0,sep=",")
    print(f.shape)
    for index,row in f.iterrows():
        cor=str(int(row["Left_1"]))+":"+str(int(row["Right_2"]))
        TIRlength[cor]=int(row["Length"])
    return TIRlength

def get_seq(seq):
    out = ''
    split = seq.split('\n')
    for sp in split:
        if any([i.isdigit() for i in sp]):
            continue
        out += sp
    return out

def compare(s1,s2):
    d=0
    for i in range(0,len(s1)):
        if(s1[i]!=s2[i]):
            d+=1
    return d


def SeqtoCompare(rec,TIRlength,Family):
    seq = str(rec.seq)
    seq = get_seq(seq)
    seqID = str(rec.id)
    startp = int(seqID.split("_")[1])
    endp=int(seqID.split("_")[2])
    cor = str(startp) + ":" + str(endp)
    tir_l = TIRlength[cor]
    tir_limit=TIR[Family]
    tsd_limit=TSD[Family]
    if tir_l > tir_limit:
        s1 = seq[20-tsd_limit:20] + seq[20:20 + tir_l - tir_limit]
        s2 = seq[-20 - (tir_l - tir_limit):-(20-tsd_limit)]
    else:
        s1 = seq[20-tsd_limit:20]
        s2 = seq[-20:-(20-tsd_limit)]
    return s1,s2

def CompareTSDs(S1,S2,Family):
    l=len(S1)
    tsd={}
    tsd_limit=TSD[Family]
    for i in range(0,l-tsd_limit+1):
        seq1=S1[l-i-tsd_limit:l-i]
        seq2=S2[i:i+tsd_limit]
        d=compare(seq1,seq2)
        if d<TSD[Family]*0.2:
            tsd[i]=d
    return tsd

def Re_writeSeq(tsdDIC,rec,outputname,TIRlength,Family):
    if (len(tsdDIC)!=0):
        for i in tsdDIC:
            tsdStart=i
            seq = str(rec.seq)
            seq = get_seq(seq)
            seqID = str(rec.id)
            startp = int(seqID.split("_")[1])
            endp = int(seqID.split("_")[2])
            chro=seqID.split("_")[0]
            cor = str(startp) + ":" + str(endp)
            tir_l = TIRlength[cor]
            tir_limit=TIR[Family]
            tsd_limit=TSD[Family]
            if tir_l > tir_limit:
                s1 = seq[20-tsd_limit:20] + seq[20:20 + tir_l - tir_limit]
                s2 = seq[-20 - (tir_l - tir_limit):-(20-tsd_limit)]
                l = len(s1)
                TSD1 = s1[l - tsdStart - tsd_limit:l - tsdStart]
                TSD2 = s2[tsdStart:tsdStart + tsd_limit]
                comSeq = TSD1 + s1[l - tsdStart:l] + seq[20 + tir_l - tir_limit:-20 - (tir_l - tir_limit)] + s2[
                                                                                                             0:tsdStart] + TSD2
                coor1 = startp + tir_l - tir_limit - tsdStart
                coor2 = endp - (tir_l - tir_limit) + tsdStart
                newID = chro +"_"+str(coor1) + "_" + str(coor2) + "_" + str(coor1 - tsd_limit) + "_" + str(coor2 + tsd_limit)
                out = open(outputname, "a+")
                out.write(">" + newID+"_"+Family + "\n" + comSeq + "\n")
                out.close()
            else:
                s1 = seq[20-tsd_limit:20]
                s2 = seq[-20:-(20-tsd_limit)]
                TSD1=s1
                TSD2=s2
                comSeq=TSD1+seq[20:-20]+TSD2
                coor1=startp
                coor2=endp
                newID = chro+"_"+str(coor1) + "_" + str(coor2) + "_" + str(coor1 - tsd_limit) + "_" + str(coor2 + tsd_limit)
                out = open(outputname, "a+")
                out.write(">" + newID+"_"+Family + "\n" + comSeq + "\n")
                out.close()



def GetDictionary(uniqueHomolist):
    dic={}
    f=pd.read_table(uniqueHomolist,header=None,sep="\t")
    for index, row in f.iterrows():
        id=row[0]
        family=row[1].split("_")[-1]
        dic[id]=family
    return dic

for chromosome in range(1,13):
    dic=GetDictionary("%s_IRF_withhomology_80_unique"%(chromosome))
    TIRlength=getLengthDic("%s_IFR_selected.csv"%(chromosome))
    records=list(SeqIO.parse("chr%s_HomoSeq.fa"%(chromosome),"fasta"))
    v=[k for k in dic]
    n = 0
    l = 0
    i=0
    reco=[rec for rec in records if rec.id in v]
    for rec in reco:
        i+=1
        print(i)
        id=str(rec.id)
        family=dic[id]
        S1, S2 = SeqtoCompare(rec, TIRlength, family)
        tsdDIC = CompareTSDs(S1, S2, family)
        if (len(tsdDIC) != 0):
            n += 1
            l += len(tsdDIC)
            Re_writeSeq(tsdDIC, rec, "chr%s_AllcheckedHomo.fa"%(chromosome), TIRlength, family)
    print(n)
    print(l)

