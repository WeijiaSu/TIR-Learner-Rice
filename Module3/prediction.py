import pandas as pd
from Bio import SeqIO
import os
import multiprocessing
from multiprocessing import Pool



def getPredictedSeq(chr):
    f = pd.read_table("prediction_chr_%s.csv"%(chr), header=0, sep=",")
    records=list(SeqIO.parse("chr%s_NonHomoSeq.fa"%(chr),"fasta"))
    for i in ["DTA","DTC","DTH","DTM","DTT"]:
        print(str(chr)+": "+i)
        sub=f.loc[f["prediction"]==i]
        preNamelist=list(sub["ID"])
        SeqIO.write((seq for seq in records if str(seq.id) in preNamelist), "chr%s_pred_%s_Seq.fa"%(chr,i),"fasta")



if __name__ == '__main__':
     l=list(range(1,13))
     pool = multiprocessing.Pool(13)
     pool.map(getPredictedSeq,l)
     pool.close()
     pool.join()
