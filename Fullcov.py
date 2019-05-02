import pandas as pd
import multiprocessing
from multiprocessing import Pool
import os


def getFile():
    for i in ["DTA","DTC","DTH","DTM","DTT"]:
        f=open("Blast_Rice_%s"%(i),"r+")
        o=open("Blast_Rice_%s_p"%(i),"a+")
        lines=f.readlines()
        for line in lines:
            if (line[0]!="#"):
                o.write(line)
        f.close()
        o.close()
        mv = "mv Blast_Rice_%s_p Blast_Rice_%s" %(i,i)
        os.system (mv)
       
getFile()


def ProcessHomology(chr):
    for i in ["DTA","DTC","DTH","DTM","DTT"]:
        blast="Blast_Rice_%s"%(i)
        f=pd.read_table(blast,header=None,sep="\t")
        f=f.loc[f[1]==chr]
        f=f.loc[(f[11]==100) & (f[3]>=80)]
        f=f.sort_values([1,8,9,11,3],ascending=[True,True,True,True,True])
        f=f.drop_duplicates([1,8,9],keep="last")
        f.to_csv("%s_%s_select.csv" % (chr, i), header=None, index=None, sep="\t")


if __name__ == '__main__':
     lines=list(range(1,13))
     pool = multiprocessing.Pool(12)
     pool.map(ProcessHomology,lines)
     pool.close()
     pool.join()
