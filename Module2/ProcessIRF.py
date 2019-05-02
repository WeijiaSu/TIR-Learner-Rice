import pandas as pd
import os
from Bio import SeqIO


def ProcessDatfile(l):
    for i in l:
        filename=i+".fa.1.5.10.80.40.10.10000.10000.dat"
        f=pd.read_table(filename,header=None,sep=" ",skiprows=[0,1,2,3,4,5,6,7,8,9,10,11,12])
        if (f.shape[1]==19):
            f = f.drop([9, 10, 11, 12, 13, 14, 15, 16,17,18], axis=1)
            f = f.loc[f[7] >= 80]
            f = f.loc[f[5] == f[2]]
            f = f.sort_values([0, 1, 3, 4, 7])
            f = f.loc[(f[6] >= 50) & ((f[6] + f[5] + f[2]) <= 10000)]
            f.columns = ["Left_1", "Left_2", "Length", "Right_1", "Right_2", "Length", "Loop", "Matches", "Indels"]
            f.to_csv(filename + ".csv", index=None, sep=",")
        else:
            f = f.drop(0, axis=1)
            f.columns = list(range(0, 17))
            f = f.drop([9, 10, 11, 12, 13, 14, 15, 16], axis=1)
            f = f.loc[f[7] >= 80]
            f = f.loc[f[5] == f[2]]
            f = f.sort_values([0, 1, 3, 4, 7])
            f = f.loc[(f[6] >= 50) & ((f[6] + f[5] + f[2]) <= 10000)]
            f.columns = ["Left_1", "Left_2", "Length", "Right_1", "Right_2", "Length", "Loop", "Matches", "Indels"]
            f.to_csv(filename + ".csv", index=None, sep=",")



def ProcessIRF_combine(l,chr):
    heaher = ["Left_1", "Left_2", "Length", "Right_1", "Right_2", "Length.1", "Loop", "Matches", "Indels"]
    Total = pd.DataFrame(columns=heaher)
    for i in l:
        filename=i+".fa.1.5.10.80.40.10.10000.10000.dat.csv"
        f = pd.read_table(filename, header=0, sep=",")
        Total=Total.append(f,ignore_index=True)
    Total=Total.sort_values(["Left_1","Left_2","Length","Loop","Matches"],ascending=[True,True,True,True,True])
    Total=Total.drop_duplicates("Left_1",keep="last")
    outname=str(chr)+"_IFR_selected.csv"
    Total.to_csv(outname,index=None,sep=",")


os.chdir("/work/LAS/thomasp-lab/weijia/research/Rice/IRF/")
for chr in range(1,13):
    l=["Rice_Chr%s"%(chr)]
    ProcessDatfile(l)
    ProcessIRF_combine(l,chr)



