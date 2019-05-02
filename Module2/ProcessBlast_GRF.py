import pandas as pd
import os

def getFile(chr,fam):
    f=open("%s_GRF_%s_withhomology"%(chr,fam),"r+")
    o=open("%s_GRF_%s_withhomology_p"%(chr,fam),"a+")
    lines=f.readlines()
    for line in lines:
        if (line[0]!="#"):
            o.write(line)
    f.close()
    o.close()
    mv = "mv %s_GRF_%s_withhomology_p %s_GRF_%s_withhomology"%(chr,fam,chr,fam)

    os.system (mv)

for i in range(1,13):
    for j in ["DTA","DTC","DTH","DTM","DTT"]:
        getFile(i,j)


def getBlast(chr,fam):
    f=pd.read_table("%s_GRF_%s_withhomology"%(chr,fam),header=None,sep="\t")
    print(f.shape)
    f=f.loc[(f[11]>=80)&(f[3]>=80)]
    print(f[0:10])
    f=f.sort_values([0,11,3],ascending=[True,True,True])
    f=f.drop_duplicates([0],keep="last")
    f.to_csv("%s_GRF_%s_withhomology_80_unique"%(chr,fam),header=None,index=None,sep="\t")
    print(f.shape)

for i in range(1,13):
    for j in ["DTA","DTC","DTH","DTM","DTT"]:
        getBlast(i,j)
