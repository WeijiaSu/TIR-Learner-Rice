import pandas as pd
import os

def getFile(chr):
    f=open("%s_IRF_withhomology"%(chr),"r+")
    o=open("%s_IRF_withhomology_p"%(chr),"a+")
    lines=f.readlines()
    for line in lines:
        if (line[0]!="#"):
            o.write(line)
    f.close()
    o.close()
    mv = "mv %s_IRF_withhomology_p %s_IRF_withhomology"%(chr,chr)

    os.system (mv)

for i in range(1,13):
    getFile(i)


def getBlast(chr):
    f=pd.read_table("%s_IRF_withhomology"%(chr),header=None,sep="\t")
    print(f.shape)
    f=f.loc[(f[11]>=80)&(f[3]>=80)]
    print(f[0:10])
    f=f.sort_values([0,11,3],ascending=[True,True,True])
    f=f.drop_duplicates([0],keep="last")
    f.to_csv("%s_IRF_withhomology_80_unique"%(chr),header=None,index=None,sep="\t")
    print(f.shape)

for i in range(1,13):
    getBlast(i)
