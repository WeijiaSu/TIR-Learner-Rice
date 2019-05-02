import pandas as pd
from Bio import SeqIO
import os




for i in range(1,13):
    file="chr%s_AllcheckedHomo.fa"%(i)
    records=list(SeqIO.parse(file,"fasta"))
    for rec in records:
        id=str(rec.id)
        fam=id.split("_")[-1]
        for family in ["DTA","DTC","DTH","DTM","DTT"]:
            if fam==family:
                output=open("chr%s_%s_allCheckIRFHomo.fa"%(i,fam),"a+")
                output.write(">"+str(rec.id)+"\n"+str(rec.seq)+"\n")
