import pandas as pd
from Bio import SeqIO
import os


d={}
d["DTA"]="3"
d["DTC"]="4"
d["DTH"]="5"
d["DTM"]="2"
d["DTT"]="1"

def writeTogff(fastafile,Family,output,source):
    record=list(SeqIO.parse(fastafile,"fasta"))
    out=open(output,"a+")
    for rec in record:
        ID=str(rec.id).split("_")

        if (source=="SB"):
            seq_ID = ID[0]
            p1 = int(ID[1])
            p2 = int(ID[2])
            type = Family
            strand = ""
            if (p1 < p2):
                start = p1
                end = p2
                strand = "."
            else:
                start = p2
                end = p1
                strand = "."
            length = end - start + 1
            if length < 50:
                pass
            attribut = ID[-6] + "_" + ID[-5] + "_" + ID[-4] + "_" + ID[-3] + "_" + ID[-2] + "_" + ID[-1] + "_" + str(
                length)
            out.write(seq_ID + "\t" + source + "\t" + type + "\t" + str(start) + "\t" + str(
                end) + "\t" + "." + "\t" + strand + "\t" + "." + "\t" + attribut + "\n")
        else:
            seq_ID=ID[-11]
            p1=int(ID[-10])
            p2=int(ID[-9])
            type = Family
            strand = ""
            if (p1 < p2 ):
                start = p1
                end = p2
                strand = "+"
            else:
                start = p2
                end = p1
                strand = "-"
            length = end - start + 1
            if length < 50:
                pass
            attribut = ID[-6] + "_" + ID[-5] + "_" + ID[-4] + "_" + ID[-3] + "_" + ID[-2] + "_" + ID[-1] + "_" + str(
                length)
            out.write(seq_ID + "\t" + source + "\t" + type + "\t" + str(start) + "\t" + str(
                end) + "\t" + "." + "\t" + strand + "\t" + "." + "\t" + attribut + "\n")

    out.close()
#
for i in ["DTA", "DTC", "DTH", "DTM", "DTT"]:
     for j in range(1,13):
          Family = i
          fastafile="chr%s_irfml_%s_withTIRTSD_Seq.fa" % (j,i)
          if (os.path.isfile(fastafile) != True):
              continue
          output = "chr%s_%s_IRFML.gff3" % (j,i)
          writeTogff(fastafile, Family, output,"SB")

cat ="cat chr*_IRFML.gff3 > Rice_IRFML.gff3"
os.system(cat)

