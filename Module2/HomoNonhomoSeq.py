import pandas as pd
from Bio import SeqIO
import os

#for i in range(1,13):
#    cat = "cat %s_GRF_DT*_withhomology_80_unique > %s_GRF_withhomology_80_unique"%(i,i)
#    os.system(cat)

def GetSeq(IRF_file,GRF_file,IRF_Homo_list,GRF_Homo_list,chr):
    f=pd.read_table(IRF_Homo_list,header=None,sep="\t")
    homolist=list(f[0])
    records=list(SeqIO.parse(IRF_file, "fasta"))
    SeqIO.write((seq for seq in records if str(seq.id) in homolist), "chr%s_HomoSeq.fa" % (chr), "fasta")
    SeqIO.write((seq for seq in records if str(seq.id) not in homolist), "chr%s_NonHomoSeq.fa" % (chr), "fasta")
    
   # f=pd.read_table(GRF_Homo_list,header=None,sep="\t")
   # homolist=list(f[0])
   # records=list(SeqIO.parse(GRF_file, "fasta"))
   # SeqIO.write((seq for seq in records if str(seq.id) in homolist), "chr%s_HomoSeqGRF.fa" % (chr), "fasta")
   # SeqIO.write((seq for seq in records if str(seq.id) not in homolist), "chr%s_NonHomoSeqGRF.fa" % (chr), "fasta")
   # cat = "chr%s_HomoSeqIRF.fa chr%s_HomoSeqGRF.fa > chr%s_HomoSeq.fa"%(chr,chr,chr)
   # os.system(cat)
   # cat = "chr%s_NonHomoSeqIRF.fa chr%s_NonHomoSeqIRF.fa > chr%s_NonHomoSeq.fa"%(chr,chr,chr)
   # os.system(cat)

for i in range(1,13):
    GetSeq("IRF_chr%s.fa"%(i),"GRFmite_Chr%s.fa"%(i),"%s_IRF_withhomology_80_unique"%(i),"%s_GRF_withhomology_80_unique"%(i),i)

