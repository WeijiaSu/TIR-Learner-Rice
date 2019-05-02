import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import subprocess
import multiprocessing
from multiprocessing import Pool

def getLTR(rec):
    s=str(rec.seq)
    ID=str(rec.id)
    seq1=s[0:200]
    seq2=s[-200:]
    seq1 = SeqRecord(Seq(seq1),id="seq1"+ID)
    seq2 = SeqRecord(Seq(seq2),id="seq2"+ID)
    SeqIO.write(seq1, "seq1_"+ID+".fasta", "fasta")
    SeqIO.write(seq2, "seq2_"+ID+".fasta", "fasta")
    blast = "blastn -query %s -subject %s -outfmt '7 qseqid sseqid length pident gaps mismatch qstart qend sstart send evalue qcovhsp' -out %s" % ("seq1_"+ID+".fasta", "seq2_"+ID+".fasta", ID+"_blast.tsv")
    os.system(blast)
    rm ="rm %s"%("seq2_"+ID+".fasta")
    os.system(rm)
    rm="rm %s"%("seq1_"+ID+".fasta")
    os.system(rm)
    

if __name__ == '__main__':
     records = list(SeqIO.parse("Rice_SB200F.fa","fasta"))
     pool = multiprocessing.Pool(16)
     pool.map(getLTR,records)
     pool.close()
     pool.join()

cat="cat *.tsv > Rice_IRFML_flankBlast"
os.system(cat)
rm ="rm *.tsv"
os.system(rm)
