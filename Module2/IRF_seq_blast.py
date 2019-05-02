import pandas as pd
from Bio import SeqIO
import os
import multiprocessing
from multiprocessing import Pool

def TEDB_Chr(chr):
    tedb = "nonStandard_all.fa"
    seq = "IRF_chr%s.fa"%(chr)
    outname = str(chr) + "_" + "IRF" + "_" + "withhomology"
    blast = "blastn -query %s -subject %s -outfmt '7 qseqid sseqid length pident gaps mismatch qstart qend sstart send evalue qcovhsp' -out %s" % (seq, tedb, outname)
    os.system(blast)


if __name__ == '__main__':
     l=list(range(1,13))
     pool = multiprocessing.Pool(12)
     pool.map(TEDB_Chr,l)
     pool.close()
     pool.join()
