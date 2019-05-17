#!/bin/bash

cd Module1/
# change the name of genome $genome--> YourGenomeName
for i in A C H M T;do blastn -query "Rice_DT"$i"_withTIR.fa" -subject ../$genome -outfmt '7 qseqid sseqid length pident gaps mismatch qstart qend sstart send evalue qcovhsp' -out "Blast_DT"$i; done

python3 Fullcov.py

for i in A C H M T; do python3 "extractTE_DT"$i".py" ;done

for i in A C H M T; do python3 "CheckTIRTSD_DT"$i"_best.py" ;done

python3 WriteToGff.py

cd ..
cd Module2/

python3 ProcessIRF.py

for i in {1..12}; do python3 "extractSeq_"$i".py" ;done

RemoveLowComplex.py

python3 IRF_seq_blast.py

python3 ProcessBlast.py

python3 HomoNonhomoSeq.py

python3 HomoSW_TSD.py

for i in {A C H M T}; do python3 "IRF_checkTIRTSD_DT"$i"_best.py"; done

python3 WriteTogff_IRFHB.py

cd ..
cd Module3/
for i in {1..12}; do python3 "IRF_dataset_chr"$i".py"; done 
python3 GetNonHomoPreSet.py

python3 ML_Ensemble.py

python3 prediction.py

python3 NonHomopre_SWTSD.py

for i in {A C H M T}; do python3 "IRF_checkTIRTSD_DT"$i"_best.py"; done

python3 WriteTogff_IRFHB.py

cd ..

for i in {1..3}; do mkdir "result_module_"$i; done

cp Module1/Rice_RefHB.gff3 result_module_1/result_module_1.gff3

cp Module2/Rice_IRFHB.gff3 result_module_2/result_module_2.gff3

cp Module3/Rice_IRFML.gff3 result_module_3/result_module_3.gff3

mkdir combine

python3 FinalOutput.py

cp combine/FinalAnn.gff3 .
cp combine/FinalAnn_noInternal.gff3 .
mv combine tem

rm Module1/*.fa
rm Module1/*.csv
rm Module2/*.fa
rm Module2/*.csv
rm Module3/*.fa
rm Module3/*.csv


