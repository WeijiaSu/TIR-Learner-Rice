# TIR-Learner-Rice

v1.11

# ** ** TIR-Learner: a new ensemble method for TIR Transposable Element annotation
## Introduction
TIR-Learner is a ensemble pipeline for Terminal Inverted Repeat (TIR) transposable elements annotation.  TIR-Learner combines homology-based detection using a newly curated TIR element library with a structure-based method to identify TIR elements de novo.  The structure-based component includes Machine Learning (ML) algorithms to classify the output sequences into five major TIR superfamilies.
## Modules
There are three modules in TIR-Learner.

**The first module**, called Homology Based, uses a library of TE reference sequences as query for BLAST with genomic sequence.  Only hits with 100% coverage were kept to generate a dataset of homologous TE candidates.  These candidate TEs were further checked to confirm the presence of both TIRs and TSDs, requiring at least 80% similarity as previously described.  The output of Module 1 is a set of high-quality intact TEs with significant homology to known TEs.  Importantly, partial TE copies  and TE internal fragments would not be included in this set.

**The second module** utilizes both sequence homology and structural features to identify TEs with partial but incomplete similarity to existing TEs.  First, genomic DNA was analyzed by software termed [GRF](https://github.com/weijiaweijia/GenericRepeatFinder) to detect sequences ranging in size from 50 bp to 10,000 bp that are flanked by inverted repeats (putative TIRs).  These candidate segments were then searched for homology using the reference TE library described in Module 1.  Sequences showing at least 80% similarity and 80% coverage of a reference TE were retained as putative TIR elements.

**The third module** (Structure-based de novo) processes the candidate sequences that were excluded from the second module because they did not have 80% similarity or coverage of a reference TE.  These sequences were analyzed by Machine Learning methods to identify sequence motifs conserved in each of the five major TIR superfamilies.  The Machine Learning component then classifies each candidate sequence into one of the five TIR superfamilies, or a nonTIR class.

![Picture1](https://user-images.githubusercontent.com/32049018/59524179-d3bc0a80-8e98-11e9-8136-97edc28e74bf.png)



## Usage
1. Install blast+ program following the instruction: https://www.ncbi.nlm.nih.gov/books/NBK279690/

2. Install GRF following the instruction: (https://github.com/weijiaweijia/GenericRepeatFinder)

3. Make sure the following python3 packages are installed:
(1)BioPython
(2)Pandas
(3)Sklearn (0.19.0)

4. Run TIR-Learner

4.1 Download and unzip TIR-Learner package.
Note: unzip all files in the folder including RefLib.zip in the current folder and Rice_model.sav.zip in the Module3 folder.

Use the script TIR-Learner.sh to run the program, change the paramters accordingly.

## Results
TIR-learner generates two output files.  One annotation file in gff3 format and one fasta file.  These two files are stored in the folder TIR-Learner-results.

## Example

TIR-Learner.sh

```shell

#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
#
#SBATCH --time=50:03:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mail-user=weijia@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --partition=freecomputing


#SBATCH --job-name='TIR-Learner-RiceAll'

#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module use /opt/rit/modules
module load python/3
module load cd-hit
module load gcc

######################################################################################################################################################
######################################################################################################################################################
######################################################## Change your parameters accordingly  #########################################################

# The absolute path to the genome file
genomeFile="/work/LAS/thomasp-lab/weijia/research/Rice/Genome/Rice.fa"                       ########################################################
# Name of the genome (prefix of intermediate files)
genomeName="Rice"                                                                           ##########################################################
# The absolute path to the TIR-Learner pipeline
path="/work/LAS/thomasp-lab/weijia/research/TIR-Learner1.11"                                ##########################################################
# The absolute path to your taget folder (current folder)
dir="/work/LAS/thomasp-lab/weijia/research/test_allRice"                                    ##########################################################
# Number of processor
t=16                                                                                        ##########################################################
# Path to the GRF program
grfp="/work/LAS/thomasp-lab/weijia/research/software/GenericRepeatFinder/bin"               ##########################################################
######################################################################################################################################################

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################





echo "############################################################ Pre-Processing ###########################################################"


python3 $path/pre.py -g $genomeFile -name $genomeName

for i in Module1 Module2 Module3; do mkdir $i; done

######################################################################################################################################################
############################################################################Module   1 ###############################################################
######################################################################################################################################################
echo "############################################################ Module   1  ###########################################################"

cd $dir"/Module1"
mkdir $genomeName
mkdir temp

echo "Module 1, Step 1: Blast Genome against Reference Library"

python3 $path/Module1/Blast_Ref.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module1"
cp $genomeName/*blast* temp/

echo "Module 1, Step 2: Select 100% coverage entries from Blast results"
python3 $path/Module1/Fullcov.py  -name $genomeName -p $path -d $dir"/Module1"

echo "Module1, Step 3: Making blastDB and get candidate sequences"

python3 $path/Module1/GetSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module1"

echo "Module1, Step 4: Check TIR and TSD"

python3 $path/Module1/CheckTIRTSD.py -name $genomeName -p $path -t $t -d $dir"/Module1"

echo "Module1, Step 5: Write to Gff3"

python3 $path/Module1/WriteToGff_M1.py -name $genomeName -p $path -t $t -d $dir"/Module1"

echo "Module1, Step 6: Check Low Complexity"

python3 $path/Module1/Lowcomp_M1.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module1"

echo "############################################################ Module   1 Finished ###########################################################"
echo "############################################################ Module   2  ###########################################################"


cd $dir"/Module2/"
mkdir $genomeName
mkdir temp

echo "Module 2, Step 1: Split Genome and  Run GRF program to find Inverted Repeats"
python3 $path/Module2/RunGRF.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2" -grfp $grfp
cp -r $genomeName/$genomeName* temp/

echo "Module 2, Step 2: Process GRF results"
python3 $path/Module2/ProcessGRFmite.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2"
cp -r $genomeName/*-p temp/

echo "Module 2 , Step 3 : GRF result blast reference sequences"
python3 $path/Module2/Blast_ref.py -name $genomeName -p $path -t $t -d $dir"/Module2"
cp $genome/*_RefLib temp/



echo "Module 2 , Step 4: Get sequences from 80% similarity"
python3 $path/Module2/GetSeq_M2.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2"

echo "Module 2 , Step 5: Check TIR and TSD"
python3 $path/Module2/CheckTIRTSD_M2.py -name $genomeName -p $path -t $t -d $dir"/Module2"

echo "Module 2 , Step 6: Write to gff3"
python3 $path/Module2/WriteToGff_M2.py -name $genomeName -p $path -t $t -d $dir"/Module2"

echo "Module 2 , Step 7: Remove Low complexity"
python3 $path/Module2/Lowcomp_M2.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2"

echo "############################################################ Module   2 Finished ###########################################################"
echo "############################################################ Module   3  ###########################################################"



cd $dir"/Module3/"

mkdir $genomeName
mv ../Module2/$genomeName/*_nonHomo.fa $genomeName/
mkdir temp

echo "Module 3, Step 1: Get dataset"
python3 $path/Module3/getDataset.py -name $genomeName -p $path -t $t -d $dir"/Module3"
cp $genomeName/*.csv temp/
cp $genomeName/*nonHomo.fa temp/

echo "Module 3, Step 2: ML prediction"
python3 $path/Module3/ML_Ensemble.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3"

echo "Module 3, Step 3: Check TIR/TSD"
python3 $path/Module3/CheckTIRTSD_M3.py -name $genomeName -p $path -t $t -d $dir"/Module3"

echo "Module 3, Step 4: Write to Gff"
python3 $path/Module3/WriteToGff_M3.py -name $genomeName -p $path -t $t -d $dir"/Module3"

echo "Module 3, Step 5: Remove Low complexity"
python3 $path/Module3/Lowcomp_M3.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3"

echo "############################################################ Module   3 Finished ###########################################################"

echo "############################################################ Post Processing  ###########################################################"

cd $dir

mkdir $genomeName


echo "Get Final GFF" 
python3 $path/CombineAll.py -name $genomeName -p $path -t $t -d $dir

mv *.gff3 $genomeName
rm *Low

echo "Get fasta file"
python3 $path/GetAllSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir


mkdir TIR-Learner-Result
mv $genomeName/*FinalAnn.gff3 TIR-Learner-Result/
mv $genomeName/*FinalAnn.fa TIR-Learner-Result/
rm -r $genomeName"_combine"
#rm -r $genomeName


echo "############################################################ TIR-Learner is finished! ###########################################################"



```

