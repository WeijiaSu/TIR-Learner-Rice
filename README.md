# TIR-Learner-Rice
# ** ** TIR-Learner: a new ensemble method for TIR Transposable Element annotation
## Introduction
TIR-Learner is a ensemble pipeline for Terminal Inverted Repeat (TIR) transposable elements annotation.  TIR-Learner combines homology-based detection using a newly curated TIR element library with a structure-based method to identify TIR elements de novo.  The structure-based component includes Machine Learning (ML) algorithms to classify the output sequences into five major TIR superfamilies.
## Modules
There are three modules in TIR-Learner.

**The first module**, called Homology Based, uses a library of TE reference sequences as query for BLAST with genomic sequence.  Only hits with 100% coverage were kept to generate a dataset of homologous TE candidates.  These candidate TEs were further checked to confirm the presence of both TIRs and TSDs, requiring at least 80% similarity as previously described.  The output of Module 1 is a set of high-quality intact TEs with significant homology to known TEs.  Importantly, partial TE copies  and TE internal fragments would not be included in this set.

**The second module** utilizes both sequence homology and structural features to identify TEs with partial but incomplete similarity to existing TEs.  First, genomic DNA was analyzed by software termed [GRF](https://github.com/weijiaweijia/GenericRepeatFinder) to detect sequences ranging in size from 50 bp to 10,000 bp that are flanked by inverted repeats (putative TIRs).  These candidate segments were then searched for homology using the reference TE library described in Module 1.  Sequences showing at least 80% similarity and 80% coverage of a reference TE were retained as putative TIR elements.

**The third module** (Structure-based de novo) processes the candidate sequences that were excluded from the second module because they did not have 80% similarity or coverage of a reference TE.  These sequences were analyzed by Machine Learning methods to identify sequence motifs conserved in each of the five major TIR superfamilies.  The Machine Learning component then classifies each candidate sequence into one of the five TIR superfamilies, or a nonTIR class.

![Picture1](https://user-images.githubusercontent.com/32049018/58972658-aa59fb00-8783-11e9-978b-24b8d68a24bb.png)



## Usage

1. Install GRF following the instruction: (https://github.com/weijiaweijia/GenericRepeatFinder)

2. Make sure the following python3 packages are installed:
(1)BioPython
(2)Pandas
(3)Sklearn

3. Run TIR-Learner

3.1 Download and unzip TIR-Learner package.

3.2 Copy preProcess.sh to your current folder and change the parameters accordingly, and run preProcess.sh

This step will pre-screen the genome sequence file to check if there are special characters in the sequence headers.  This script will also configure the program, such as creat seperate folder for each module and copy the masker script that is needed for each module.

3.3 Go to folder Module1 in your target folder.  Run Module1.sh (change parameters accordingly)

3.4 Go to folder Module2 in your target folder.  Run Module2.sh (change parameters accordingly)

3.5 Go to folder Module3 in your target folder.  Run Module3.sh (change parameters accordingly)
Note: Module1 is independent to Module2 and Module3, therefore, you can run Module1 and Module2 in the same time.  However, Module3 should begin after Module2 is finished.

3.6 Go to the target folder and run postProcessing.sh 
This script will combine and filter the results from the three modules, this will also remove and clean the temperary files that are generated during each module.

## Results
TIR-learner generates two output files.  One annotation file in gff3 format and one fasta file.  These two files are stored in the folder TIR-Learner-results.

## Example

Module1.sh

```shell

#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
#
#SBATCH --time=02:03:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mail-user=weijia@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --partition=freecomputing


#SBATCH --job-name='Module1'

#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module use /opt/rit/modules
export PATH=$PATH:/work/LAS/thomasp-lab/weijia/research/software/ncbi-blast-2.6.0+/bin

module load python/3

genomeFile="/work/LAS/thomasp-lab/weijia/research/test1.9/Genome/LNNJ01.1.fsa_nt"
genomeName="LNNJ"
path="/work/LAS/thomasp-lab/weijia/research/TIR-Learner1.9"
t=16
dir="/work/LAS/thomasp-lab/weijia/research/test1.9/Module1"

mkdir $genomeName
mkdir temp
echo "Module 1, Step 1: Blast Genome against Reference Library"

python3 $path/Module1/Blast_Ref.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir
cp $genomeName/*blast* temp/

echo "Module 1, Step 2: Select 100% coverage entries from Blast results"
python3 $path/Module1/Fullcov.py  -name $genomeName -p $path -d $dir

echo "Module1, Step 3: Making blastDB and get candidate sequences"

python3 $path/Module1/GetSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir

echo "Module1, Step 4: Check TIR and TSD"

python3 $path/Module1/CheckTIRTSD.py -name $genomeName -p $path -t $t -d $dir

echo "Module1, Step 5: Write to Gff3"

python3 $path/Module1/WriteToGff_M1.py -name $genomeName -p $path -t $t -d $dir

echo "Module1, Step 6: Check Low Complexity"

python3 $path/Module1/Lowcompl.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir
```
Module2.sh

```shell
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
#
#SBATCH --time=03:03:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mail-user=weijia@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --partition=freecomputing


#SBATCH --job-name='Module2'

#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module use /opt/rit/modules
export PATH=$PATH:/work/LAS/thomasp-lab/weijia/research/software/ncbi-blast-2.6.0+/bin

module load python/3
module load cd-hit
module load gcc


genomeFile="/work/LAS/thomasp-lab/weijia/research/test1.9/Genome/LNNJ01.1.fsa_nt"
genomeName="LNNJ"
path="/work/LAS/thomasp-lab/weijia/research/TIR-Learner1.9"
t=16
dir="/work/LAS/thomasp-lab/weijia/research/test1.9/Module2"
grfp="/work/LAS/thomasp-lab/weijia/research/software/GenericRepeatFinder/bin"


mkdir $genomeName
mkdir temp

echo "Module 2, Step 1: Split Genome and  Run GRF program to find Inverted Repeats"
#python3 $path/Module2/RunGRF.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir -grfp $grfp
cp -r $genomeName/$genomeName* temp/

echo "Module 2, Step 2: Process GRF results"
python3 $path/Module2/ProcessGRFmite.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir
cp -r $genomeName/*-p temp/

echo "Module 2 , Step 3 : GRF result blast reference sequences"
python3 $path/Module2/Blast_ref.py -name $genomeName -p $path -t $t -d $dir
cp $genome/*_RefLib temp/



echo "Module 2 , Step 4: Get sequences from 80% similarity"
python3 $path/Module2/GetSeq_M2.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir

echo "Module 2 , Step 5: Check TIR and TSD"
python3 $path/Module2/CheckTIRTSD_M2.py -name $genomeName -p $path -t $t -d $dir

echo "Module 2 , Step 6: Write to gff3"
python3 $path/Module2/WriteToGff_M2.py -name $genomeName -p $path -t $t -d $dir

echo "Module 2 , Step 7: Remove Low complexity"
python3 $path/Module2/Lowcomp_M2.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir
```
Module3.sh
```shell
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
#
#SBATCH --time=02:03:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mail-user=weijia@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --partition=freecomputing


#SBATCH --job-name='Module3'

#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module use /opt/rit/modules
export PATH=$PATH:/work/LAS/thomasp-lab/weijia/research/software/ncbi-blast-2.6.0+/bin

module load python/3

genomeFile="/work/LAS/thomasp-lab/weijia/research/test1.9/Genome/LNNJ01.1.fsa_nt"
genomeName="LNNJ"
path="/work/LAS/thomasp-lab/weijia/research/TIR-Learner1.9"
t=16
dir="/work/LAS/thomasp-lab/weijia/research/test1.9/Module3"


mkdir $genomeName
mv ../Module2/$genomeName/*_nonHomo.fa $genomeName/
mkdir temp

echo "Module 3, Step 1: Get dataset"
python3 $path/Module3/getDataset.py -name $genomeName -p $path -t $t -d $dir
cp $genomeName/*.csv temp/
cp $genomeName/*nonHomo.fa temp/

echo "Module 3, Step 2: ML prediction"
python3 $path/Module3/ML_Ensemble.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir

echo "Module 3, Step 3: Check TIR/TSD"
python3 $path/Module3/CheckTIRTSD_M3.py -name $genomeName -p $path -t $t -d $dir

echo "Module 3, Step 4: Write to Gff"
python3 $path/Module3/WriteToGff_M3.py -name $genomeName -p $path -t $t -d $dir

echo "Module 3, Step 5: Remove Low complexity"
python3 $path/Module3/Lowcomp_M3.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir

```

