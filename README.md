# TIR-Learner-Rice
# ** ** TIR-Learner: a new ensemble method for TIR Transposable Element annotation
## Introduction
TIR-Learner is a ensemble pipeline for Terminal Inverted Repeat (TIR) transposable elements annotation.  TIR-Learner combines homology-based detection using a newly curated TIR element library with a structure-based method to identify TIR elements de novo.  The structure-based component includes Machine Learning (ML) algorithms to classify the output sequences into five major TIR superfamilies.
## Modules
There are three modules in TIR-Learner.

The first module, called Homology Based, uses a library of TE reference sequences as query for BLAST with genomic sequence.  Only hits with 100% coverage were kept to generate a dataset of homologous TE candidates.  These candidate TEs were further checked to confirm the presence of both TIRs and TSDs, requiring at least 80% similarity as previously described.  The output of Module 1 is a set of high-quality intact TEs with significant homology to known TEs.  Importantly, partial TE copies  and TE internal fragments would not be included in this set.

The second module utilizes both sequence homology and structural features to identify TEs with partial but incomplete similarity to existing TEs.  First, genomic DNA was analyzed by software termed [GRF](https://github.com/weijiaweijia/GenericRepeatFinder) to detect sequences ranging in size from 50 bp to 10,000 bp that are flanked by inverted repeats (putative TIRs).  These candidate segments were then searched for homology using the reference TE library described in Module 1.  Sequences showing at least 80% similarity and 80% coverage of a reference TE were retained as putative TIR elements.

The third module (Structure-based de novo) processes the candidate sequences that were excluded from the second module because they did not have 80% similarity or coverage of a reference TE.  These sequences were analyzed by Machine Learning methods to identify sequence motifs conserved in each of the five major TIR superfamilies.  The Machine Learning component then classifies each candidate sequence into one of the five TIR superfamilies, or a nonTIR class.

![Picture1](https://user-images.githubusercontent.com/32049018/58972658-aa59fb00-8783-11e9-978b-24b8d68a24bb.png)



## How to use
1. Install GRF following the instruction: (https://github.com/weijiaweijia/GenericRepeatFinder)

2. Make sure the following python3 packages are installed:
(1)BioPython
(2)Pandas
(3)Sklearn

3 Run TIR-Learner

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
