# Gut-Metagenome-Pipeline-Based-on-Nanopore-Sequencing


This pipeline was used for raw sequence quailty control, metagenome assembly, binning of metagenome-assembled genomes (MAGs), detection of structural variations and and dynamic methylomes based on the NGS and nanopore sequencing data.

This pipeline was described in the following publication:
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

ONT and Illumina sequencing data generated from this study is deposited at National Microbiology Data Center (https://nmdc.cn/) with accession number of NMDC10017965.


## SOFTWARE

MetaWRAP 1.2

Guppy v.3.3.0

minimap2 v2.17-r941

Canu 1.7

Flye 2.8.1-b1676

OPERA-MS v0.9.0

Spades v3.13.0

fastANI 1.32

MUM&Co v2.4.2

IGV 2.6.2

Prokka 1.13

emapper.py 1.0.3

seqkit v0.10.2

dRep v2.6.2

gtdbtk 1.3.0

RNammer -1.2

samtools 0.13.0

ProphageHunter

CD-hit v4.7

MAFFT v7.450

CAT

IQTREE2

iTOL v5

CRSPRDetect v2.4

blast v2.6.0

Tombo v1.5.1

Nanopolish v0.11.1

mCaller v0.3

Nanodisco v1.0.0




## Structure of this repository


0_Zymo_community_standard_data_assembly.sh

This pipeline was used for evaluation of assembly quality of Canu 1.7, Flye 2.8.1-b1676, OPERA-MS and Spades based on Zymo_community_Standards_data.

0_1_self_quant_bins.sh

One script was used for the estimation of abundance of MAGs based on the NGS reads.


1_Gut_metagenome_assembly&binning.sh

Part I read quality control, hybird assembly, binning of contigs, dreplication and classficatio of MAG.


2_Gut_metagenome_SV_analysis.sh

Part II Structure detection of MAGs and related genes prediction and annotation.


3_Prophage_prediction.sh

Part III. Prediction of prophage, redundancy removal, gene annotation, and construction of evolutionary trees.


4_CRISPR_spacers_prediction.sh

Part IV. Predicting CRISPR spacers from the contigs of all samples


5_Methylation.sh

Part V. Determination of DNA methylation using tombo software in genomes with more than 10 samples present.


Directory: metadata

This directory include the matadata of this experiment and some useful intermediate file.
