#########################################################################
# File Name: self_quant_bins.sh
# Author: Chen Liang
# mail: chen318liang@163.com
# Created Time: Mon 27 May 2019 04:46:57 PM CST
#########################################################################
#!/bin/bash

threads="20"
assembly="/data/chenliang/Zymo_Community_Standards_data_assembly/spades_hybird_assembly/Zymo-GridION-EVEN_spades_hybird_assembly/scaffolds.fasta"
out="quant_bins"
reads_1="/data/chenliang/MC.Hiseq/MC_1.fastq"
reads_2="/data/chenliang/MC.Hiseq/MC_2.fastq"
tmp=${reads_1##*/}
sample=${tmp%_*}
bin_folder="bin_refinement/metawrap_70_10_bins/"


salmon index -p $threads -t $assembly -i ${out}/assembly_index
salmon quant -i ${out}/assembly_index --libType IU -1 $reads_1 -2 $reads_2 -o ${out}/alignment_files/${sample}.quant --meta -p $threads

home=$(pwd)
cd ${out}/alignment_files/
/software_users/chenliang/miniconda3/envs/metawrap/bin/metawrap-scripts/summarize_salmon_files.py
cd $home
mkdir ${out}/quant_files
#for f in $(ls ${out}/alignment_files/ | grep .quant.counts); do mv ${out}/alignment_files/$f ${out}/quant_files/; done

#n=$(ls ${out}/quant_files/ | grep counts | wc -l)
/software_users/chenliang/miniconda3/envs/metawrap/bin/metawrap-scripts/split_salmon_out_into_bins.py ${out}/quant_files/ $bin_folder $assembly > ${out}/bin_abundance_table.tab
