## by Liang CHEN

# this pipeline was used for evaluation of assembly quality of Canu 1.7, Flye 2.8.1-b1676, OPERA-MS and Spades based on Zymo_community_Standards_data

## /data/chenliang/Zymo_Community_Standards_data_assembly

source activate assembly

## /data/chenliang/Zymo_Community_Standards_data_assembly
NGSdatadir="/data/chenliang/MC.Hiseq"
nanodatadir="/data/chenliang/Zymo_Community_Standards_data"
s=(Zymo-GridION-EVEN)
home=$(pwd)

## canu assembly
mkdir canu_assembly 
for i in ${s[@]};do
  time canu -p genome.60m -d canu_assembly/${i}_canu_assembly maxThreads=40 gnuplotTested=true genomeSize=60m \
   correctedErrorRate=0.035 -nanopore-raw  ${nanodatadir}/${i}-BB-SN.fq
  quast canu_assembly/${i}_canu_assembly/genome.60m.contigs.fasta -o canu_assembly/${i}_canu_assembly/quast_results/
done

# flye assembly
mkdir flye_assembly
for i in ${s[@]};do
    ## i(intreaction)次数太高会比较耗时,时间成倍增加，建议2-3次
 time flye --nano-raw ${nanodatadir}/${i}-BB-SN.fq --out-dir flye_assembly/${i}_flye_assembly -i 10 --meta --threads 40
 quast flye_assembly/${i}_flye_assembly/assembly.fasta -o flye_assembly/${i}_flye_assembly/quast_results/
done

##hybrid assembly using MetaSPAdes 
mkdir spades_hybird_assembly 
for i in ${s[@]};do
  time spades.py --meta -1 ${NGSdatadir}/MC_1.fastq -2 ${NGSdatadir}/MC_2.fastq --nanopore  ${nanodatadir}/${i}-BB-SN.fq \
   --threads 40 --memory 100 -o spades_hybird_assembly/${i}_spades_hybird_assembly
   quast spades_hybird_assembly/${i}_spades_hybird_assembly/scaffolds.fasta -o spades_hybird_assembly/${i}_spades_hybird_assembly/quast_results/
done

##only Illumina reads assembly using MetaSPAdes
mkdir spades_NGS_assembly 
for i in ${s[@]};do
  time spades.py --meta -1 ${NGSdatadir}/MC_1.fastq -2 ${NGSdatadir}/MC_2.fastq \
   --threads 40 --memory 100 -o spades_NGS_assembly/${i}_spades_NGS_assembly
   quast spades_NGS_assembly/${i}_spades_NGS_assembly/scaffolds.fasta -o spades_NGS_assembly/${i}_spades_NGS_assembly/quast_results/
done

##OPERA-MS assembly
conda create -n operams -c compbiocore perl-switch perl==5.26.2
conda activate operams
conda install -c bioconda perl-file-which perl-statistics-basic perl-statistics-r
git clone https://github.com/CSB5/OPERA-MS.git
cd OPERA-MS
make
perl OPERA-MS.pl check-dependency
perl OPERA-MS.pl install-db

source activate operams
cd Zymo_Community_Standards_data_assembly/Zymo-GridION-EVEN_operams_assembly
nohup time perl /software_users/chenliang/OPERA-MS/OPERA-MS.pl sample_config.config > OPERA-MS.nohup &
quast operams_output/contigs.fasta -o operams_output/quast_results/
cd $home



##calculate ANI for the assemnled contigs
cat /data/chenliang/ZymoBIOMICS.STD.refseq.v2/Genomes/*_complete_genome.fasta >/data/chenliang/ZymoBIOMICS.STD.refseq.v2/Genomes/8.bacteria.complete_genome.fasta

fastANI -r /data/chenliang/ZymoBIOMICS.STD.refseq.v2/Genomes/8.bacteria.complete_genome.fasta \
 -q flye_assembly/Zymo-GridION-EVEN_flye_assembly/assembly.fasta -o flye_assembly/Zymo-GridION-EVEN_flye_assembly/assembly.fasta.8.bacteria.ANI

fastANI -r /data/chenliang/ZymoBIOMICS.STD.refseq.v2/Genomes/8.bacteria.complete_genome.fasta \
 -q canu_assembly/Zymo-GridION-EVEN_canu_assembly/genome.60m.contigs.fasta -o canu_assembly/Zymo-GridION-EVEN_canu_assembly/genome.60m.contigs.fasta.8.bacteria.ANI

fastANI -r /data/chenliang/ZymoBIOMICS.STD.refseq.v2/Genomes/8.bacteria.complete_genome.fasta \
 -q spades_hybird_assembly/Zymo-GridION-EVEN_spades_hybird_assembly/contigs.fasta -o Zymo-GridION-EVEN-BB-SN_binning_metawrap_70_10_bins/all.contigs.8.bacteria.ANI

fastANI -r /data/chenliang/ZymoBIOMICS.STD.refseq.v2/Genomes/8.bacteria.complete_genome.fasta \
 -q spades_NGS_assembly/Zymo-GridION-EVEN_spades_NGS_assembly/contigs.fasta -o spades_NGS_assembly/Zymo-GridION-EVEN_spades_NGS_assembly/contigs.fasta.8.bacteria.ANI

fastANI -r /data/chenliang/ZymoBIOMICS.STD.refseq.v2/Genomes/8.bacteria.complete_genome.fasta \
 -q Zymo-GridION-EVEN_operams_assembly/operams_output/contigs.fasta -o Zymo-GridION-EVEN_operams_assembly/operams_output/contigs.fasta.8.bacteria.ANI



##binning
source activate metawrap

##for canu assembly
for i in ${s[@]};do
metawrap binning -o canu_assembly/${i}-BB-SN_binning -t 40 -a canu_assembly/${i}_canu_assembly/genome.60m.contigs.fasta \
  --metabat2 --maxbin2 --concoct --single-end ${nanodatadir}/${i}-BB-SN.fastq

metawrap bin_refinement -o canu_assembly/${i}-BB-SN_binning/bin_refinement -t 40 -m 100\
  -A canu_assembly/${i}-BB-SN_binning/metabat2_bins/ \
  -B canu_assembly/${i}-BB-SN_binning/maxbin2_bins/ \
  -C canu_assembly/${i}-BB-SN_binning/concoct_bins/ \
  -c 70 -x 10

done

##for flye assembly
for i in ${s[@]};do
metawrap binning -o flye_assembly/${i}-BB-SN_binning -t 40 -a flye_assembly/${i}_flye_assembly/scaffolds.fasta \
  --metabat2 --maxbin2 --concoct --single-end ${nanodatadir}/${i}-BB-SN.fastq

metawrap bin_refinement -o flye_assembly/${i}-BB-SN_binning/bin_refinement -t 40 -m 100\
  -A flye_assembly/${i}-BB-SN_binning/metabat2_bins/ \
  -B flye_assembly/${i}-BB-SN_binning/maxbin2_bins/ \
  -C flye_assembly/${i}-BB-SN_binning/concoct_bins/ \
  -c 70 -x 10

done

##for hybrid assembly using MetaSPAdes 
for i in ${s[@]};do
metawrap binning -o spades_hybird_assembly/${i}-BB-SN_binning -t 60 -a spades_hybird_assembly/${i}_spades_hybird_assembly/scaffolds.fasta \
  --metabat2 --maxbin2 --concoct --single-end ${nanodatadir}/${i}.NGS.nano.fastq

metawrap bin_refinement -o spades_hybird_assembly/${i}-BB-SN_binning/bin_refinement -t 60 \
  -A spades_hybird_assembly/${i}-BB-SN_binning/metabat2_bins/ \
  -B spades_hybird_assembly/${i}-BB-SN_binning/maxbin2_bins/ \
  -C spades_hybird_assembly/${i}-BB-SN_binning/concoct_bins/ \
  -c 70 -x 10

done

##for only Illumina reads assembly using MetaSPAdes
for i in ${s[@]};do
metawrap binning -o spades_NGS_assembly/${i}-BB-SN_binning -t 60 -a spades_NGS_assembly/${i}_spades_NGS_assembly/scaffolds.fasta \
  --metabat2 --maxbin2 --concoct --single-end ${nanodatadir}/${i}.NGS.nano.fastq

metawrap bin_refinement -o spades_NGS_assembly/${i}-BB-SN_binning/bin_refinement -t 60 \
  -A spades_NGS_assembly/${i}-BB-SN_binning/metabat2_bins/ \
  -B spades_NGS_assembly/${i}-BB-SN_binning/maxbin2_bins/ \
  -C spades_NGS_assembly/${i}-BB-SN_binning/concoct_bins/ \
  -c 70 -x 10

done

##for OPERA-MS assembly
for i in ${s[@]};do

metawrap binning -o Zymo-GridION-EVEN_operams_assembly/${i}-BB-SN_binning -t 40 -a Zymo-GridION-EVEN_operams_assembly/operams_output/contigs.fasta \
  --metabat2 --maxbin2 --concoct --single-end ${nanodatadir}/${i}-BB-SN.fastq

metawrap bin_refinement -o Zymo-GridION-EVEN_operams_assembly/${i}-BB-SN_binning/bin_refinement -t 40 -m 100\
  -A Zymo-GridION-EVEN_operams_assembly/${i}-BB-SN_binning/metabat2_bins/ \
  -B Zymo-GridION-EVEN_operams_assembly/${i}-BB-SN_binning/maxbin2_bins/ \
  -C Zymo-GridION-EVEN_operams_assembly/${i}-BB-SN_binning/concoct_bins/ \
  -c 70 -x 10

done


##calculation of code density

source activate metawrap
for i in $(ls canu_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/*fa); do
  sample=`echo ${i} |cut -d'/' -f5 | sed 's/.fa//'`
  mkdir canu_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/prokka/
  awk '{if(/>/){print substr($1,1,37)}else{print $0}}' ${i} >canu_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/prokka/${sample}.fa
  #prodigal -i hybird_temp/prodigal/${sample}_contig.fa -c -m -g 11 -p meta -f sco -q -a hybird_temp/prodigal/${sample} -d hybird_temp/prodigal/${sample} -o hybird_temp/prodigal/${sample}
  time prokka canu_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/prokka/${sample}.fa --outdir canu_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/prokka/${sample} \
  --prefix ${sample} --metagenome --kingdom Archaea,Bacteria,Mitochondria,Viruses \
  --force --cpus 60

  genome_len=`seqkit stat -T canu_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/prokka/${sample}.fa | sed '1d' |cut -f 5`
  code_len=`cut -f 3 canu_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/prokka/${sample}/${sample}.tsv |sed '1d' |awk '{sum+=$1} END {print sum}' `
  echo $code_len $genome_len >>canu_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/prokka/code_density.txt
done

for i in $(ls flye_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_5_bins/*fa); do
  sample=`echo ${i} |cut -d'/' -f5 | sed 's/.fa//'`
  mkdir flye_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_5_bins/prokka/

  awk '{if(/>/){print substr($1,1,37)}else{print $0}}' ${i} >flye_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_5_bins/prokka/${sample}.fa
  #prodigal -i hybird_temp/prodigal/${sample}_contig.fa -c -m -g 11 -p meta -f sco -q -a hybird_temp/prodigal/${sample} -d hybird_temp/prodigal/${sample} -o hybird_temp/prodigal/${sample}
  time prokka flye_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_5_bins/prokka/${sample}.fa --outdir flye_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_5_bins/prokka/${sample} \
  --prefix ${sample} --metagenome --kingdom Archaea,Bacteria,Mitochondria,Viruses \
  --force --cpus 60

  genome_len=`seqkit stat -T flye_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_5_bins/prokka/${sample}.fa | sed '1d' |cut -f 5`
  code_len=`cut -f 3 flye_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_5_bins/prokka/${sample}/${sample}.tsv |sed '1d' |awk '{sum+=$1} END {print sum}' `
  echo $code_len $genome_len >>flye_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_5_bins/prokka/code_density.txt
done

for i in $(ls spades_hybird_assembly/Zymo-GridION-EVEN-BB-SN_NGS_binning/bin_refinement/metawrap_70_10_bins/*fa); do
  sample=`echo ${i} |cut -d'/' -f5 | sed 's/.fa//'`
  mkdir spades_hybird_assembly/Zymo-GridION-EVEN-BB-SN_NGS_binning/bin_refinement/metawrap_70_10_bins/prokka/

  awk '{if(/>/){print substr($1,1,37)}else{print $0}}' ${i} >spades_hybird_assembly/Zymo-GridION-EVEN-BB-SN_NGS_binning/bin_refinement/metawrap_70_10_bins/prokka/${sample}.fa
  #prodigal -i hybird_temp/prodigal/${sample}_contig.fa -c -m -g 11 -p meta -f sco -q -a hybird_temp/prodigal/${sample} -d hybird_temp/prodigal/${sample} -o hybird_temp/prodigal/${sample}
  time prokka spades_hybird_assembly/Zymo-GridION-EVEN-BB-SN_NGS_binning/bin_refinement/metawrap_70_10_bins/prokka/${sample}.fa --outdir spades_hybird_assembly/Zymo-GridION-EVEN-BB-SN_NGS_binning/bin_refinement/metawrap_70_10_bins/prokka/${sample} \
  --prefix ${sample} --metagenome --kingdom Archaea,Bacteria,Mitochondria,Viruses \
  --force --cpus 60

  genome_len=`seqkit stat -T spades_hybird_assembly/Zymo-GridION-EVEN-BB-SN_NGS_binning/bin_refinement/metawrap_70_10_bins/prokka/${sample}.fa | sed '1d' |cut -f 5`
  code_len=`cut -f 3 spades_hybird_assembly/Zymo-GridION-EVEN-BB-SN_NGS_binning/bin_refinement/metawrap_70_10_bins/prokka/${sample}/${sample}.tsv |sed '1d' |awk '{sum+=$1} END {print sum}' `
  echo $code_len $genome_len >>spades_hybird_assembly/Zymo-GridION-EVEN-BB-SN_NGS_binning/bin_refinement/metawrap_70_10_bins/prokka/code_density.txt
done

for i in $(ls spades_NGS_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/*fa); do
  sample=`echo ${i} |cut -d'/' -f5 | sed 's/.fa//'`
  mkdir spades_NGS_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/prokka/

  awk '{if(/>/){print substr($1,1,37)}else{print $0}}' ${i} >spades_NGS_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/prokka/${sample}.fa
  #prodigal -i hybird_temp/prodigal/${sample}_contig.fa -c -m -g 11 -p meta -f sco -q -a hybird_temp/prodigal/${sample} -d hybird_temp/prodigal/${sample} -o hybird_temp/prodigal/${sample}
  time prokka spades_NGS_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/prokka/${sample}.fa --outdir spades_NGS_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/prokka/${sample} \
  --prefix ${sample} --metagenome --kingdom Archaea,Bacteria,Mitochondria,Viruses \
  --force --cpus 60

  genome_len=`seqkit stat -T spades_NGS_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/prokka/${sample}.fa | sed '1d' |cut -f 5`
  code_len=`cut -f 3 spades_NGS_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/prokka/${sample}/${sample}.tsv |sed '1d' |awk '{sum+=$1} END {print sum}' `
  echo $code_len $genome_len >>spades_NGS_assembly/Zymo-GridION-EVEN-BB-SN_binning_NGS/bin_refinement/metawrap_70_10_bins/prokka/code_density.txt
done



## Abundance estimation of bins from hybrid assembly using MetaSPAdes based on nanopore reads and NGS reads separately

# for nanopore reads
source activate nanosoft
cd Zymo_Community_Standards_data_assembly/spades_hybird_assembly/Zymo-GridION-EVEN-BB-SN_binning/bin_refinement/metawrap_70_10_bins
for i in ${s[@]};do
    for j in `ls *fa`;do
        minimap2 -a ${j} ${nanodatadir}/${i}-BB-SN.fastq >${j}.sam
        /software_users/chenliang/miniconda3/envs/metawrap/bin/samtools sort -T tmp ${j}.sam -o ${j}.sorted.bam
        /software/software/miniconda3/envs/assembly/lib/python2.7/site-packages/quast_libs/bedtools/bin/genomeCoverageBed -bg -ibam ${j}.sorted.bam > ${j}.bed.coverage
        awk -v j=${j} '{sum+=$4} END {print j " abu=", sum/NR}' ${j}.bed.coverage > ${j}_abu.txt
        rm ${j}.sam
        rm ${j}.sorted.bam
    done
done
cd $home

# for NGS reads
cd spades_hybird_assembly/Zymo-GridION-EVEN-BB-SN_binning 
bash self_quant_bins.sh
cd $home
