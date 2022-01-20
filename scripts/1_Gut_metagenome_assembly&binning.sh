##############################################################################################################################
## Part I read quality control, hybird assembly, binning of contigs, dreplication and classficatio of MAGs
##############################################################################################################################

home=$(pwd)

## 1. remove host contamination for nanopore reads

source activate nanosoft
mkdir -p result/dehost
mkdir -p nano_dehost
sed '1d' sample_design.txt |while read line; do
    sample=`echo ${line} |cut -d' ' -f1`
    minimap2 -ax map-ont -t 20 /software_users/chenliang/database/hg38/hg38.fa nano_cleandata/${sample}.qcat.fastq >nano_dehost/${sample}.sam
    awk '($2==4) {print $1}' nano_dehost/${sample}.sam >nano_dehost/${sample}.unmapped.names
    cp nano_dehost/${sample}.unmapped.names result/dehost
    seqkit grep -f nano_dehost/${sample}.unmapped.names nano_cleandata/${sample}.qcat.fastq > nano_dehost/nano.${sample}.dehost.fq
    gzip nano_dehost/nano.${sample}.dehost.fq
done


## 2. QC of NGS reads, hybird assembly using metaspades and binning of contigs

source activate metawrap
mkdir NGS_rawdata
mkdir read_qc
mkdir hybird_temp
mkdir hybird_temp/binning/

sed '1d' sample_design.txt |while read line; do
    sample=`echo ${line} |cut -d' ' -f1`

    gunzip -c Cleandata/META20IMWJ*_${sample}/META20IMWJ*_${sample}_R1.fq.gz >NGS_rawdata/${sample}_1.fastq
    gunzip -c Cleandata/META20IMWJ*_${sample}/META20IMWJ*_${sample}_R2.fq.gz >NGS_rawdata/${sample}_2.fastq

    # redas QC
    metawrap read_qc -t 60  --skip-pre-qc-report --skip-post-qc-report -1 NGS_rawdata/${sample}_1.fastq -2 NGS_rawdata/${sample}_2.fastq -o read_qc/${sample}

    # hybird pre-assembly
    spades.py --meta -1 read_qc/${sample}/final_pure_reads_1.fastq -2 read_qc/${sample}/final_pure_reads_2.fastq \
     --nanopore nano_dehost/nano.${sample}.dehost.fq.gz --threads 12 --memory 800 -o hybird_temp/spades/${sample}_assembly

    metawrap binning -o hybird_temp/binning/${sample}_binning -t 12 -a hybird_temp/spades/${sample}_assembly/contigs.fasta \
      --metabat2 --maxbin2 --concoct read_qc/${sample}/final_pure_reads_1.fastq read_qc/${sample}/final_pure_reads_2.fastq

    metawrap bin_refinement -o hybird_temp/binning/${sample}_binning/bin_refinement -t 12 \
      -A hybird_temp/binning/${sample}_binning/metabat2_bins/ \
      -B hybird_temp/binning/${sample}_binning/maxbin2_bins/ \
      -C hybird_temp/binning/${sample}_binning/concoct_bins/ \
      -c 70 -x 10

    # Re-assemble the consolidated bin set
    metawrap reassemble_bins -o hybird_temp/binning/${sample}_binning/bin_reassembly -1 read_qc/${sample}/final_pure_reads_1.fastq \
    -2 read_qc/${sample}/final_pure_reads_2.fastq -t 12 -m 800 -c 70 -x 10 \
    -b hybird_temp/binning/${sample}_binning/bin_refinement/metawrap_70_10_bins

    rm hybird_temp/binning/${sample}_binning/work_files/final_pure_reads.sam
    gzip read_qc/${sample}/*fastq

done


## 3. get the bins information of reassembled bins

mkdir result/metawrap_bins/hybird
sed '1d' sample_design.txt |while read line; do
  sample=`echo ${line} |cut -d' ' -f1`
  sed '1d' hybird_temp/binning/${sample}_binning/bin_reassembly/reassembled_bins.stats |while read line1; do
    fa=`echo ${line1} |cut -d' ' -f1`
    completeness=`echo ${line1} |cut -d' ' -f2`
    contamination=`echo ${line1} |cut -d' ' -f3`
    GC=`echo ${line1} |cut -d' ' -f4`
    lineage=`echo ${line1} |cut -d' ' -f5`
    N50=`echo ${line1} |cut -d' ' -f6`
    size=`echo ${line1} |cut -d' ' -f7`
    echo ${sample}_${completeness}_${contamination}_${fa}.fa,${completeness},${contamination},${GC},${lineage},${N50},${size} >>result/metawrap_bins/hybird/hybird_genomeInfo.csv
    cp hybird_temp/binning/${sample}_binning/bin_reassembly/reassembled_bins/${fa}.fa result/metawrap_bins/hybird/${sample}_${completeness}_${contamination}_${fa}.fa
  done
done
sed -i '1 i genome,completeness,contamination,GC,lineage,N50,size' result/metawrap_bins/hybird/hybird_genomeInfo.csv


## 4. dereplication of all reassembled bins

source activate drep
dRep dereplicate hybird_temp/drep/hybird -g result/metawrap_bins/hybird/*fa -pa 0.9 -sa 0.95 -nc 0.30 -cm larger -comp 70 -con 10 -p 80 --genomeInfo result/metawrap_bins/hybird/hybird_genomeInfo.csv
mkdir -p result/drep_bins/hybird
cp -r hybird_temp/drep/hybird/dereplicated_genomes/ result/drep_bins/hybird


## 5. classification of dereplicated MAGs

source activate gtdbtk
gtdbtk classify_wf --genome_dir result/drep_bins/hybird/dereplicated_genomes --out_dir hybird_temp/classify --cpus 40 -x fa --prefix hybird
gtdbtk infer --msa_file hybird_temp/classify/hybird.bac120.user_msa.fasta --out_dir hybird_temp/classify/bac120_infer_out --cpu 40
gtdbtk infer --msa_file hybird_temp/classify/hybird.ar122.user_msa.fasta --out_dir hybird_temp/classify/ar122_infer_out --cpu 40

source activate gtdbtk
gtdbtk de_novo_wf --genome_dir result/drep_bins/hybird/dereplicated_genomes \
  --out_dir hybird_temp/de_novo_classify \
  --extension fa \
  --bacteria \
  --outgroup_taxon p__Patescibacteria \
  --prefix de_novo \
  --cpus 40


## 6. prediction of rRNA genes

mkdir result/rnammer
mkdir result/rnammer/hybird
for i in $(ls result/drep_bins/hybird/dereplicated_genomes/*fa);do
  fa=`echo ${i} | cut -d'/' -f 5 |sed 's/.fa//'`
  mkdir result/rnammer/hybird/${fa}
  rnammer -S bac -multi -f result/rnammer/hybird/${fa}/rRNA.fasta -h result/rnammer/hybird/${fa}/rRNA.hmmreport -xml result/rnammer/hybird/${fa}/rRNA.xml \
  -gff result/rnammer/hybird/${fa}/rRNA.gff2 result/drep_bins/hybird/dereplicated_genomes/${fa}.fa
  cat result/rnammer/hybird/${fa}/rRNA.gff2 | tail -n +7 |sed '$d' | awk -v var="${fa}" '{print var "\t" $0}' >>result/rnammer/hybird/all_hybird_rnammer_result.txt
done


## 7. comparision with previous published Unified Human Gastrointestinal Genome (UHGG) collection

mkdir mgnify_genomes/ && cd mgnify_genomes/

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' genomes-nr_metadata.tsv representative_species.txt >representative_species.info

mkdir representative_species
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$24} NR>FNR{print "wget "a[$1]" -P representative_species"}' genomes-nr_metadata.tsv representative_species.txt >download.sh
bash download.sh
gunzip representative_species/*gz

# get fasta files of genome from gff files
for i in $(ls representative_species/*gff);do
  sample=`echo $i |sed 's/.gff//'`
#  line=$(awk '/##FASTA/ {print NR}' $i)
  awk -v line=$(awk '/##FASTA/ {print NR}' $i) '{if (NR>line) {print $0}}' $i | cat >$sample.fa
done 

cp ../result/drep_bins/hybird/dereplicated_genomes/*fa representative_species/

awk 'BEGIN{FS=OFS=","} NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../result/metawrap_bins/hybird/hybird_genomeInfo.csv ../hybird_temp/drep/hybird/data_tables/Widb.csv >genomeInfo.csv
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$9","$10","$8","$19","$7","$5} NR>FNR{print $1".fa,"a[$1]}' genomes-nr_metadata.tsv representative_species.txt >>genomeInfo.csv

## dereplication
source activate drep
dRep dereplicate drep -g representative_species/*fa -pa 0.9 -sa 0.95 -nc 0.30 -cm larger -l 50000 -comp 50 -con 25 -p 40 --genomeInfo genomeInfo.csv

cd $home


