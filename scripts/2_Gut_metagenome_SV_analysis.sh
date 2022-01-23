##############################################################################################################################
## Part II Structure detection of MAGs and related genes prediction and annotation
##############################################################################################################################

home=$(pwd)

## 1. MUMandCo call SV 

conda activate mummer4
mkdir hybird_temp/drep/hybird/MUMandCo
OUTPUT_DIR="hybird_temp/drep/hybird/MUMandCo"

# filter the clusters present in >10 samples 
for i in `sed '1d' hybird_temp/drep/hybird/data_tables/Cdb.csv |cut -f2 -d',' |sort |uniq -c | sort -nr | awk '($1>10) {print}' |awk '{print $2}'`;do
  mkdir $OUTPUT_DIR/${i}
  for j in `sed '1d' hybird_temp/drep/hybird/data_tables/Cdb.csv | awk -v var="${i}" 'BEGIN {FS=","} ($2==var) {print $1}'`;do
    cp result/metawrap_bins/hybird/$j $OUTPUT_DIR/${i}
  done
  ref=`sed '1d' hybird_temp/drep/hybird/data_tables/Widb.csv |awk -v var="${i}" 'BEGIN {FS=","} ($8==var) {print $1}'`
  rm $OUTPUT_DIR/${i}/${ref}
  for m in $(ls $OUTPUT_DIR/${i}/*fa);do 
    m=`echo ${m} |cut -d'/' -f 6`
    gsize=`seqkit stat -T result/metawrap_bins/hybird/${ref} | sed '1d' |cut -f 5`
    cd $OUTPUT_DIR/${i}/
    # revise the code of mumandco_v2.4.2.sh (line 171) awk '{if($5>1000) print $0}' 
    bash /software_users/chenliang/MUMandCo/mumandco_v2.4.2.sh -r /data/nano_meta_ref/result/metawrap_bins/hybird/${ref} \
    -q /data/nano_meta_ref/$OUTPUT_DIR/${i}/${m} -g ${gsize} -o ${ref}-${m}

    cd /data/nano_meta_ref/
  done
  rm $OUTPUT_DIR/${i}/*fa
done


## 2.  test of SV accuracy 

mkdir igv_test

source activate nanosoft
ref="result/metawrap_bins/hybird/TD33_98.94_0.0_bin.20.orig.fa"
minimap2 -ax splice  -t 48 ${ref} temp/minimap2/nano.CD100.dehost.fq.gz > igv_test/nano.CD100vsTD33.dehost.aln.sam
samtools sort igv_test/nano.CD100vsTD33.dehost.aln.sam -o igv_test/nano.CD100vsTD33.dehost.aln.sorted.bam -@ 40
samtools index -@ 40 igv_test/nano.CD100vsTD33.dehost.aln.sorted.bam igv_test/nano.CD100vsTD33.dehost.aln.sorted.bam.bai

#grep "deletion" igv_test/DEL100_test_output/DEL100_test.SVs_all.tsv | cut -f2,7,8 |awk '{if($2<1000) {print $1"\t"$2"\t"$3+1000"\t"$3-$2} else {print $1"\t"$2-1000"\t"$3+1000"\t"$3-$2}}' |awk '($4>500) {print}' >test1/DEL100_test_output/deletion.ref.bed
cat hybird_temp/drep/hybird/MUMandCo/522_1/TD33_98.94_0.0_bin.20.orig.fa-CD100_92.19_2.004_bin.11.strict.fa_output/ref.deletion.bed | cut -f1-3 \
|awk '{if($2<1000) {print $1"\t"$2"\t"$3+1000"\t"$3-$2} else {print $1"\t"$2-1000"\t"$3+1000"\t"$3-$2}}' | cut -f1-3 >igv_test/ref_deletion.bed
conda activate assembly
make_IGV_snapshots.py  -ht 1000 -o igv_test/MUMandCo_ref_deletion_igv_snap -bin /software_users/chenliang/IGV-snapshot-automator/bin/IGV_2.3.81/igv.jar \
-g result/drep_bins/hybird/dereplicated_genomes/TD33_98.94_0.0_bin.20.orig.fa -r igv_test/ref_deletion.bed igv_test/nano.CD100vsTD33.dehost.aln.sorted.bam


minimap2 -ax splice  -t 48 ${ref} temp/minimap2/nano.TD33.dehost.fq.gz > igv_test/nano.TD33vsTD33.dehost.aln.sam
samtools sort igv_test/nano.TD33vsTD33.dehost.aln.sam -o igv_test/nano.TD33vsTD33.dehost.aln.sorted.bam -@ 40
samtools index -@ 40 igv_test/nano.TD33vsTD33.dehost.aln.sorted.bam igv_test/nano.TD33vsTD33.dehost.aln.sorted.bam.bai

conda activate assembly
make_IGV_snapshots.py  -ht 1000 -o igv_test/MUMandCo_ref_deletion_igv_snap_2 -bin /software_users/chenliang/IGV-snapshot-automator/bin/IGV_2.3.81/igv.jar \
-g result/drep_bins/hybird/dereplicated_genomes/TD33_98.94_0.0_bin.20.orig.fa -r igv_test/ref_deletion.bed igv_test/nano.CD100vsTD33.dehost.aln.sorted.bam \
igv_test/ref_deletion.bed igv_test/nano.TD33vsTD33.dehost.aln.sorted.bam



source activate nanosoft
ref="result/metawrap_bins/hybird/CD100_92.19_2.004_bin.11.strict.fa"
minimap2 -ax splice  -t 48 ${ref} temp/minimap2/nano.TD33.dehost.fq.gz > igv_test/nano.TD33vsCD100.dehost.aln.sam
samtools sort igv_test/nano.TD33vsCD100.dehost.aln.sam -o igv_test/nano.TD33vsCD100.dehost.aln.sorted.bam -@ 40
samtools index -@ 40 igv_test/nano.TD33vsCD100.dehost.aln.sorted.bam igv_test/nano.TD33vsCD100.dehost.aln.sorted.bam.bai

#grep "insertion" igv_test/DEL100_test_output/DEL100_test.SVs_all.tsv | cut -f2,7,8 |awk '{if($2<1000) {print $1"\t"$2"\t"$3+1000"\t"$3-$2} else {print $1"\t"$2-1000"\t"$3+1000"\t"$3-$2}}' |awk '($4>500) {print}'>test1/DEL100_test_output/insertion.query.bed
cat hybird_temp/drep/hybird/MUMandCo/522_1/TD33_98.94_0.0_bin.20.orig.fa-CD100_92.19_2.004_bin.11.strict.fa_output/query.insertion.bed | cut -f1-3 \
|awk '{if($2<1000) {print $1"\t"$2"\t"$3+1000"\t"$3-$2} else {print $1"\t"$2-1000"\t"$3+1000"\t"$3-$2}}' | cut -f1-3 >igv_test/query_insertion.bed
conda activate assembly
make_IGV_snapshots.py  -ht 1000 -o igv_test/MUMandCo_query_insertion_igv_snap -bin /software_users/chenliang/IGV-snapshot-automator/bin/IGV_2.3.81/igv.jar \
-g result/metawrap_bins/hybird/CD100_92.19_2.004_bin.11.strict.fa -r igv_test/query_insertion.bed igv_test/nano.TD33vsCD100.dehost.aln.sorted.bam

minimap2 -ax splice  -t 48 ${ref} temp/minimap2/nano.CD100.dehost.fq.gz > igv_test/nano.CD100vsCD100.dehost.aln.sam
samtools sort igv_test/nano.CD100vsCD100.dehost.aln.sam -o igv_test/nano.CD100vsCD100.dehost.aln.sorted.bam -@ 40
samtools index -@ 40 igv_test/nano.CD100vsCD100.dehost.aln.sorted.bam igv_test/nano.CD100vsCD100.dehost.aln.sorted.bam.bai

conda activate assembly
make_IGV_snapshots.py  -ht 1000 -o igv_test/MUMandCo_query_insertion_igv_snap_2 -bin /software_users/chenliang/IGV-snapshot-automator/bin/IGV_2.3.81/igv.jar \
-g result/metawrap_bins/hybird/CD100_92.19_2.004_bin.11.strict.fa -r igv_test/query_insertion.bed igv_test/nano.CD100vsCD100.dehost.aln.sorted.bam \
igv_test/nano.TD33vsCD100.dehost.aln.sorted.bam


## 3. summary of predicted SVs

# statistic of numbers of predicted SVs
for i in $(ls hybird_temp/drep/hybird/MUMandCo/*/*_output/*.SVs_all.tsv);do
  ref=`echo $i |cut -d'/' -f6 |sed 's/_output//' |cut -d'-' -f1`
  query=`echo $i |cut -d'/' -f6 |sed 's/_output//' |cut -d'-' -f2`

  num1=`sed '1d' ${i} |grep "deletion"  | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l `

  num2=`sed '1d' ${i} |grep "insertion" | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' | cut -f2,5,6,7,8 |awk '{print $1"\t"$4"\t"$5"\t"$2"\t"$3}' |wc -l `

  num3=`sed '1d' ${i} |grep "duplication" | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l `

  num4=`sed '1d' ${i} |grep transloc | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l`
  
  num5=`sed '1d' ${i} |grep "inversion" | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l`

  total=`expr $num1 + $num2 + $num3 + $num4 + $num5`

  echo -e "$ref\t$query\t$num1\t$num2\t$num3\t$num4\t$num5\t$total" >>hybird_temp/drep/hybird/hybird_MUMandCo_SVs_num_all.summary
done
sed -i '1 i ref\tquery\tdeletion_num\tinsertion_num\tduplication_num\ttransloc_num\tinversion_num\ttotal' hybird_temp/drep/hybird/hybird_MUMandCo_SVs_num_all.summary

# statistic of length of all predicted SVs
cat hybird_temp/drep/hybird/MUMandCo/*/*_output/*.SVs_all.tsv | grep "deletion" |awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 >>hybird_temp/drep/hybird/hybird_MUMandCo_SVs_all.summary
cat hybird_temp/drep/hybird/MUMandCo/*/*_output/*.SVs_all.tsv | grep "insertion" | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' | cut -f2,5,6,7,8 |awk '{print $1"\t"$4"\t"$5"\t"$2"\t"$3}' >>hybird_temp/drep/hybird/hybird_MUMandCo_SVs_all.summary
cat hybird_temp/drep/hybird/MUMandCo/*/*_output/*.SVs_all.tsv | grep "inversion" | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 >>hybird_temp/drep/hybird/hybird_MUMandCo_SVs_all.summary
sed -i '1 i contigs\tstart\tend\tlength\ttype' hybird_temp/drep/hybird/hybird_MUMandCo_SVs_all.summary


# statistic of length of >500bp predicted SVs
mkdir -p result/SV_res
s=(query.insertion ref.deletion ref.inversion)
for i in ${s[@]};do
  for j in `ls hybird_temp/drep/hybird/MUMandCo/*/*/${i}.bed`; do
    ref=`echo $j |cut -d'/' -f6 |cut -d'-' -f1 |cut -d'_' -f1`
    query=`echo $j |cut -d'/' -f6 |cut -d'-' -f2 |cut -d'_' -f1`
    awk -v ref="${ref}" -v query="${query}" '{print ref"\t"query"\t"$3-$2+1}' ${j} >>result/SV_res/${i}.length
  done
  sed -i '1 i ref\tquery\tevent_legth' result/SV_res/${i}.length
done
s=(query.insertion ref.deletion ref.inversion)
for i in ${s[@]};do
  for j in `ls hybird_temp/drep/hybird/MUMandCo/*/*/${i}.bed`; do
    ref=`echo $j |cut -d'/' -f6 |cut -d'-' -f1 |cut -d'_' -f1`
    query=`echo $j |cut -d'/' -f6 |cut -d'-' -f2 |cut -d'_' -f1`
    event_num=`wc -l $j |cut -d' ' -f1`
    echo $ref'\t'$query'\t'$event_num >>result/SV_res/${i}.num
  done
  sed -i '1 i ref\tquery\tevent_num' result/SV_res/${i}.num
done
# statistic of numbers of >500bp predicted SVs
rm hybird_temp/drep/hybird/hybird_MUMandCo.summary
for i in $(ls hybird_temp/drep/hybird/MUMandCo/*/*_output/*.SVs_all.tsv);do
  ref=`echo $i |cut -d'/' -f6 |sed 's/_output//' |cut -d'-' -f1`
  query=`echo $i |cut -d'/' -f6 |sed 's/_output//' |cut -d'-' -f2`

  num1=`sed '1d' ${i} |grep "deletion" | awk '($5>500) {print}' | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l `

  num2=`sed '1d' ${i} |grep "insertion" |awk '($5>500) {print}' | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' | cut -f2,5,6,7,8 |awk '{print $1"\t"$4"\t"$5"\t"$2"\t"$3}' |wc -l `

  num3=`sed '1d' ${i} |grep "duplication" | awk '($5>500) {print}' | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l `

  num4=`sed '1d' ${i} |grep transloc | awk '($5>500) {print}' | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l`
  
  num5=`sed '1d' ${i} |grep "inversion" | awk '($5>500) {print}' | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l`

  total=`expr $num1 + $num2 + $num3 + $num4 + $num5`

  echo -e "$ref\t$query\t$num1\t$num2\t$num3\t$num4\t$num5\t$total" >>hybird_temp/drep/hybird/hybird_MUMandCo.summary
done
sed -i '1 i ref\tquery\tdeletion_num\tinsertion_num\tduplication_num\ttransloc_num\tinversion_num\ttotal' hybird_temp/drep/hybird/hybird_MUMandCo.summary


## 4. sequences extraction, genes prediction and annotation of predicted SVs

OUTPUT_DIR="hybird_temp/drep/hybird/MUMandCo"
## extract SV sequences
for i in `sed '1d' hybird_temp/drep/hybird/data_tables/Cdb.csv |cut -f2 -d',' |sort |uniq -c | sort -nr | awk '($1>10) {print}' |awk '{print $2}'`;do
  for j in `sed '1d' hybird_temp/drep/hybird/data_tables/Cdb.csv | awk -v var="${i}" 'BEGIN {FS=","} ($2==var) {print $1}'`;do
    cp result/metawrap_bins/hybird/$j $OUTPUT_DIR/${i}
  done
  ref=`sed '1d' hybird_temp/drep/hybird/data_tables/Widb.csv |awk -v var="${i}" 'BEGIN {FS=","} ($8==var) {print $1}'`
  rm $OUTPUT_DIR/${i}/${ref}
  for m in $(ls $OUTPUT_DIR/${i}/*fa);do 
    m=`echo ${m} |cut -d'/' -f 6`
    cd $OUTPUT_DIR/${i}/

    grep "deletion" ${ref}-${m}_output/${ref}-${m}.SVs_all.tsv | awk '($5>500) {print}' | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6  >${ref}-${m}_output/ref.deletion.bed
    cat ${ref}-${m}_output/ref.deletion.bed | while read line; do
      seqname=`echo ${line} |cut -d' ' -f1`
      start=`echo ${line} |cut -d' ' -f2`
      end=`echo ${line} |cut -d' ' -f3`
      len=`expr ${end} - ${start}`
      awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' /data/nano_meta_ref/result/metawrap_bins/hybird/${ref} \
      | awk -v seqname="${seqname}" -v start="${start}" -v len="${len}" '{if ($1==">"seqname) {print $1"_"start"_"len"\n"substr($2,start,len)}}' \
      >> ${ref}-${m}_output/ref.deletion.fa
#      awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' /data/nano_meta_ref/result/metawrap_bins/hybird/${ref} \
#      | awk -v seqname="${seqname}" -v start="${start}" -v len="${len}" '{if ($1==">"seqname) {print $1"_"start"_"len"\t""1""\t"len}}' | sed 's/>//' \
#      >> ${ref}-${m}_output/ref.deletion.bed
    done

    grep "insertion" ${ref}-${m}_output/${ref}-${m}.SVs_all.tsv |awk '($5>500) {print}' | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' | cut -f2,5,6,7,8 |awk '{print $1"\t"$4"\t"$5"\t"$2"\t"$3}' >${ref}-${m}_output/query.insertion.bed
    cat ${ref}-${m}_output/query.insertion.bed | while read line; do
      seqname=`echo ${line} |cut -d' ' -f1`
      start=`echo ${line} |cut -d' ' -f2`
      end=`echo ${line} |cut -d' ' -f3`
      len=`expr ${end} - ${start}`
      awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' /data/nano_meta_ref/$OUTPUT_DIR/${i}/${m} \
      | awk -v seqname="${seqname}" -v start="${start}" -v len="${len}" '{if ($1==">"seqname) {print $1"_"start"_"len"\n"substr($2,start,len)}}' \
      >> ${ref}-${m}_output/query.insertion.fa
#      awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' /data/nano_meta_ref/$OUTPUT_DIR/${i}/${m}  \
#      | awk -v seqname="${seqname}" -v start="${start}" -v len="${len}" '{if ($1==">"seqname) {print $1"_"start"_"len"\t""1""\t"len}}' | sed 's/>//' \
#      >> ${ref}-${m}_output/query.insertion.bed
    done

    grep "inversion" ${ref}-${m}_output/${ref}-${m}.SVs_all.tsv | awk '($5>500) {print}' | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 >${ref}-${m}_output/ref.inversion.bed
    cat ${ref}-${m}_output/ref.inversion.bed | while read line; do
      seqname=`echo ${line} |cut -d' ' -f1`
      start=`echo ${line} |cut -d' ' -f2`
      end=`echo ${line} |cut -d' ' -f3`
      len=`expr ${end} - ${start}`
      awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' /data/nano_meta_ref/result/metawrap_bins/hybird/${ref} \
      | awk -v seqname="${seqname}" -v start="${start}" -v len="${len}" '{if ($1==">"seqname) {print $1"_"start"_"len"\n"substr($2,start,len)}}' \
      >> ${ref}-${m}_output/ref.inversion.fa
#      awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' /data/nano_meta_ref/result/metawrap_bins/hybird/${ref} \
#      | awk -v seqname="${seqname}" -v start="${start}" -v len="${len}" '{if ($1==">"seqname) {print $1"_"start"_"len"\t""1""\t"len}}' | sed 's/>//' \
#      >> ${ref}-${m}_output/ref.inversion.bed
    done

    cd /data/nano_meta_ref/
  done
done

#gene prediction using prokka
source activate metawrap1.2
for i in $(find hybird_temp/drep/hybird/MUMandCo/*/*/ -name "*.*.fa");do
  type=`echo ${i} | cut -d'/' -f7| sed 's/.fa//' `
  dir=`echo ${i} |sed "s/${type}.fa//"`
  awk 'BEGIN{i=0 ; FS="," ; OFS=","}{ if(/>/){gsub($1,">"++i,$1);print $0}else{print $0}}' ${i} >${dir}/revised.${type}.fa.tem
  prokka ${dir}/revised.${type}.fa.tem --outdir ${dir}/${type} \
  --prefix ${type} --metagenome --kingdom Archaea,Bacteria,Mitochondria,Viruses \
  --force --cpus 40
done
#rename 's/.fa/.fa.tem/' hybird_temp/drep/hybird/MUMandCo/*/*/revised*

#gene annotation using prokka
source activate metagenome
s=(query.insertion ref.deletion ref.inversion)
OUTPUT_DIR="hybird_temp/drep/hybird/MUMandCo"
for j in ${s[@]};do
  mkdir ${OUTPUT_DIR}/../${j}_eggnog
  db="/software/database/"
  time emapper.py -m diamond --no_annot --usemem --no_file_comments --data_dir ${db}/eggnog \
    --cpu 40 -i ${OUTPUT_DIR}/../MUMandCo.morethan10.${j}.faa -o ${OUTPUT_DIR}/../${j}_eggnog/protein --override

  time emapper.py --annotate_hits_table ${OUTPUT_DIR}/../${j}_eggnog/protein.emapper.seed_orthologs --no_file_comments \
          -o ${OUTPUT_DIR}/../${j}_eggnog/output --cpu 40 --data_dir ${db}/eggnog --override

  sed '1 i Name\teggNOG\tEvalue\tScore\tGeneName\tGO\tKO\tBiGG\tTax\tOG\tBestOG\tCOG\tAnnotation' \
    ${OUTPUT_DIR}/../${j}_eggnog/output.emapper.annotations > ${OUTPUT_DIR}/../${j}_eggnog/output

  cut -f 1,7 ${OUTPUT_DIR}/../${j}_eggnog/output |cut -f 1 -d ',' > ${OUTPUT_DIR}/../${j}_eggnog/2ko.list
  awk 'NR==FNR{a[$1]=$2;next}{print $0"\t"a[$4]}' ${OUTPUT_DIR}/../${j}_eggnog/2ko.list ${OUTPUT_DIR}/../MUMandCo.morethan10.${j}.sample2proteinseq.txt \
   |grep -v -P '^\t' |grep -v -P '\t$' >${OUTPUT_DIR}/../${j}_eggnog/2ko.filter.list

done


## 5. Temporal dynamics of structural variations in genome of F.saccharivorans, A. hadrus and A. rectalis within each individual of our time-series cohort

# statistic of MAGs present in >10 samples 
for i in `sed '1d' hybird_temp/drep/hybird/data_tables/Cdb.csv |cut -f2 -d',' |sort |uniq -c | sort -nr | awk '($1>10) {print}' |awk '{print $2}'`;do 
  awk -v var="${i}" 'BEGIN {FS=","} ($2==var) {print}' hybird_temp/drep/hybird/data_tables/Cdb.csv >>hybird_temp/drep/hybird/hybird_more_than_10_samples.Cdb.txt
done

# the reference use the MAGs detected from the sample collected from the day before
# run_TD_10sample_MUMandCo.sh
source activate mummer4
mkdir hybird_temp/drep/hybird/TD_10sample_MUMandCo_new
OUTPUT_DIR="hybird_temp/drep/hybird/TD_10sample_MUMandCo_new"
cat /data/nano_meta_ref/hybird_temp/drep/hybird/ALL_hybird_more_than_10_samplesCdb.SV.txt | awk '/522_1/||/471_2/||/461_1/ {print}' |while read line; do
  arr=($line)
  num=`expr ${#arr[@]} - 2`
  for i in $(seq 1 $num);do 
    echo $i
    tax=`sed '1d' hybird_temp/drep/hybird/data_tables/Widb.csv |awk -v var="${arr[0]}" 'BEGIN {FS=","} ($8==var) {print $1}'`
    ref=`cat /data/nano_meta_ref/hybird_temp/drep/hybird/hybird_more_than_10_samples.Cdb.txt | grep ${arr[$i]}_ |grep ${arr[0]} |cut -d',' -f1`
    query=`cat /data/nano_meta_ref/hybird_temp/drep/hybird/hybird_more_than_10_samples.Cdb.txt | grep ${arr[$i+1]}_ |grep ${arr[0]} |cut -d',' -f1`

    gsize=`seqkit stat -T result/metawrap_bins/hybird/${ref} | sed '1d' |cut -f 5`

    cd $OUTPUT_DIR
    bash /software_users/chenliang/MUMandCo/mumandco_v2.4.2.sh -r /data/nano_meta_ref/result/metawrap_bins/hybird/${ref} \
      -q /data/nano_meta_ref/result/metawrap_bins/hybird/${query} -g ${gsize} -o ${arr[0]}_${arr[$i]}_${arr[$i+1]}

    num1=`sed '1d' ${arr[0]}_${arr[$i]}_${arr[$i+1]}_output/*.SVs_all.tsv |grep "deletion" | awk '($5>500) {print}' | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
      awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l `

    num2=`sed '1d' ${arr[0]}_${arr[$i]}_${arr[$i+1]}_output/*.SVs_all.tsv |grep "insertion" |awk '($5>500) {print}' | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
      awk '{split($2,a,"_");if($8<a[4]-10){print}}' | cut -f2,5,6,7,8 |awk '{print $1"\t"$4"\t"$5"\t"$2"\t"$3}' |wc -l `

    num3=`sed '1d' ${arr[0]}_${arr[$i]}_${arr[$i+1]}_output/*.SVs_all.tsv |grep "duplication" | awk '($5>500) {print}' | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
      awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l `

    num4=`sed '1d' ${arr[0]}_${arr[$i]}_${arr[$i+1]}_output/*.SVs_all.tsv |grep transloc | awk '($5>500) {print}' | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
      awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l`
  
    num5=`sed '1d' ${arr[0]}_${arr[$i]}_${arr[$i+1]}_output/*.SVs_all.tsv |grep "inversion" | awk '($5>500) {print}' | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
      awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l`

    total=`expr $num1 + $num2 + $num3 + $num4 + $num5`

    echo -e "${tax}\t${arr[$i+1]}\t$num1\t$num2\t$num3\t$num4\t$num5\t$total" >>TD_10sample_MUMandCo_var.summary

    cd /data/nano_meta_ref/
  done
done
sed -i '1 i tax\tsample\tdeletion_num\tinsertion_num\tduplication_num\ttransloc_num\tinversion_num\ttotal' hybird_temp/drep/hybird/TD_10sample_MUMandCo_new/TD_10sample_MUMandCo_var.summary



