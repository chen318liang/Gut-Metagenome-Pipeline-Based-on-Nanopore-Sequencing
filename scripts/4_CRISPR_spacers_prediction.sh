
#####################
Part IV. Predicting CRISPR spacers from the contigs of all samples
###################
#!/usr/bin/bash
## CRISPRDetect (https://github.com/davidchyou/CRISPRDetect_2.4)
## crisprdetectparser.py (https://github.com/hwalinga/crisprdetect-parser)

path="/data/nano_meta_ref/result/spades/hybird_temp/spades"

perl CRISPRDetect.pl -f $path/contigs.fasta -o CRISPRDetect_result \
-check_direction 0 -array_quality_score_cutoff 3 -T 20

python crisprdetectparser.py --spacers-directory spacer_dir --spacers-extension fna CRISPRDetect_result > metadata.tsv

