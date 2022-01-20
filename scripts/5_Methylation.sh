
#####################
Part V. Determination of DNA methylation using tombo software in genomes with more than 10 samples present
###################
##Specific Alternate Base Detection usting Tombo software (https://github.com/nanoporetech/tombo)
path="/yidong1/BJXWZ-202001013C-01/raw_data/genome/Nanopore"

multi_to_single_fast5 -i $path/fast5_pass/ -s fast5_pass_sep/ --recursive -t 12

# re-squiggle raw reads
tombo resquiggle fast5_pass_sep/ drep_10_genome.fasta  --processes 12

# run modified base detection
tombo detect_modifications alternative_model --fast5-basedirs fast5_pass_sep/ --alternate-bases dam \
--statistics-file-basename all.sample.alt_modified_base_detection \
--per-read-statistics-basename all.pre_reads_statistics_basename --processes 12

# output to genome browser compatible format
tombo text_output browser_files --fast5-basedirs fast5_pass_sep/ --statistics-filename \
all.sample.alt_modified_base_detection.dam.tombo.stats --browser-file-basename \
all.sample.alt_modified_base_detection.dam.sample_alt_model --file-types dampened_fraction coverage

