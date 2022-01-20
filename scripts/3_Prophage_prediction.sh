#####################
Part III. Prediction of prophage, redundancy removal, gene annotation, and construction of evolutionary trees
###################
# prophage dientification using the machine-learning-based tool ProphageHunter
# Webserver address
https://pro-hunter.genomics.cn/

# Candidate propages clustering by CDhit
CDhit-est -i Candidate_propahge.fasta -c 0.95 -o Candidate_propahge_clusted.fasta

# Viral genomic CDS prediction using multiphate2 (https://github.com/carolzhou/multiPhATE2)
#Set up the configuration file as required and then run cmd as follow
python3 multiPhate.py multiPhate.config


# Create phylogenetic tree using iqtree2 (https://github.com/iqtree/iqtree2)
#Join_MCP_TLS.faa is joined protein sequence with major capsid protein and terminal large subunit 
iqtree2 -s Join_MCP_TLS.faa -m MFP -B 1000 --bnni -T 40

