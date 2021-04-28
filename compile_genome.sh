# Bash Script for Hybrid Human-influenza Cal09 CoV2 Pipeline
# Based off of script at https://github.com/heznanda/scrnaseq-hybrid-cov2/blob/master/hybrid_pipeline_bash_script.sh
# Modified by Frances Shepherd, 26-April-2021

# 1. Retrieve Cal09 specific genome files and create annotation file (done by hand). 
# Cal_HA.fa
# Cal_M.fa
# Cal_NA.fa
# Cal_NP.fa
# Cal_NS.fa
# Cal_PA.fa
# Cal_PB1.fa
# Cal_PB2.fa
# Annotation file: IAV_Cal09_spliced.gtf

# 2. Download human genome reference files
# fasta nucleodic acid
wget http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz #fasta
wget http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz  # gff
# uncompress files
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh38.103.gtf.gz 

# 3 Concat genome and annotation files
cat *.fa > GRCh38.103_Cal09.fa

cat *.gtf > GRCh38.103_Cal09.gtf 