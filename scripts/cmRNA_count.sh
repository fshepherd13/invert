#!/bin/bash

SAMPLES=$1 #Assign BAM file provided by snakemake rule to variable (path/to/sample_fwd.bam)
FILE=$2 #Assign the sample name provided by snakemake rule (sample)
OUTPUT=$3 #Assign output file path provided by snakemake rule (path/to/sample_cmratio.txt)

samtools view $SAMPLES Cal_PB2:2325-2341  -o $(dirname $OUTPUT)/PB2.bam
samtools view $SAMPLES Cal_PB1:2326-2341  -o $(dirname $OUTPUT)/PB1.bam
samtools view $SAMPLES Cal_PA:2217-2233  -o $(dirname $OUTPUT)/PA.bam
samtools view $SAMPLES Cal_HA:1761-1777  -o $(dirname $OUTPUT)/HA.bam
samtools view $SAMPLES Cal_NP:1549-1565 -o $(dirname $OUTPUT)/NP.bam
samtools view $SAMPLES Cal_NA:1443-1458 -o $(dirname $OUTPUT)/NA.bam
samtools view $SAMPLES Cal_M:1011-1027 -o $(dirname $OUTPUT)/M.bam
samtools view $SAMPLES Cal_NS:874-890  -o $(dirname $OUTPUT)/NS.bam
samtools merge -f $(dirname $OUTPUT)/${FILE}_fwd3prime.bam $(dirname $OUTPUT)/{PB2,PB1,PA,HA,NP,NA,M,NS}.bam
samtools index $(dirname $OUTPUT)/${FILE}_fwd3prime.bam

cPB2=$(samtools view $(dirname $OUTPUT)/PB2.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAACG"|wc -l) 
cPB1=$(samtools view $(dirname $OUTPUT)/PB1.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAATG"|wc -l)
cPA=$(samtools view $(dirname $OUTPUT)/PA.bam| awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAGT"|wc -l)
cNP=$(samtools view $(dirname $OUTPUT)/NP.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAATA"|wc -l)
cHA=$(samtools view $(dirname $OUTPUT)/HA.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAACA"|wc -l)
cNA=$(samtools view $(dirname $OUTPUT)/NA.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAACT"|wc -l)
cM=$(samtools view $(dirname $OUTPUT)/M.bam| awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAACT"|wc -l)
cNS=$(samtools view $(dirname $OUTPUT)/NS.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAACA"|wc -l)

mPB2=$(samtools view $(dirname $OUTPUT)/PB2.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
mPB1=$(samtools view $(dirname $OUTPUT)/PB1.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAAA"|wc -l)
mPA=$(samtools view $(dirname $OUTPUT)/PA.bam| awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
mNP=$(samtools view $(dirname $OUTPUT)/NP.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
mHA=$(samtools view $(dirname $OUTPUT)/HA.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
mNA=$(samtools view $(dirname $OUTPUT)/NA.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
mM=$(samtools view $(dirname $OUTPUT)/M.bam| awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
mNS=$(samtools view $(dirname $OUTPUT)/NS.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)

echo -e "gene\tmRNA\tcRNA" > $OUTPUT
echo -e "PB2\t$mPB2\t$cPB2" >> $OUTPUT
echo -e "PB1\t$mPB1\t$cPB1" >> $OUTPUT
echo -e "PA\t$mPA\t$cPA" >> $OUTPUT
echo -e "HA\t$mHA\t$cHA" >> $OUTPUT
echo -e "NP\t$mNP\t$cNP" >> $OUTPUT 
echo -e "NA\t$mNA\t$cNA" >> $OUTPUT
echo -e "M\t$mM\t$cM" >> $OUTPUT
echo -e "NS\t$mNS\t$cNS" >> $OUTPUT

awk 'NR==1{$5="mRNA:total_pos_RNA"}NR>1{if($2+$3==0) $4 = "N/A"; else $4=$2/($2+$3)}{print}' $OUTPUT > testfile.tmp && mv testfile.tmp $OUTPUT 