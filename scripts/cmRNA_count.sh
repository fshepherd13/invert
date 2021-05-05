#!/bin/bash

SAMPLES=$1 #Assign BAM file provided by snakemake rule to variable
FILE=$2 #Assign the sample name provided by snakemake rule
OUTPUT=$3 #Assign output file path provided by snakemake rule

	#mkdir ${FILE}_cmRNA
     #   cd ./${FILE}_cmRNA
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

	cPB2=$(samtools view PB2.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAACG"|wc -l) 
	cPB1=$(samtools view PB1.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAATG"|wc -l)
	cPA=$(samtools view PA.bam| awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAGT"|wc -l)
	cNP=$(samtools view NP.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAATA"|wc -l)
	cHA=$(samtools view HA.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAACA"|wc -l)
	cNA=$(samtools view NA.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAACT"|wc -l)
	cM=$(samtools view M.bam| awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAACT"|wc -l)
	cNS=$(samtools view NS.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAACA"|wc -l)
	
	mPB2=$(samtools view PB2.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
	mPB1=$(samtools view PB1.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAAA"|wc -l)
	mPA=$(samtools view PA.bam| awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
	mNP=$(samtools view NP.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
	mHA=$(samtools view HA.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
	mNA=$(samtools view NA.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
	mM=$(samtools view M.bam| awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
	mNS=$(samtools view NS.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)

	echo ${FILE} > $OUTPUT
	echo " " "mRNA" "cRNA" >> $OUTPUT
	echo "PB2" $mPB2 $cPB2 >> $OUTPUT
	echo "PB1" $mPB1 $cPB1 >> $OUTPUT
	echo "PA" $mPA $cPA >> $OUTPUT
	echo "HA" $mHA $cHA >> $OUTPUT
	echo "NP" $mNP $cNP >> $OUTPUT 
	echo "NA" $mNA $cNA >> $OUTPUT
	echo "M" $mM $cM >> $OUTPUT
	echo "NS" $mNS $cNS >> $OUTPUT

    #awk '{$4=$2/($2+$3); $5=$3/($2+$3); print}' $OUTPUT