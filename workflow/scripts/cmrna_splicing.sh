#!/bin/bash

SAMPLES=$1 #Assign BAM file provided by snakemake rule to variable (path/to/sample_fwd.bam)
FILE=$2 #Assign the sample name provided by snakemake rule (sample)
OUTPUT=$3 #Assign output file path provided by snakemake rule (path/to/sample_splicing_count.txt)


## Calculate the number of reads mapped to SA/SD on the left. 
#Number of reads from  both unspliced and spliced transcript (i.e. position at the 5' side of the splice junction)
T1_M=$(samtools depth -a -d 0 -r EF467824.1:51-51 $SAMPLES)
T1_NS=$(samtools depth -a -d 0 -r EF467817.1:56-56 $SAMPLES)	
#Number of reads from unspliced transcript (i.e. the position at the 3' side of the splice junction, where no spliced transcript is present).
U1_M=$(samtools depth -a -d 0 -r EF467824.1:52-52 $SAMPLES)
U1_NS=$(samtools depth -a -d 0 -r EF467817.1:57-57 $SAMPLES)

## Calculate the number of reads mapped to SA/SD on the right
#Number of reads from  both unspliced and spliced transcript (i.e. position at the 3' side of the splice junction)
T2_M=$(samtools depth -a -d 0 -r EF467824.1:740-740 $SAMPLES)
T2_NS=$(samtools depth -a -d 0 -r EF467817.1:529-529 $SAMPLES)
#Number of reads from unspliced transcript (i.e. position at the 5' side of the splice junction, where no spliced transcript is present).
U2_M=$(samtools depth -a -d 0 -r EF467824.1:739-739 $SAMPLES)
U2_NS=$(samtools depth -a -d 0 -r EF467817.1:528-528 $SAMPLES)

## Print out results to temporary file
echo $T1_M	$U1_M	$T2_M	$U2_M   > $(dirname $OUTPUT)/tmp.txt
echo $T1_NS	$U1_NS	$T2_NS	$U2_NS >> $(dirname $OUTPUT)/tmp.txt

#Awk to put it into usable format
awk 'BEGIN {print "Gene\tDepth_total_left\tDepth_unspliced_left\tDepth_total_right\tDepth_unspliced_right";}
 {print $1,"\t",$3,"\t",$6,"\t",$9,"\t",$12;}' $(dirname $OUTPUT)/tmp.txt > $OUTPUT

#Remove tmp file
rm $(dirname $OUTPUT)/tmp.txt