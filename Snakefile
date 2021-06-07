import pandas as pd

#Create tmp directory to hold trimmed files? Can be cleared after pipeline runs to make room
configfile: "./config.yaml"
IN_DIR = config["in_dir"]
SAMPLES = pd.read_table(config["sample_file"])['Sample']

rule all:
    input: 
        expand("invert_results/{sample}/{sample}_final_results.csv", sample=SAMPLES),
        expand("cufflinks/{sample}/transcripts.gtf", sample=SAMPLES),
        expand("star_mapping/{sample}/Aligned.sortedByCoord.out.bam", sample = SAMPLES),
        expand("invert_results/{sample}/bam/{sample}_fwd.bam", sample=SAMPLES),
        expand("invert_results/{sample}/bam/{sample}_rev.bam", sample=SAMPLES),
        expand("invert_results/{sample}/cmrna/{sample}_cmratio.txt", sample=SAMPLES),
        expand("invert_results/{sample}/splice_counts/{sample}_splicing_count.txt", sample=SAMPLES),
        expand("tmp/{sample}_{direction}_trimmed.fastq.gz", sample=SAMPLES, direction = ["R1", "R2"]),
        expand("qc/{sample}_{direction}_trimmed_fastqc.html", sample=SAMPLES, direction = ["R1", "R2"]),


rule trimmomatic_pe:
    message:
        """
        Pre-processing raw reads with trimmomatic. Trimming low quality reads and adapter sequences. Running QC on trimmed reads.
        """
    input:
        r1 = f"{IN_DIR}/{{sample}}_R1_001.fastq.gz",
        r2 = f"{IN_DIR}/{{sample}}_R2_001.fastq.gz"
    params:
        trimmer = config["parameters"]["trim"],
        extra = ""
    output:
        r1 = "tmp/{sample}_R1_trimmed.fastq.gz",
        r2 = "tmp/{sample}_R2_trimmed.fastq.gz",
        r1_unpaired = "tmp/{sample}_R1_unpaired_trimmed.fastq.gz",
        r2_unpaired = "tmp/{sample}_R2_unpaired_trimmed.fastq.gz"
    threads:
        6
    wrapper:
        "0.74.0/bio/trimmomatic/pe"

rule fastqc:
    input:
        r1="tmp/{sample}_R1_trimmed.fastq.gz",
        r2="tmp/{sample}_R2_trimmed.fastq.gz"
    output:
        r1="qc/{sample}_R1_trimmed_fastqc.html",
        r2="qc/{sample}_R2_trimmed_fastqc.html"
    shell:
        """
        #!/bin/bash
        fastqc {input.r1} {input.r2} -q -o qc/
        """

rule map_reads:
    message:
        """
        Mapping trimmed reads from {wildcards.sample} to host genome
        """
    input:
        fq1 = "tmp/{sample}_R1_trimmed.fastq.gz",
        fq2 = "tmp/{sample}_R2_trimmed.fastq.gz"
    params:
        index= config["genome_index"],
        extra=f"--sjdbGTFfile ref_files/cal09_human/GRCh38.103_Cal09.gtf \
            --sjdbOverhang 149 \
            --outFilterType BySJout \
            --outFilterMultimapNmax 10 \
            --alignSJoverhangMin 5 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
            --outSAMtype BAM SortedByCoordinate \
            --runMode alignReads"
    output:
        "star_mapping/{sample}/Aligned.sortedByCoord.out.bam"
    log:
        "logs/star/{sample}.log"
    threads: 20
    wrapper:
        "0.74.0/bio/star/align"

rule cufflinks:
    input:
        bam="star_mapping/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        gtf="cufflinks/{sample}/transcripts.gtf",
        iso="cufflinks/{sample}/isoforms.fpkm_tracking",
        genes="cufflinks/{sample}/genes.fpkm_tracking"
    params:
        gtf=config["annotations"]["iav_only"]
    threads: 6
    shell:
        """
        mkdir -p cufflinks/{wildcards.sample}
        cufflinks \
            {input.bam} \
            --num-threads {threads} \
            -G {params.gtf} \
            -o cufflinks/{wildcards.sample} \
            --library-type fr-firststrand
        """

rule strand_separation:
    input:
        bam="star_mapping/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        fwd1_bam = "invert_results/{sample}/bam/{sample}_fwd1.bam",
        fwd2_bam = "invert_results/{sample}/bam/{sample}_fwd2.bam",
        fwd_bam = "invert_results/{sample}/bam/{sample}_fwd.bam",
        rev1_bam = "invert_results/{sample}/bam/{sample}_rev1.bam",
        rev2_bam = "invert_results/{sample}/bam/{sample}_rev2.bam",
        rev_bam = "invert_results/{sample}/bam/{sample}_rev.bam"
    shell:
        """
        #Forward strand
        # 1. alignments of the second in pair if they map to the forward strand
        samtools view -b -f 128 -F 16 {input.bam} > {output.fwd1_bam}

        # 2. alignments of the first in pair if they map to the reverse  strand 
        samtools view -b -f 80 {input.bam} > {output.fwd2_bam}
        
        # Combine alignments that originate on the forward strand.
        samtools merge -f {output.fwd_bam} {output.fwd1_bam} {output.fwd2_bam}
        samtools index {output.fwd_bam}

        # Reverse strand
        # 1. alignments of the second in pair if they map to the reverse strand
        samtools view -b -f 144 {input.bam} > {output.rev1_bam}

        # 2. alignments of the first in pair if they map to the forward strand
        samtools view -b -f 64 -F 16 {input.bam} > {output.rev2_bam}

        # Combine alignments that originate on the reverse strand.
        samtools merge -f {output.rev_bam} {output.rev1_bam} {output.rev2_bam}
        samtools index {output.rev_bam}
        """

rule cmRNA_count:
    input:
        bam="invert_results/{sample}/bam/{sample}_fwd.bam"
    output:
        ratios="invert_results/{sample}/cmrna/{sample}_cmratio.txt"
    shell:
        """
        mkdir -p $(dirname {output.ratios}/)

        scripts/cmrna_count.sh {input.bam} {wildcards.sample} {output.ratios}
        """

rule splice_ratios:
    input:
        bam="invert_results/{sample}/bam/{sample}_fwd.bam"
    output:
        splice_counts="invert_results/{sample}/splice_counts/{sample}_splicing_count.txt"
    shell:
        """
        mkdir -p $(dirname {output.splice_counts}/)
        scripts/cmrna_splicing.sh {input.bam} {wildcards.sample} {output.splice_counts}
        """

rule kinetics_calculations:
    input:
        ratios="invert_results/{sample}/cmrna/{sample}_cmratio.txt",
        splice_counts="invert_results/{sample}/splice_counts/{sample}_splicing_count.txt",
        expression_levels="cufflinks/{sample}/genes.fpkm_tracking"
    output:
        "invert_results/{sample}/{sample}_final_results.csv"
    shell:
        """
        python scripts/vcmRNA_calculation.py {input.ratios} {input.splice_counts} {input.expression_levels} {output}
        """
