#Create tmp directory to hold trimmed files? Can be cleared after pipeline runs to make room
configfile: "./config.yaml"
IN_DIR = config["in_dir"]
SAMPLES = config["samples"]

rule all:
    input: 
        expand("invert_results/{sample}/bam/{sample}_fwd.bam", sample=SAMPLES),
        expand("invert_results/{sample}/bam/{sample}_rev.bam", sample=SAMPLES),
        expand("invert_results/{sample}/cmrna/{sample}_cmratio.txt", sample=SAMPLES)

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
        2
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
        Mapping trimmed reads to host genome
        """
    input:
        r1 = "tmp/{sample}_R1_trimmed.fastq.gz",
        r2 = "tmp/{sample}_R2_trimmed.fastq.gz"
    params:
        annotation = config["annotation_file"], #Combined host-virus annotation file
        star_index = config["genome_index"] #Directory host-virus index file generated by STAR
    output:
        "star_mapping/{sample}Aligned.sortedByCoord.out.bam"
    shell:
        """
        mkdir -p star_mapping/
        STAR \
            --runThreadN 16 \
            --sjdbGTFfile {params.annotation} \
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
            --outFileNamePrefix star_mapping/{wildcards.sample} \
            --outSAMtype BAM SortedByCoordinate \
            --runMode alignReads \
            --genomeDir {params.star_index} \
            --readFilesCommand gunzip -c \
            --readFilesIn {input.r1} {input.r2}
        """

rule cufflinks:
    input:
        bam="star_mapping/{sample}Aligned.sortedByCoord.out.bam"
    output:
        gtf="cufflinks/{sample}/transcripts.gtf",
        iso="cufflinks/{sample}/isoforms.fpkm_tracking",
        genes="cufflinks/{sample}/genes.fpkm_tracking"
    params:
        gtf=config["annotation_file"],
        fasta=config["ref_genome"]
    threads: 4
    shell:
        """
        mkdir -p cufflinks/{wildcards.sample}
        cufflinks \
            {input.bam} \
            --num-threads {threads} \
            -g {params.gtf} \
            -o cufflinks/{wildcards.sample} \
            --library-type fr-firststrand
        """

rule strand_separation:
    input:
        bam="star_mapping/{sample}Aligned.sortedByCoord.out.bam"
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

        scripts/cmRNA_count.sh {input.bam} {wildcards.sample} {output.ratios}
        """