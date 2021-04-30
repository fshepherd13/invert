#Create tmp directory to hold trimmed files? Can be cleared after pipeline runs to make room
configfile: "./config.yaml"
IN_DIR = config["in_dir"]
SAMPLES = config["samples"]

rule all:
    input: 
        expand("assembly/{sample}/{sample}_transcripts.gtf", sample=SAMPLES)

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

rule map_reads:
    message:
        """
        Mapping trimmed reads to host genome
        """
    input:
        r1 = "tmp/{sample}_R1_trimmed.fastq.gz",
        r2 = "tmp/{sample}_R2_trimmed.fastq.gz"
    params:
        annotation = config["annotation_file"]
    output:
        "{sample}_Aligned.sortedByCoord.out.bam"
    shell:
        """
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
            --outFileNamePrefix {wildcards.sample}_ \
            --outSAMtype BAM SortedByCoordinate \
            --runMode alignReads \
            --genomeDir ./index \
            --readFilesIn {input.r1} {input.r2}
        """

rule cufflinks:
    input:
        bam="{sample}_Aligned.sortedByCoord.out.bam"
    output:
        gtf="assembly/{sample}/{sample}_transcripts.gtf",
    params:
        gtf=config["annotation_file"]
    threads: 4
    shell:
        """
        cufflinks \
            {input.bam} \
            --num-threads {threads} \
            -g {params.gtf} \
            -o assembly/{wildcards.sample} \
        """

