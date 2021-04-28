#Create tmp directory to hold trimmed files? Can be cleared after pipeline runs to make room
configfile = config.yaml

rule all:
    input: "{sample}_Aligned.sortedByCoord.out.bam"

rule trimmomatic_pe:
    message:
        """
        Pre-processing raw reads with trimmomatic. Trimming low quality reads and adapter sequences. Running QC on trimmed reads.
        """"
    input:
        r1 = expand("{in_dir}/{{sample}}_R1_001.fastq.gz", in_dir=config["in_dir"]),
        r2 = expand("{in_dir}/{{sample}}_R2_001.fastq.gz", in_dir=config["in_dir"])
    params:
        extra = config["paramters"]["illuminaclip"],
        trimmer = ["SLIDINGWINDOW:4:20"]
    output:
        r1 = "tmp/{sample}_R1_trimmed.fastq.gz",
        r2 = "tmp/{sample}_R2_trimmed.fastq.gz",
        r1_unpaired = "tmp/{sample}_R1_unpaired_trimmed.fastq.gz",
        r2_unpaired = "tmp/{sample}_R2_unpaired_trimmed.fastq.gz"
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

