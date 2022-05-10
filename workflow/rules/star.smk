rule star:
    input:
        fq1 = rules.trim.output.r1_paired,
        fq2 = rules.trim.output.r2_paired
    params:
        index = config["genome_index"],
        annotation = config["annotations"]["combined"],
        outdir = '../results/star/{sample}'
    output:
        bam = '../results/star/{sample}/{sample}_Aligned.sortedByCoord.out.bam'
    log:
        "logs/star/{sample}.log"
    conda:
        "../envs/star.yaml"
    threads: 
        16
    shell:
        """
        STAR --runThreadN {threads} \
         --sjdbGTFfile {params.annotation} \
         --sjdbOverhang 149 \
         --outFilterType BySJout \
         --outFilterMultimapNmax 20 \
         --readFilesCommand zcat \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outFilterScoreMinOverLread 0 \
         --outFilterMatchNminOverLread 0 \
         --outFilterMatchNmin 0 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
         --outSAMtype BAM SortedByCoordinate \
         --runMode alignReads \
         --genomeDir {params.index} \
         --readFilesIn {input.fq1} {input.fq2} \
         --outFileNamePrefix {params.outdir}/{wildcards.sample}_ &> {log}

        rm -rf '../results/star/*/*_STARtmp'
        """

rule star_sort:
    input:
        '../results/star/{sample}/{sample}_Aligned.sortedByCoord.out.bam'
    output:
        bam = '../results/star/{sample}/{sample}_sorted_Q20.bam',
        index = '../results/star/{sample}/{sample}_sorted_Q20.bam.bai'
    params:
        indir = "../results/star/{sample}"
    conda:
        "../envs/invert.yaml"
    shell:
        '''
        samtools flagstat {input} > ../results/star/{wildcards.sample}/{wildcards.sample}_stats.txt
        samtools view -b -h -q 20 -o {output.bam} {input}
        samtools index {params.indir}/{wildcards.sample}_sorted_Q20.bam
        samtools flagstat {params.indir}/{wildcards.sample}_sorted_Q20.bam > {params.indir}/{wildcards.sample}_sorted_Q20_stats.txt
        '''