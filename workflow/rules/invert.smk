rule strand_separation:
    input:
        bam='../results/star/{sample}/{sample}_sorted_Q20.bam'
    output:
        fwd1_bam = temp("../results/invert/{sample}/bam/{sample}_fwd1.bam"),
        fwd2_bam = temp("../results/invert/{sample}/bam/{sample}_fwd2.bam"),
        fwd_bam = "../results/invert/{sample}/bam/{sample}_fwd.bam",
        rev1_bam = temp("../results/invert/{sample}/bam/{sample}_rev1.bam"),
        rev2_bam = temp("../results/invert/{sample}/bam/{sample}_rev2.bam"),
        rev_bam = "../results/invert/{sample}/bam/{sample}_rev.bam"
    conda:
        "../envs/invert.yaml"
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
        bam="../results/invert/{sample}/bam/{sample}_fwd.bam"
    output:
        ratios="../results/invert/{sample}/cmrna/{sample}_cmratio.txt"
    conda:
        "../envs/invert.yaml"
    shell:
        """
        mkdir -p $(dirname {output.ratios}/)

        scripts/cmrna_count.sh {input.bam} {wildcards.sample} {output.ratios}
        """

rule splice_ratios:
    input:
        bam="../results/invert/{sample}/bam/{sample}_fwd.bam"
    output:
        splice_counts="../results/invert/{sample}/splice_counts/{sample}_splicing_count.txt"
    conda:
        "../envs/invert.yaml"
    shell:
        """
        mkdir -p $(dirname {output.splice_counts}/)
        scripts/cmrna_splicing.sh {input.bam} {wildcards.sample} {output.splice_counts}
        """

rule kinetics_calculations:
    input:
        ratios="../results/invert/{sample}/cmrna/{sample}_cmratio.txt",
        splice_counts="../results/invert/{sample}/splice_counts/{sample}_splicing_count.txt",
        expression_levels="../results/cufflinks/genes.fpkm_tracking"
    output:
        "../results/invert/{sample}/{sample}_final_results.csv"
    shell:
        """
        python scripts/vcmRNA_calculation.py {input.ratios} {input.splice_counts} {input.expression_levels} {wildcards.sample} {output}
        """
