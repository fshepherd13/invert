rule cufflinks:
    input:
        bam=rules.star.output
    output:
        gtf="../results/cufflinks/{sample}/transcripts.gtf",
        iso="../results/cufflinks/{sample}/isoforms.fpkm_tracking",
        genes="../results/cufflinks/{sample}/genes.fpkm_tracking"
    log:
        "logs/cufflinks/{sample}.log"
    conda:
        "../envs/cufflinks.yaml"
    params:
        gtf=config["annotations"]["iav_only"]
    threads: 24
    shell:
        """
        mkdir -p ../results/cufflinks/{wildcards.sample}
        cufflinks \
            {input.bam} \
            --num-threads {threads} \
            -G {params.gtf} \
            -o ../results/cufflinks/{wildcards.sample} \
            --library-type fr-firststrand
        """