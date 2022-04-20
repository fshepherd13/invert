rule cufflinks:
    input:
        bam=expand('../results/star/{sample}/{sample}_sorted_Q20.bam', sample=SAMPLES)
    output:
        iso="../results/cufflinks/isoforms.fpkm_tracking",
        genes="../results/cufflinks/genes.fpkm_tracking"
    log:
        "logs/cufflinks/cufflinks.log"
    conda:
        "../envs/cufflinks.yaml"
    params:
        gtf=config["annotations"]["combined"],
        sample_list=get_cufflinks_labels(SAMPLES)
    threads:
        8
    shell:
        """
        mkdir -p ../results/cufflinks/
        cuffdiff -p {threads} \
            -T \
            -library-type fr-firststrand \
            -o ../results/cufflinks/ \
            -L {params.sample_list} \
            {params.gtf} \
            {input.bam} &> {log}
        """