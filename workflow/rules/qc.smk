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
