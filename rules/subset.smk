rule subset_bam:
    input:
        get_input_bam
    output:
        "subset/{sample}.bam"
    log:
        "logs/subset/{sample}.log"
    params:
        "-hb -L " + config['restrict-regions']
    wrapper:
        "0.45.1/bio/samtools/view"

rule bam_to_fastq:
    input:
        "subset/{sample}.bam"
    output:
        get_bam_to_fastq
    conda:
        "../envs/subset.yaml"
    shell:
        "picard SamToFastq I={input} OUTPUT_PER_RG=true OUTPUT_DIR=bam_to_fastq"
