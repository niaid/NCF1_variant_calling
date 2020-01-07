rule subset_bam:
    input:
        get_input_bam
    output:
        "subset/{sample}.bam"
    log:
        "logs/subset/{sample}.log"
    