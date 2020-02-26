#!/usr/bin/env python

rule subset_bam:
    input:
        get_start_bam
    output:
        bam = "subset/{sample}-{unit}.bam"
    params:
        get_rg_subset_param
    wrapper:
        "0.49.0/bio/samtools/view"

rule subset_read_list:
    input:
        bam = "subset/{sample}-{unit}.bam"
    output:
        "subset/{sample}-{unit}.bam_reads.txt"
    conda:
        "../envs/subset.yaml"
    shell:
        "samtools view {input.bam} | cut -f1 > {output}"

rule keep_paired_reads:
    input:
        "subset/{sample}-{unit}.bam_reads.txt"
    output:
        "subset/{sample}-{unit}.bam_reads_paired.txt"
    run:
        read_dict = {}
        with open(input[0]) as f:
            for line in f:
                read = line.rstrip()
                if not read_dict.get(read):
                    read_dict[read] = 0
                read_dict[read] += 1
        with open(output[0], 'w') as out:
            for read in read_dict.keys():
                read_count = read_dict[read]
                if read_count == 2:
                    out.write(read + '\n')

rule filter_sam_reads:
    input:
        bam = "subset/{sample}-{unit}.bam",
        txt = "subset/{sample}-{unit}.bam_reads_paired.txt"
    output:
        "subset_paired/{sample}-{unit}.bam"
    conda:
        "../envs/subset.yaml"
    shell:
        "picard FilterSamReads I={input.bam} O={output} READ_LIST_FILE={input.txt} FILTER=includeReadList"

rule samtofastq:
    input:
        "subset_paired/{sample}-{unit}.bam"
    output:
        fastq1 = "samtofastq/{sample}-{unit}.1.fq",
        fastq2 = "samtofastq/{sample}-{unit}.2.fq"
    log:
        "logs/picard/samtofastq/{sample}-{unit}.log"
    wrapper:
        "0.49.0/bio/picard/samtofastq"
