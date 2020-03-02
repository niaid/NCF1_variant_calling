#!/usr/bin/env python
if "restrict-regions" in config["processing"]:
    rule compose_regions:
        input:
            config["processing"]["restrict-regions"]
        output:
            "called/regions.bed"
        shell:
            "cp {input} {output}"


rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"],
        regions="called/regions.bed" if config["processing"].get("restrict-regions") else []
    output:
        gvcf=protected("called/{method}/{sample}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{method}/{sample}.log"
    params:
        extra=get_call_variants_params
    wrapper:
        "0.27.1/bio/gatk/haplotypecaller"


rule combine_calls:
    input:
        ref=config["ref"]["genome"],
        gvcfs=expand("called/{{method}}/{sample}.g.vcf.gz", sample=samples.index)
    output:
        gvcf="called/{method}/all.g.vcf.gz"
    log:
        "logs/gatk/{method}/combinegvcfs.log"
    wrapper:
        "0.27.1/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf="called/{method}/all.g.vcf.gz"
    output:
        vcf=temp("genotyped/{method}/all.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "logs/gatk/{method}/genotypegvcfs.log"
    wrapper:
        "0.27.1/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        ref=get_fai(), # fai is needed to calculate aggregation over contigs below
        vcfs=lambda w: expand("genotyped/{method}/all.vcf.gz", method = calling_methods),
    output:
        vcf="nemo_genotypes/all.vcf.gz"
    log:
        "logs/picard/merge-genotyped.log"
    wrapper:
        "0.40.2/bio/picard/mergevcfs"
