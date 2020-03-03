#!/usr/bin/env python
if "restrict-regions" in config["processing"]:
    rule compose_regions:
        input:
            config["processing"]["restrict-regions"]
        output:
            "called/regions.bed"
        shell:
            "cp {input} {output}"

rule call_known_variants:
    input:
        bams = get_all_sample_bams,
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"],
        regions="called/regions.bed",
        known_sites = config["known_sites"]
    output:
        vcf = "genotyped/known/all.vcf.gz"
    params:
        extra = get_call_known_variants_params,
        java_opts = ""
    conda:
        "../envs/gatk.yaml"
    run:
        known = "--dbsnp " + input.known
        bams = list(map("-I {}".format, input.bams))
        shell(
            "gatk --java-options '{params.java_opts}' HaplotypeCaller {params.extra} "
            "-R {input.ref} {bams} "
            "-O {output.vcf} {known} {log}"
        )
    


rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"],
        regions="called/regions.bed" if config["processing"].get("restrict-regions") else []
    output:
        gvcf=protected("called/ploidy/{sample}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/ploidy/{sample}.log"
    params:
        extra=get_call_ploidy_variants_params
    wrapper:
        "0.27.1/bio/gatk/haplotypecaller"


rule combine_calls:
    input:
        ref=config["ref"]["genome"],
        gvcfs=expand("called/ploidy/{sample}.g.vcf.gz", sample=samples.index)
    output:
        gvcf="called/ploidy/all.g.vcf.gz"
    log:
        "logs/gatk/ploidy/combinegvcfs.log"
    wrapper:
        "0.27.1/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf="called/ploidy/all.g.vcf.gz"
    output:
        vcf=temp("genotyped/ploidy/all.vcf.gz")
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
