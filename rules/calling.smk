#!/usr/bin/env python
if "restrict-regions" in config["processing"]:
    rule compose_regions:
        input:
            config["processing"]["restrict-regions"]
        output:
            "called/regions.bed"
        shell:
            "cp {input} {output}"

rule register_gatk3:
    input:
        jar = config["ref"]["gatk3_jar"]
    output:
        touch(".gatk-registered")
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk3-register {input.jar}"

rule call_known_variants:
    input:
        gatk_registered = ".gatk-registered",
        bams = get_all_sample_bams,
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"],
        regions="called/regions.bed",
        known_sites = config["known_sites"]
    output:
        vcf = "genotyped/known/all.vcf.gz"
    params:
        extra = get_call_known_variants_params
    log:
        "logs/gatk/haplotypecaller/known/all.known_sites.log"
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/known_sites.py"


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
        vcf="genotyped/ploidy/all.vcf.gz"
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "logs/gatk/ploidy/genotypegvcfs.log"
    wrapper:
        "0.27.1/bio/gatk/genotypegvcfs"

rule make_diploid_variants:
    input:
        vcf="genotyped/{method}/all.vcf.gz"
    output:
        vcf="diploid/{method}/all.vcf"
    run:
        make_diploid_vcf(input.vcf, output.vcf)

rule filter_unique_variants:
    input:
        known = "diploid/known/all.vcf",
        ploidy = "diploid/ploidy/all.vcf"
    output:
        vcf = "filter_unique/ploidy.vcf"
    run:
        known_variant_dict = {}
        with open(input.known) as f:
            for line in f:
                if line[0] != '#':
                    line_list = line.split()
                    (chrom, pos, x, ref, alt) = line_list[:5]
                    known_variant_dict[(chrom, pos, ref, alt)] = 1
        with open(input.ploidy) as f, open(output.vcf, 'w') as out:
            for line in f:
                if line[0] != '#':
                    line_list = line.split()
                    (chrom, pos, x, ref, alt) = line_list[:5]
                    if not known_variant_dict.get((chrom, pos, ref, alt)):
                        out.write(line)
                else:
                    out.write(line)

rule merge_variants:
    input:
        ref=get_fai(), # fai is needed to calculate aggregation over contigs below
        vcfs= ["filter_unique/ploidy.vcf", "diploid/known/all.vcf"]
    output:
        vcf="nemo_genotypes/all.vcf.gz"
    log:
        "logs/picard/merge-genotyped.log"
    wrapper:
        "0.40.2/bio/picard/mergevcfs"
