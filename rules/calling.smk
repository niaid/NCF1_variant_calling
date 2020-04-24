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

rule call_putative_variants:
    input:
        bam = "merge_recal/{sample}.bam",
        ref = config["ref"]["genome"],
        vcf = config["params"]["putative"]["vcf_for_header"]
    output:
        "putative_variants/{sample}.vcf"
    params:
        chrom = config["params"]["putative"]["chrom"],
        start = config["params"]["putative"]["start"],
        end = config["params"]["putative"]["end"],
        min_alt = config["params"]["putative"]["min_alt"]
    conda:
        "../envs/putative.yaml"
    shell:
        "python scripts/putative_variants.py -b {input.bam} -o {output} -v {input.vcf} "
        "-R {input.ref} -c {params.chrom} -s {params.start} -e {params.end} -m {params.min_alt}"

rule merge_putative_variants:
    input:
        expand("putative_variants/{sample}.vcf", sample=samples.index)
    output:
        "combined_putative/all.vcf"
    run:
        def get_vcf_header(vcf):
            vcf_header = ''
            with open(vcf) as f:
                line = f.readline()
                while not line.startswith('#CHROM'):
                    vcf_header += line
                    line = f.readline()
            vcf_header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
            return vcf_header

        vcf_header = get_vcf_header(input[0])
        variant_dict = {}
        for vcf in input:
            with open(vcf) as f:
                head = f.readline()
                while not head.startswith('#CHROM') and head != '':
                    head = f.readline()
                line = f.readline()
                while line != '':
                    pos = int(line.split()[1])
                    variant_dict[pos] = line
                    line = f.readline()
        with open(output[0], 'w') as out:
            out.write(vcf_header)
            for pos in sorted(variant_dict.keys()):
                line = variant_dict[pos]
                out.write(line)

rule zip_putative:
    input:
        "combined_putative/all.vcf"
    output:
        vcf = "known_sites/putative/all.vcf.gz",
        tbi = "known_sites/putative/all.vcf.gz.tbi"
    conda:
        "../envs/tabix.yaml"
    shell:
        "bgzip -c {input} > {output.vcf};tabix -p vcf {output.vcf}"

rule copy_known_sites:
    input:
        vcf = config["params"]["known_sites"]["vcf"],
        tbi = config["params"]["known_sites"]["vcf"] + ".tbi"
    output:
        vcf = "known_sites/known/all.vcf.gz",
        tbi = "known_sites/known/all.vcf.gz.tbi"
    shell:
        "cp {input.vcf} {output.vcf};cp {input.tbi} {output.tbi}"


rule call_known_variants:
    input:
        gatk_registered = ".gatk-registered",
        bams = get_all_sample_bams,
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"],
        regions="called/regions.bed",
        known_sites = "known_sites/{known}/all.vcf.gz"
    output:
        vcf = "genotyped/{known}/all.vcf.gz"
    params:
        extra = get_call_known_variants_params
    log:
        "logs/gatk/haplotypecaller/{known}/all.known_sites.log"
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
        unpack(get_all_variant_vcfs)
    output:
        expand("filter_unique/{method}/all.vcf", method = calling_methods)
    run:
        all_variant_dict = {}
        for vcf in input.vcfs:
            with open(vcf) as f:
                for line in f:
                    if line[0] != '#':
                        line_list = line.split()
                        (chrom, pos, x, ref, alt) = line_list[:5]
                        all_variant_dict[(chrom, pos, ref, alt)] = 1
        input_dict = get_all_variant_vcfs()
        for method in input_dict.keys():
            in_vcf = input_dict[method]
            out_vcf = "filter_unique/" + method + "/all.vcf"
            with open(in_vcf) as f, open(out_vcf, 'w') as out:
                for line in f:
                    if line[0] == '#':
                        out.write(line)
                    else:
                        line_list = line.split()
                        (chrom, pos, x, ref, alt) = line_list[:5]
                        if all_variant_dict.get((chrom, pos, ref, alt)):
                            out.write(line)
                            del all_variant_dict[(chrom, pos, ref, alt)]

rule merge_variants:
    input:
        ref=get_fai(), # fai is needed to calculate aggregation over contigs below
        vcfs= expand("filter_unique/{method}/all.vcf", method = calling_methods)
    output:
        vcf="nemo_genotypes/all.vcf.gz"
    log:
        "logs/picard/merge-genotyped.log"
    wrapper:
        "0.40.2/bio/picard/mergevcfs"
