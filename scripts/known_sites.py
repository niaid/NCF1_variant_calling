import os

from snakemake.shell import shell

known = snakemake.input.get("known", "")
if known:
    known = "--dbsnp " + known
sites_vcf = snakemake.input.get("known_sites", "")
if sites_vcf:
    sites_vcf = "--alleles " + sites_vcf

extra = snakemake.params.get("extra", "")
java_opts = snakemake.params.get("java_opts", "")
bams = snakemake.input.bams
if isinstance(bams, str):
    bams = [bams]
bams = list(map("-I {}".format, bams))

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
shell(
    "gatk3 {java_opts} -T HaplotypeCaller {extra} {sites_vcf} "
    "-R {snakemake.input.ref} {bams} "
    "-o {snakemake.output.vcf} {known} {log}"
)
