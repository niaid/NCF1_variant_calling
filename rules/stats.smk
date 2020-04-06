rule vcf_to_tsv:
    input:
        "annotated/all.vcf.gz"
    output:
        report("tables/calls.tsv.gz", caption="../report/calls.rst", category="Calls")
    conda:
        "../envs/rbt.yaml"
    shell:
        "bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%ANN\t[ %GT]\t[ %DP]\t[ %AD]\n' {input} | "
        "gzip > {output}"


rule plot_stats:
    input:
        "tables/calls.tsv.gz"
    output:
        depths=report("plots/depths.svg", caption="../report/depths.rst", category="Plots"),
        freqs=report("plots/allele-freqs.svg", caption="../report/freqs.rst", category="Plots")
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot-depths.py"
