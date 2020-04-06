#!/usr/bin/env python

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


rule update_tsv_header:
    input:
        "tables/calls.tsv.gz"
    output:
        "tables/calls_update.tsv.gz"
    run:
        with gzip.open(input[0]) as f, gzip.open(output[0], 'wb') as out:
            head = f.readline()
            head_list = head.decode().split()[1:]
            new_head1 = []
            new_head2 = []
            for i in head_list:
                fields = i.split(']')[1].split(':')
                if len(fields) == 1:
                    new_head1.append('VARIANT')
                    new_head2.append(fields[0])
                else:
                    new_head1.append(fields[0])
                    new_head2.append(fields[1])
            out.write('\t'.join(new_head1) + '\n')
            out.write('\t'.join(new_head2) + '\n')
            line = f.readline()
            while line.decode() != '':
                out.write(line.decode())
                line = f.readline()


rule plot_stats:
    input:
        "tables/calls_update.tsv.gz"
    output:
        depths=report("plots/depths.svg", caption="../report/depths.rst", category="Plots"),
        freqs=report("plots/allele-freqs.svg", caption="../report/freqs.rst", category="Plots")
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot-depths.py"
