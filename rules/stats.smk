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


rule update_tsv:
    input:
        "tables/calls.tsv.gz"
    output:
        "tables/calls_update.tsv.gz"
    run:
        sampToColsDict = {}
        annotations = ['GT', 'DP', 'AD']
        with gzip.open(input[0]) as f, gzip.open(output[0], 'wb') as out:
            head = f.readline()
            head_list = head.decode().split()[1:]
            new_head1 = []
            new_head2 = []
            for i in range(len(head_list)):
                field = head_list[i]
                fields = field.split(']')[1].split(':')
                if len(fields) == 1:
                    new_head1.append('VARIANT')
                    new_head2.append(fields[0])
                else:
                    samp = fields[0]
                    if not sampToColsDict.get(samp):
                        sampToColsDict[samp] = {}
                    sampToColsDict[samp][fields[1]] = i
            for samp in sorted(sampToColsDict.keys()):
                for i in range(len(annotations)):
                    new_head1.append(samp)
                    new_head2.append(annotations[i])
            out.write(('\t'.join(new_head1) + '\n').encode('utf-8'))
            out.write(('\t'.join(new_head2) + '\n').encode('utf-8'))
            line = f.readline()
            while line.decode() != '':
                line_list = line.decode().split()
                new_line_list = line_list[:6]
                for samp in sorted(sampToColsDict.keys()):
                    for annotation in annotations:
                        annotCol = sampToColsDict[samp][annotation]
                        new_line_list.append(line_list[annotCol])
                out.write(('\t'.join(new_line_list) + '\n').encode('utf-8'))
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
