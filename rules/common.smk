#!/usr/bin/env python
import pandas as pd
import pysam
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.7.1")

report: "../report/workflow.rst"

calling_methods = ['ploidy', 'known']

###### Function to parse bam file #####
def get_rg_from_bam(bam):
    '''
    (str) -> [dict]
    returns a list of readgroup dicts for the given bam
    '''
    
    samfile = pysam.AlignmentFile(bam, "rb")
    read_groups = samfile.header['RG']
    samfile.close()
    return read_groups


###### Config file and sample sheets #####
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")
d = { 'sample':[],
  'bam': [],
  'unit':[],
  'ID': [],
  'LB':[],
  'platform':[],
  'PU': []}
for index, row in samples.iterrows():
    sample = row['sample']
    bam = row['bam']
    readgroups = get_rg_from_bam(bam)
    for i in range(len(readgroups)):
        d['sample'].append(sample)
        d['bam'].append(bam)
        d['unit'].append(str(i+1))
        d['ID'].append(readgroups[i]['ID'])
        d['LB'].append(readgroups[i]['LB'])
        d['PU'].append(readgroups[i]['PU'])
        d['platform'].append('ILLUMINA')
units = pd.DataFrame(data=d).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")


##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"])


##### Helper functions #####

def get_start_bam(wildcards):
    """get input bam given sample-units"""
    return samples.loc[(wildcards.sample, wildcards.unit), ["bam"]].dropna().bam

def get_rg_subset_param(wildcards):
    """get the param to use for samtools view, namely the RG string"""
    rg = units.loc[(wildcards.sample, wildcards.unit), ["ID"]].dropna().ID
    bed = config["processing"]["restrict-regions"]
    return "-hbr " + rg + "-L " + bed


def get_fai():
    return config["ref"]["genome"] + ".fai"


# contigs in reference genome
def get_contigs():
    return pd.read_table(get_fai(),
                         header=None, usecols=[0], squeeze=True, dtype=str)

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def is_single_end(sample, unit):
    """Return True if sample-unit is single end.
    I'll always expect paired-end sequencing and just return false here.
    """
    return False


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("recal/{sample}-{unit}.bam",
                  sample=wildcards.sample,
                  unit=units.loc[wildcards.sample].unit)


def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_call_variants_params(wildcards, input):
    return (get_regions_param(regions=input.regions, default="--intervals {}".format(wildcards.contig)) +
            config["params"]["gatk"]["HaplotypeCaller"])


def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "mapped/{sample}-{unit}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "dedup/{sample}-{unit}.bam"
    if bai:
        if config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f
