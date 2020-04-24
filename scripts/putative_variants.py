#!/usr/bin/env python2.7

import pysam
import gzip
import sys
import argparse


def get_ref_base(chrom, position, pysamFa):
    for base in pysamFa.fetch(chrom, position - 1, position):
        return base
    
def get_del_ref(chrom, position, length, pysamFa):
    ref = ''
    for base in pysamFa.fetch(chrom, position - 1, position + length):
        ref += base
    return ref


def makeBaseCountDict(samfile, chrom, testPos, ref):
    '''
    (str, int) -> dict
    '''
    pysamFa = pysam.FastaFile(ref)
    ref_base = get_ref_base(chrom, testPos, pysamFa)
    baseCountDict = {'A':0, 'C':0, 'G':0, 'T':0, 'N': 0, 'indel':0}
    indel_dict = {}
    for pileupcolumn in samfile.pileup(chrom, testPos - 5, testPos + 5):
        if int(pileupcolumn.pos) == testPos - 1:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    baseCountDict[base] += 1
                if pileupread.indel < 0:
                    baseCountDict['indel'] += 1
                    if indel_dict.get(pileupread.indel):
                        (del_ref, del_alt, count) = indel_dict[pileupread.indel]
                        count += 1
                        indel_dict[pileupread.indel] = (del_ref, del_alt, count)
                    else:
                        del_ref = get_del_ref(chrom, testPos, abs(pileupread.indel), pysamFa)
                        del_alt = del_ref[0]
                        indel_dict[pileupread.indel] = (del_ref, del_alt, 1)
                elif pileupread.indel > 0:
                    baseCountDict['indel'] += 1
                    
                    if indel_dict.get(pileupread.indel):
                        (ins_ref, ins_alt, count) = indel_dict[pileupread.indel]
                        count += 1
                        indel_dict[pileupread.indel] = (ins_ref, ins_alt, count)
                    else:
                        insertion = ''
                        for i in range(pileupread.indel + 1):    
                            insertion += pileupread.alignment.query_sequence[pileupread.query_position + i]
                        ins_ref =  insertion[0]
                        indel_dict[pileupread.indel] = (ins_ref, insertion, 1)
    pysamFa.close()            
    return (baseCountDict, ref_base, indel_dict)


def get_variant_sites(samfile, ref, chrom, start, end, min_alt):
    site_dict = {}
    for pos in range(start, end + 1):
        (baseCountDict, ref_base, indel_dict) = makeBaseCountDict(samfile, chrom, pos, ref)
        c = 0
        for base in baseCountDict.keys():
            count = baseCountDict[base]
            if count >= min_alt:
                c += 1
        if c >= 2:
            site_dict[pos] = {}
            site_dict[pos]['ref'] = ref_base
            for base in baseCountDict.keys():
                count = baseCountDict[base]
                if count >= min_alt:
                    site_dict[pos][base] = count
            site_dict[pos]['indel_dict'] = indel_dict
            
    return site_dict


def output_vcf(vcf_header, out_vcf, samfile, ref, chrom, start, end, min_alt):
    base_list = ['A', 'C', 'G', 'T', 'indel']
    with open(out_vcf, 'w') as output:
        output.write(vcf_header)
        variant_dict = get_variant_sites(samfile, ref, chrom, start, end, min_alt)
        for pos in sorted(variant_dict.keys()):
            pos_dict = variant_dict[pos]
            ref_base = pos_dict['ref']
            alt_list = []
            for base in base_list:
                if base != ref_base and pos_dict.get(base):
                    alt_list.append(base)
            for alt in alt_list:
                if alt == 'indel':
                    indel_dict = pos_dict['indel_dict']
                    for key in indel_dict.keys():
                        (REF, ALT, count) = indel_dict[key]
                        output.write('\t'.join([chrom, str(pos), '.', REF, ALT, '.', '.', '.']) + '\n')
                else:
                    output.write('\t'.join([chrom, str(pos), '.', ref_base, alt, '.', '.', '.']) + '\n')

def get_vcf_header(vcf):
    vcf_header = ''
    if vcf.endswith('.gz'):
        f = gzip.open(vcf)
    else:
        f = open(vcf)
    line = f.readline()
    while not line.startswith('#CHROM') and line != '':
        vcf_header += line
        line = f.readline()
    vcf_header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    f.close()
    return vcf_header



def get_args():
    '''
    return the arguments from parser
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam', type=str, required=True, help='REQUIRED. Path to bam file')
    parser.add_argument('-o', '--out', type=str, required=True, help='REQUIRED. Path to output vcf file')
    parser.add_argument('-v', '--vcf', type=str, required=True, help='REQUIRED. Path to vcf file made from reference to use to get the vcf header')
    parser.add_argument('-R', '--reference', type=str, required=True, help='REQUIRED. Path to reference fasta file')
    parser.add_argument('-c', '--chrom', type=str, required=True, help='REQUIRED. Chromosome')
    parser.add_argument('-s', '--start', type=int, required=True, help='REQUIRED. Start position')
    parser.add_argument('-e', '--end', type=int, required=True, help='REQUIRED. End position')
    parser.add_argument('-m', '--min_alt', type=int, required=True, default=2, help='REQUIRED. Default = 2. Minimum number of ALT alleles to emit a variant position.')
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    samfile = pysam.AlignmentFile(args.bam, 'rb')
    vcf_header = get_vcf_header(args.vcf)
    output_vcf(vcf_header, args.out, samfile, args.reference, args.chrom, args.start, args.end, args.min_alt)


if __name__ == "__main__":
    main()
