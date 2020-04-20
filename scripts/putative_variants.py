import pysam

def get_ref_base(chrom, position, pysamFa):
    for base in pysamFa.fetch(chrom, position - 1, position):
        return base
    
def get_del_ref(chrom, position, length, pysamFa):
    ref = ''
    for base in pysamFa.fetch(chrom, position - 1, position + length):
        ref += base
    return ref


def makeBaseCountDict(samfile, chrom, testPos, ref = '/hpcdata/bio_data/current/GATK/bundle/2.8/b37/BWAIndex/human_g1k_v37_decoy.fasta'):
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


def get_variant_sites(samfile, chrom, start, end, min_alt):
    site_dict = {}
    for pos in range(start, end + 1):
        (baseCountDict, ref_base, indel_dict) = makeBaseCountDict(samfile, chrom, pos)
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


vcf_header = '##fileformat=VCFv4.2\n##FILTER=<ID=LowQual,Description="Low quality">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">\n##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller  --sample-ploidy 4 --output OtA2633.vcf --intervals /hpcdata/bcbb/karlinser/Projects/NEMO/IKBKG.bed --input /hpcdata/bcbb/karlinser/Projects/NEMO/2019-10-22_WES_MergeSam/OtA2633.bam --reference /hpcdata/bio_data/current/GATK/bundle/2.8/b37/human_g1k_v37.fasta  --emit-ref-confidence NONE --gvcf-gq-bands 1 --gvcf-gq-bands 2 --gvcf-gq-bands 3 --gvcf-gq-bands 4 --gvcf-gq-bands 5 --gvcf-gq-bands 6 --gvcf-gq-bands 7 --gvcf-gq-bands 8 --gvcf-gq-bands 9 --gvcf-gq-bands 10 --gvcf-gq-bands 11 --gvcf-gq-bands 12 --gvcf-gq-bands 13 --gvcf-gq-bands 14 --gvcf-gq-bands 15 --gvcf-gq-bands 16 --gvcf-gq-bands 17 --gvcf-gq-bands 18 --gvcf-gq-bands 19 --gvcf-gq-bands 20 --gvcf-gq-bands 21 --gvcf-gq-bands 22 --gvcf-gq-bands 23 --gvcf-gq-bands 24 --gvcf-gq-bands 25 --gvcf-gq-bands 26 --gvcf-gq-bands 27 --gvcf-gq-bands 28 --gvcf-gq-bands 29 --gvcf-gq-bands 30 --gvcf-gq-bands 31 --gvcf-gq-bands 32 --gvcf-gq-bands 33 --gvcf-gq-bands 34 --gvcf-gq-bands 35 --gvcf-gq-bands 36 --gvcf-gq-bands 37 --gvcf-gq-bands 38 --gvcf-gq-bands 39 --gvcf-gq-bands 40 --gvcf-gq-bands 41 --gvcf-gq-bands 42 --gvcf-gq-bands 43 --gvcf-gq-bands 44 --gvcf-gq-bands 45 --gvcf-gq-bands 46 --gvcf-gq-bands 47 --gvcf-gq-bands 48 --gvcf-gq-bands 49 --gvcf-gq-bands 50 --gvcf-gq-bands 51 --gvcf-gq-bands 52 --gvcf-gq-bands 53 --gvcf-gq-bands 54 --gvcf-gq-bands 55 --gvcf-gq-bands 56 --gvcf-gq-bands 57 --gvcf-gq-bands 58 --gvcf-gq-bands 59 --gvcf-gq-bands 60 --gvcf-gq-bands 70 --gvcf-gq-bands 80 --gvcf-gq-bands 90 --gvcf-gq-bands 99 --indel-size-to-eliminate-in-ref-model 10 --use-alleles-trigger false --disable-optimizations false --just-determine-active-regions false --dont-genotype false --max-mnp-distance 0 --dont-trim-active-regions false --max-disc-ar-extension 25 --max-gga-ar-extension 300 --padding-around-indels 150 --padding-around-snps 20 --kmer-size 10 --kmer-size 25 --dont-increase-kmer-sizes-for-cycles false --allow-non-unique-kmers-in-ref false --num-pruning-samples 1 --recover-dangling-heads false --do-not-recover-dangling-branches false --min-dangling-branch-length 4 --consensus false --max-num-haplotypes-in-population 128 --error-correct-kmers false --min-pruning 2 --debug-graph-transformations false --kmer-length-for-read-error-correction 25 --min-observations-for-kmer-to-be-solid 20 --likelihood-calculation-engine PairHMM --base-quality-score-threshold 18 --pair-hmm-gap-continuation-penalty 10 --pair-hmm-implementation FASTEST_AVAILABLE --pcr-indel-model CONSERVATIVE --phred-scaled-global-read-mismapping-rate 45 --native-pair-hmm-threads 4 --native-pair-hmm-use-double-precision false --debug false --use-filtered-reads-for-annotations false --bam-writer-type CALLED_HAPLOTYPES --dont-use-soft-clipped-bases false --capture-assembly-failure-bam false --error-correct-reads false --do-not-run-physical-phasing false --min-base-quality-score 10 --smith-waterman JAVA --use-new-qual-calculator false --annotate-with-num-discovered-alleles false --heterozygosity 0.001 --indel-heterozygosity 1.25E-4 --heterozygosity-stdev 0.01 --standard-min-confidence-threshold-for-calling 10.0 --max-alternate-alleles 6 --max-genotype-count 1024 --num-reference-samples-if-no-call 0 --genotyping-mode DISCOVERY --genotype-filtered-alleles false --contamination-fraction-to-filter 0.0 --output-mode EMIT_VARIANTS_ONLY --all-site-pls false --min-assembly-region-size 50 --max-assembly-region-size 300 --assembly-region-padding 100 --max-reads-per-alignment-start 50 --active-probability-threshold 0.002 --max-prob-propagation-distance 50 --interval-set-rule UNION --interval-padding 0 --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true --create-output-variant-md5 false --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --sites-only-vcf-output false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --disable-tool-default-read-filters false --minimum-mapping-quality 20 --disable-tool-default-annotations false --enable-all-annotations false",Version=4.0.8.1,Date="October 30, 2019 10:07:54 AM EDT">\n##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">\n##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">\n##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">\n##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">\n##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">\n##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">\n##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test to detect strand bias">\n##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">\n##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">\n##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">\n##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">\n##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">\n##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">\n##INFO=<ID=RAW_MQ,Number=1,Type=Float,Description="Raw data for RMS Mapping Quality">\n##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">\n##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">\n##contig=<ID=1,length=249250621>\n##contig=<ID=2,length=243199373>\n##contig=<ID=3,length=198022430>\n##contig=<ID=4,length=191154276>\n##contig=<ID=5,length=180915260>\n##contig=<ID=6,length=171115067>\n##contig=<ID=7,length=159138663>\n##contig=<ID=8,length=146364022>\n##contig=<ID=9,length=141213431>\n##contig=<ID=10,length=135534747>\n##contig=<ID=11,length=135006516>\n##contig=<ID=12,length=133851895>\n##contig=<ID=13,length=115169878>\n##contig=<ID=14,length=107349540>\n##contig=<ID=15,length=102531392>\n##contig=<ID=16,length=90354753>\n##contig=<ID=17,length=81195210>\n##contig=<ID=18,length=78077248>\n##contig=<ID=19,length=59128983>\n##contig=<ID=20,length=63025520>\n##contig=<ID=21,length=48129895>\n##contig=<ID=22,length=51304566>\n##contig=<ID=X,length=155270560>\n##contig=<ID=Y,length=59373566>\n##contig=<ID=MT,length=16569>\n##contig=<ID=GL000207.1,length=4262>\n##contig=<ID=GL000226.1,length=15008>\n##contig=<ID=GL000229.1,length=19913>\n##contig=<ID=GL000231.1,length=27386>\n##contig=<ID=GL000210.1,length=27682>\n##contig=<ID=GL000239.1,length=33824>\n##contig=<ID=GL000235.1,length=34474>\n##contig=<ID=GL000201.1,length=36148>\n##contig=<ID=GL000247.1,length=36422>\n##contig=<ID=GL000245.1,length=36651>\n##contig=<ID=GL000197.1,length=37175>\n##contig=<ID=GL000203.1,length=37498>\n##contig=<ID=GL000246.1,length=38154>\n##contig=<ID=GL000249.1,length=38502>\n##contig=<ID=GL000196.1,length=38914>\n##contig=<ID=GL000248.1,length=39786>\n##contig=<ID=GL000244.1,length=39929>\n##contig=<ID=GL000238.1,length=39939>\n##contig=<ID=GL000202.1,length=40103>\n##contig=<ID=GL000234.1,length=40531>\n##contig=<ID=GL000232.1,length=40652>\n##contig=<ID=GL000206.1,length=41001>\n##contig=<ID=GL000240.1,length=41933>\n##contig=<ID=GL000236.1,length=41934>\n##contig=<ID=GL000241.1,length=42152>\n##contig=<ID=GL000243.1,length=43341>\n##contig=<ID=GL000242.1,length=43523>\n##contig=<ID=GL000230.1,length=43691>\n##contig=<ID=GL000237.1,length=45867>\n##contig=<ID=GL000233.1,length=45941>\n##contig=<ID=GL000204.1,length=81310>\n##contig=<ID=GL000198.1,length=90085>\n##contig=<ID=GL000208.1,length=92689>\n##contig=<ID=GL000191.1,length=106433>\n##contig=<ID=GL000227.1,length=128374>\n##contig=<ID=GL000228.1,length=129120>\n##contig=<ID=GL000214.1,length=137718>\n##contig=<ID=GL000221.1,length=155397>\n##contig=<ID=GL000209.1,length=159169>\n##contig=<ID=GL000218.1,length=161147>\n##contig=<ID=GL000220.1,length=161802>\n##contig=<ID=GL000213.1,length=164239>\n##contig=<ID=GL000211.1,length=166566>\n##contig=<ID=GL000199.1,length=169874>\n##contig=<ID=GL000217.1,length=172149>\n##contig=<ID=GL000216.1,length=172294>\n##contig=<ID=GL000215.1,length=172545>\n##contig=<ID=GL000205.1,length=174588>\n##contig=<ID=GL000219.1,length=179198>\n##contig=<ID=GL000224.1,length=179693>\n##contig=<ID=GL000223.1,length=180455>\n##contig=<ID=GL000195.1,length=182896>\n##contig=<ID=GL000212.1,length=186858>\n##contig=<ID=GL000222.1,length=186861>\n##contig=<ID=GL000200.1,length=187035>\n##contig=<ID=GL000193.1,length=189789>\n##contig=<ID=GL000194.1,length=191469>\n##contig=<ID=GL000225.1,length=211173>\n##contig=<ID=GL000192.1,length=547496>\n##contig=<ID=NC_007605,length=171823>\n##contig=<ID=hs37d5,length=35477943>\n##source=HaplotypeCaller\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'


def output_vcf(vcf_header, out_vcf, samfile, chrom, start, end, min_alt):
    base_list = ['A', 'C', 'G', 'T', 'indel']
    with open(out_vcf, 'w') as output:
        output.write(vcf_header)
        variant_dict = get_variant_sites(samfile, chrom, start, end, min_alt)
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

