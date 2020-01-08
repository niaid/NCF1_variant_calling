import pysam

def get_rg(bam):
    samfile = pysam.AlignmentFile(bam, "rb")
    read_groups = samfile.header['RG']
    samfile.close()
    return read_groups


def outputUnitsTsv():
    with open('bams.tsv') as f, open('units.tsv', 'w') as output:
        output.write('sample\tunit\tplatform\tfq1\tfq2\n')
        head = f.readline()
        line = f.readline()
        while line != '':
            (samp, bam) = line.split()
            read_groups = get_rg(bam)
            for i in range(len(read_groups)):
                rg = read_groups[i]
                unit = str(i + 1)
                platform = rg['PL']
                pu = rg['PU']
                fq1 = pu + '_1.fastq'
                fq2 = pu + '_2.fastq'
                output.write('\t'.join([samp, unit, platform, fq1, fq2]) + '\n')
            line = f.readline()

def main():
    args = sys.argv[1:]
    if len(args) != 1:
        print('error: usage output_units.py /path/to/file.bam')
        sys.exit(1)
    else:
        outputUnitsTsv(args[0])
        


if __name__ == "__main__":
    main()
