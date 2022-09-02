
import pysam
import argparse

parser = argparse.ArgumentParser(description='To calculate the Read-through index of pNET-seq/RNA-seq data.',
                                epilog='Author: Fengli Zhao\tzhaofl@sustech.edu.cn      20220114')

parser.add_argument('--gene', type=str, help='File of gene list to be analyzed. [Non]', default='Non')
parser.add_argument('--interval', type=str, required=True, help='The gene interval file of genes.')
parser.add_argument('--tpm', type=str, required=True, help='TPM file of genes.')
parser.add_argument('--bam', type=str, required=True, help='The single-base BAM file.')

parser.add_argument('--Clen', type=int, help='The gene-body region to be used were ??bp Upstream to TES. [150]', default=150)
parser.add_argument('--Dlen', type=int, help='The read-through region were ??bp Downstream from TES. [300]', default=300)
parser.add_argument('--RTlen', type=int, help='The length of read-through region to be used (At most, it extends to the TSS of downstream gene). [500]', default=500)

parser.add_argument('--prefix', type=str, required=True, help='The prefix of output files.')

args = parser.parse_args()

gene = args.gene
interval = args.interval
tpm = args.tpm
bam = args.bam
Clen = args.Clen
Dlen = args.Dlen
RTlen = args.RTlen
prefix = args.prefix

#########################################################################################################

def readcount_sb (bam,chrom,strand0,start0,end0):
    read_num = 0; read_strand = ''
    INbam = pysam.AlignmentFile(bam, 'rb', threads=4)
    for read in INbam.fetch(chrom, start0, end0):
        if read.is_reverse:
            read_strand = '-'
        else:
            read_strand = '+'
        if read_strand == strand0:
            read_num += 1
    INbam.close()
    return read_num

#########################################################################################################

Genes = []
if gene != 'Non':
    INgene = open(gene, 'r')
    while True:
        line = INgene.readline()
        if not line:
            break
        Genes.append(line.strip())
    INgene.close()

Strand = {}; Interval = {}
INfr = open(interval, 'r')
while True:
    line = INfr.readline()
    if not line:
        break
    array = line.strip().split()
    Strand[array[0]] = array[1]
    Interval[array[0]] = int(array[2])
INfr.close()


INtpm = open(tpm, 'r')
OUTput = open(prefix + '.ReadThroughIndex', 'w')
jj = 0
while True:
    line = INtpm.readline()
    if not line:
        break
    array = line.strip().split()
    if jj == 0:
        OUTput.write(line.strip() + '\tStrand\tGeneBody_stat_posi\tLength\tRead\tReadThrough_region\tRT_length\tRT_read\tRT_index\n')
    else:
        if float(array[-1]) < 1:
            continue
        elif len(Genes) > 0 and array[0] not in Genes:
            continue
        elif array[0] in Strand.keys() and Interval[array[0]] > 400: # min Interval, 400bp
            beg_g = 0; end_g = 0; beg_rt = 0; end_rt = 0; rt_len = RTlen
            if Interval[array[0]] < 500:
                rt_len = Interval[array[0]] - 300
            if Strand[array[0]] == '+':
                beg_rt = int(array[3]) + Dlen + 1; end_rt = beg_rt + rt_len - 1
                beg_g = int(array[2]); end_g = int(array[3]) - Clen
            else:
                beg_rt = int(array[2]) - Dlen - rt_len + 1; end_rt = int(array[2]) - Dlen
                beg_g = int(array[2]) + Clen - 1; end_g = int(array[3])
                if beg_rt < 1:
                    beg_rt = 1

            len_g = end_g - beg_g + 1
            rt_len = end_rt - beg_rt + 1
            if len_g < 100: # Min-len of gene-body to be counted, 100 bp
                continue
            else:
                read_g = readcount_sb(bam, array[1], Strand[array[0]], beg_g - 1, end_g)
                if read_g == 0:
                    continue
                else:
                    read_rt = readcount_sb(bam, array[1], Strand[array[0]], beg_rt - 1, end_rt)
                    RT_index = format((read_rt / rt_len) / (read_g / len_g), '.4f')

                    OUTput.write(line.strip() + '\t' + Strand[array[0]] + '\t' + str(beg_g) + '-' + str(end_g) + '\t' + str(len_g) + '\t' + str(read_g) + '\t' \
                                 + str(beg_rt) + '-' + str(end_rt) + '\t' + str(rt_len) + '\t' + str(read_rt) + '\t' + RT_index + '\n')
    jj += 1
INtpm.close()
OUTput.close()

