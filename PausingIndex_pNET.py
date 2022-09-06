
import os
import pysam
import numpy as np
import scipy.stats as stats
from PvalueAdjust import BH_qvalues
import argparse

parser = argparse.ArgumentParser(description='To calculate the pausing index of pNET-seq data.',
                                epilog='Author: Fengli Zhao\tzhaofl@sustech.edu.cn      20220424 20220401 20220113 20211112 20210908 20210901')

parser.add_argument('--gene', type=str, help='File of gene list to be analyzed. [Non]', default='Non')
parser.add_argument('--f_r', type=str, required=True, help='The strand file of genes.')
parser.add_argument('--GRstr', type=str, help='Whether the read and gene are the same strand or not. [Yes]', default='Yes', choices=['Yes', 'Not']) # 20220424
parser.add_argument('--chrlen', type=str, required=True, help='The Chr-len file, without header.')
parser.add_argument('--tpm', type=str, required=True, help='TPM file of genes.')
parser.add_argument('--bam', type=str, required=True, help='The single-base BAM file.')
parser.add_argument('--pause', type=str, required=True, help='Which end does Pol-II pause.', choices=['TSS', 'TES', 'TSS-TES'])

parser.add_argument('--Alen', type=int, help='Upstream extension length to TSS. [150]', default=150)
parser.add_argument('--Blen', type=int, help='Downstream extension length from TSS. [150]', default=150)
parser.add_argument('--Clen', type=int, help='Upstream extension length to TES. [150]', default=150)
parser.add_argument('--Dlen', type=int, help='Downstream extension length from TES. [150]', default=150)

parser.add_argument('--win', type=int, help='Window size. [50]', default=50)
parser.add_argument('--step', type=int, help='Step size. [5]', default=5)
parser.add_argument('--flt', type=str, help='The threshold of p_val or q_val to filter. [q:0.05]', default='q:0.05')
parser.add_argument('--prefix', type=str, required=True, help='The prefix of output files.')

args = parser.parse_args()

gene = args.gene
f_r = args.f_r
GRstr = args.GRstr
chrlen = args.chrlen
tpm = args.tpm
bam = args.bam
pause = args.pause
Alen = args.Alen; Blen = args.Blen
Clen = args.Clen; Dlen = args.Dlen
win = args.win
step = args.step
flt = args.flt
prefix = args.prefix

#########################################################################################################

def readcount_sb (bam,chrom,start0,end0,strand0): # 20220401: strand-specific
    openbam = pysam.AlignmentFile(bam, mode='rb', threads=4)
    strand_r = ''; num_r = 0
    for read in openbam.fetch(chrom, start0, end0):
        if read.is_reverse:
            strand_r = '-'
        else:
            strand_r = '+'
        if strand_r == strand0 and GRstr == 'Yes': # 20220424
            num_r += 1
        elif strand_r != strand0 and GRstr == 'Not': # 20220424
            num_r += 1
    openbam.close()
    return num_r

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

ChrLen = {} # 20211111
INlen = open(chrlen, 'r')
while True:
    line = INlen.readline()
    if not line:
        break
    array = line.strip().split()
    ChrLen[array[0]] = int(array[1])
INlen.close()

Strand = {}
INfr = open(f_r, 'r')
while True:
    line = INfr.readline()
    if not line:
        break
    array = line.strip().split()
    Strand[array[0]] = array[1]
INfr.close()

Pause_arr = pause.split('-')
for pause0 in Pause_arr:

    jj = 0; P_vals = []
    INtpm = open(tpm, 'r')
    OUTput1 = open(prefix + '.' + pause0 + '.win.res', 'w')
    OUTput2 = open(prefix + '.' + pause0 + '.PauseIndex', 'w')
    while True:
        line = INtpm.readline()
        if not line:
            break
        array = line.strip().split()
        if jj == 0:
            OUTput1.write(line.strip() + '\tStrand\tGeneBody_stat_posi\tLength\tRead\tPause\tPause_posi\tWin_posi\tWin_read\n')
            OUTput2.write(line.strip() + '\tStrand\tGeneBody_stat_posi\tLength\tRead\tPause\tPause_win_posi\tPause_win_read\tPauseIndex\tp_value\n')
        else:
            if float(array[-1]) < 1:
                continue
            elif len(Genes) > 0 and array[0] not in Genes:
                continue
            elif array[0] in Strand.keys() and array[1] in ChrLen.keys(): # 20211112
                beg_g = 0; end_g = 0; beg_p = 0; end_p = 0
                if pause0 == 'TSS':
                    if Strand[array[0]] == '+':
                        beg_p = int(array[2]) - Alen + 1; end_p = int(array[2]) + Blen
                        beg_g = int(array[2]) + Blen + 150; end_g = int(array[3])
                    else:
                        beg_p = int(array[3]) - Blen + 1; end_p = int(array[3]) + Alen
                        beg_g = int(array[2]); end_g = int(array[3]) - Blen - 150
                else:
                    if Strand[array[0]] == '+':
                        beg_p = int(array[3]) - Clen + 1; end_p = int(array[3]) + Dlen
                        beg_g = int(array[2]); end_g = int(array[3]) - Clen - 150
                    else:
                        beg_p = int(array[2]) - Dlen + 1; end_p = int(array[2]) + Clen
                        beg_g = int(array[2]) + Clen + 150; end_g = int(array[3])
                if len(Pause_arr) == 2:
                    if pause0 == 'TSS':
                        if Strand[array[0]] == '+':
                            end_g = int(array[3]) - Clen - 150
                        else:
                            beg_g = int(array[2]) + Clen + 150
                    else:
                        if Strand[array[0]] == '+':
                            beg_g = int(array[2]) + Blen + 150
                        else:
                            end_g = int(array[3]) - Blen - 150
                len_g = end_g - beg_g + 1
                len_p = end_p - beg_p + 1
                if len_g < 100: # Min-len of gene-body to be counted, 100 bp
                    continue
                else:
                    read_g = readcount_sb(bam, array[1], beg_g - 1, end_g, Strand[array[0]])
                    if read_g == 0:
                        continue
                    else:
                        Beg_w_arr = []; Read_w_arr = []; idx = 0
                        win_num = int((len_p - win) / step)
                        for i in range(win_num):
                            win_beg = beg_p + i * step
                            if win_beg >= ChrLen[array[1]]: # 20211112
                                break
                            win_end = win_beg + win - 1
                            read_w = readcount_sb(bam, array[1], win_beg - 1, win_end, Strand[array[0]])
                            OUTput1.write(line.strip() + '\t' + Strand[array[0]] + '\t' + str(beg_g) + '-' + str(end_g) + '\t' + str(len_g) + '\t' + str(read_g) + '\t' \
                                          + pause0 + '\t' + str(beg_p) + '-' + str(end_p) + '\t' + str(win_beg) + '-' + str(win_end) + '\t' + str(read_w) + '\n')

                            Beg_w_arr.append(win_beg); Read_w_arr.append(read_w) # 20220113, BUG correct: Beg_w_arr.append(beg_p) -> Beg_w_arr.append(win_beg)

                        idx_arr = np.where(np.array(Read_w_arr)==max(Read_w_arr))[0].tolist()
                        if pause0 == 'TSS':
                            if Strand[array[0]] == '+':
                                idx = idx_arr[-1]
                            else:
                                idx = idx_arr[0]
                        else:
                            if Strand[array[0]] == '+':
                                idx = idx_arr[0]
                            else:
                                idx = idx_arr[-1]
                        read_p_max = Read_w_arr[idx]; beg_p_max = Beg_w_arr[idx]; end_p_max = beg_p_max + win - 1
                        P_index = format((read_p_max / win) / (read_g / len_g), '.2f')
                        p_val = stats.fisher_exact([[read_p_max, win], [read_g, len_g]], alternative = 'greater')[1]

                        OUTput2.write(line.strip() + '\t' + Strand[array[0]] + '\t' + str(beg_g) + '-' + str(end_g) + '\t' + str(len_g) + '\t' + str(read_g) + '\t' \
                                      + pause0 + '\t' + str(beg_p_max) + '-' + str(end_p_max) + '\t' + str(read_p_max) + '\t' + str(P_index) + '\t' + str(p_val) + '\n')
                        P_vals.append(str(p_val))
        jj += 1
    INtpm.close()
    OUTput1.close()
    OUTput2.close()

    p_q, pq_max = flt.split(':')

    Q_vals = BH_qvalues(P_vals, 'Not')

    INres = open(prefix + '.' + pause0 + '.PauseIndex', 'r')
    OUTres = open(prefix + '.' + pause0 + '.PauseIndex' + '.q', 'w')
    OUTres0 = open(prefix + '.' + pause0 + '.PauseIndex' + '.sign', 'w')

    ii = 0
    while True:
        line = INres.readline()
        if not line:
            break
        if ii == 0:
            OUTres.write(line.strip() + '\tq_value\n')
            OUTres0.write(line.strip() + '\tq_value\n')
        else:
            array = line.strip().split() # 20220702
            OUTres.write(line.strip() + '\t' + str(Q_vals[ii - 1]) + '\n')
            if p_q == 'q' and float(Q_vals[ii - 1]) <= float(pq_max):
                OUTres0.write(line.strip() + '\t' + str(Q_vals[ii - 1]) + '\n')
            elif p_q == 'p' and float(array[-1]) <= float(pq_max):
                OUTres0.write(line.strip() + '\t' + str(Q_vals[ii - 1]) + '\n')
        ii += 1

    INres.close()
    OUTres.close()
    OUTres0.close()

    os.remove(prefix + '.' + pause0 + '.PauseIndex')
    os.rename(prefix + '.' + pause0 + '.PauseIndex' + '.q', prefix + '.' + pause0 + '.PauseIndex')

