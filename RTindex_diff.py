
import os
import scipy.stats as stats
import argparse

parser = argparse.ArgumentParser(description='To perform the differential analysis of PauseIndex between two samples.',
                                epilog='Author: Fengli Zhao\tzhaofl@sustech.edu.cn      20220116')

parser.add_argument('--col_rt', type=str, required=True, help='The Read-through index of Col.')
parser.add_argument('--mut_rt', type=str, required=True, help='The Read-through index of mutant.')
parser.add_argument('--col_tot', type=int, required=True, help='The total read number of Col.')
parser.add_argument('--mut_tot', type=int, required=True, help='The total read number of mutant.')
parser.add_argument('--flt', type=str, help='The threshold of p_val or q_val to filter. [q:0.05]', default='q:0.05')
parser.add_argument('--prefix', type=str, required=True, help='The prefix of output files.')

args = parser.parse_args()

col_rt = args.col_rt
mut_rt = args.mut_rt
col_tot = args.col_tot
mut_tot = args.mut_tot
flt = args.flt
prefix = args.prefix

#########################################################################################################

def BH_qvalues(p_values, sort_ornot):
    if p_values == []:
        return []
    Max_rank = {}; total_0 = len(p_values); Q_val = {}; q_values = []
    if sort_ornot == 'Not':
        openTmp = open('p-tmp', 'w')
        p_val_str0 = '\n'.join(p_values)
        openTmp.write(p_val_str0 + '\n')
        openTmp.close()
        os.system('msort -k n1 p-tmp >p-tmp.msort')
        os.remove('p-tmp')
        openTmpN = open('p-tmp.msort', 'r')
        ii_0 = 0
        while True:
            lineTmpN = openTmpN.readline()
            if not lineTmpN:
                break
            ii_0 += 1
            Max_rank[float(lineTmpN.strip())] = ii_0
        openTmpN.close()
        os.remove('p-tmp.msort')
    elif sort_ornot == 'Yes':
        for ii_1 in range(total_0):
            Max_rank[float(p_values[ii_1])] = ii_1 + 1
    else:
        os._exit()

    p_uniq_sort = sorted(Max_rank.keys())
    for p_val_0 in p_uniq_sort:
        rank = Max_rank[p_val_0]
        Q_val[p_val_0] = total_0 * float(p_val_0) / rank
    for i_0 in range(len(p_uniq_sort) - 2, -1, -1):
        if Q_val[p_uniq_sort[i_0]] > Q_val[p_uniq_sort[i_0 + 1]]:
            Q_val[p_uniq_sort[i_0]] = Q_val[p_uniq_sort[i_0 + 1]]
    for mm_0 in range(total_0):
        if float(p_values[mm_0]) in Q_val.keys():
            q_values.append(float(Q_val[float(p_values[mm_0])]))

    return q_values

#########################################################################################################

factor = round(mut_tot / col_tot, 2)

gbLen = {}; rtLen = {}
Col_gb_Num = {}; Col_win_Num = {}; Col_RTi = {}
INcol = open(col_rt, 'r')
while True:
    line = INcol.readline()
    if not line:
        break
    if line[:4] != 'Gene':
        array = line.strip().split()
        gbLen[array[0]] = int(array[7]); rtLen[array[0]] = int(array[10])
        Col_gb_Num[array[0]] = int(int(array[8]) * factor)
        Col_win_Num[array[0]] = int(int(array[11]) * factor)
        Col_RTi[array[0]] = float(array[12])
INcol.close()

Mut_gb_Num = {}; Mut_win_Num = {}; Mut_RTi = {}
INmut = open(mut_rt, 'r')
while True:
    line = INmut.readline()
    if not line:
        break
    if line[:4] != 'Gene':
        array = line.strip().split()
        gbLen[array[0]] = int(array[7]); rtLen[array[0]] = int(array[10])
        Mut_gb_Num[array[0]] = int(array[8])
        Mut_win_Num[array[0]] = int(array[11])
        Mut_RTi[array[0]] = float(array[12])
INmut.close()

P_vals = []
OUTput = open(prefix + '.RTindex_diff', 'w')
OUTput.write('Gene_Id\tGeneBody_statLen\tCol_gbNum\tCol_rtWin_Num\tCol_RTindex\tMut_gbNum\tMut_rtWin_Num\tMut_RTindex\tRTindex_FC\tUpDw\tp_value\n')
for gene in sorted(gbLen.keys()):
    fc = 0; p_val = 1
    if gene not in Col_RTi.keys():
        num1 = gbLen[gene]; num2 = 0; fc = round((Mut_RTi[gene] + 0.001) / 0.001, 4)
        p_val = stats.fisher_exact([[Mut_win_Num[gene], Mut_gb_Num[gene]], [num2, num1]], alternative = 'greater')[1]
        OUTput.write(gene + '\t' + str(gbLen[gene]) + '\t' + '--\t--\t--' + '\t' + str(Mut_gb_Num[gene]) + '\t' + str(Mut_win_Num[gene]) + '\t' \
                     + str(Mut_RTi[gene]) + '\t' + str(fc) + '\tUp\t' + str(p_val) + '\n')
    elif gene not in Mut_RTi.keys():
        num1 = gbLen[gene]; num2 = 0; fc = round(0.001 / (Col_RTi[gene] + 0.001), 4)
        p_val = stats.fisher_exact([[Col_win_Num[gene], Col_gb_Num[gene]], [num2, num1]], alternative = 'greater')[1]
        OUTput.write(gene + '\t' + str(gbLen[gene]) + '\t' + str(Col_gb_Num[gene]) + '\t' + str(Col_win_Num[gene]) + '\t' + str(Col_RTi[gene]) + '\t' \
                     + '--\t--\t--' + '\t' + str(fc) + '\tDw\t' + str(p_val) + '\n')
    else:
        OUTput.write(gene + '\t' + str(gbLen[gene]) + '\t' + str(Col_gb_Num[gene]) + '\t' + str(Col_win_Num[gene]) + '\t' + str(Col_RTi[gene]) + '\t' \
                     + str(Mut_gb_Num[gene]) + '\t' + str(Mut_win_Num[gene]) + '\t' + str(Mut_RTi[gene]))
        fc = round((Mut_RTi[gene] + 0.001) / (Col_RTi[gene] + 0.001), 4)
        if fc > 1:
            p_val = stats.fisher_exact([[Mut_win_Num[gene], Mut_gb_Num[gene]], [Col_win_Num[gene], Col_gb_Num[gene]]], alternative = 'greater')[1]
            OUTput.write('\t' + str(fc) + '\tUp\t' + str(p_val) + '\n')
        elif fc < 1:
            p_val = stats.fisher_exact([[Col_win_Num[gene], Col_gb_Num[gene]], [Mut_win_Num[gene], Mut_gb_Num[gene]]], alternative = 'greater')[1]
            OUTput.write('\t' + str(fc) + '\tDw\t' + str(p_val) + '\n')
        else:
            OUTput.write('\t' + str(fc) + '\tEq\t1.00\t' + '\n')
    P_vals.append(str(p_val))
OUTput.close()

p_q, pq_max = flt.split(':')

Q_vals = BH_qvalues(P_vals, 'Not')

INres = open(prefix + '.RTindex_diff', 'r')
OUTres = open(prefix + '.RTindex_diff' + '.q', 'w')
OUTres0 = open(prefix + '.RTindex_diff' + '.sign', 'w')

ii = 0
while True:
    line = INres.readline()
    if not line:
        break
    if ii == 0:
        OUTres.write(line.strip() + '\tq_value\n')
        OUTres0.write(line.strip() + '\tq_value\n')
    else:
        OUTres.write(line.strip() + '\t' + str(Q_vals[ii - 1]) + '\n')
        if p_q == 'q' and float(Q_vals[ii - 1]) <= float(pq_max):
            OUTres0.write(line.strip() + '\t' + str(Q_vals[ii - 1]) + '\n')
        elif p_q == 'p' and float(array[-2]) <= float(pq_max):
            OUTres0.write(line.strip() + '\t' + str(Q_vals[ii - 1]) + '\n')
    ii += 1

INres.close()
OUTres.close()
OUTres0.close()

os.remove(prefix + '.RTindex_diff')
os.rename(prefix + '.RTindex_diff' + '.q', prefix + '.RTindex_diff')
