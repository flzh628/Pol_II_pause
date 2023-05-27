#!/public/home/Wuz/software/anaconda3/bin/python3

import os
import math
import pandas as pd
import scipy.stats as stats
import MultipleComparisons as MultiComp # 20220424
from plotnine import *
import argparse

parser = argparse.ArgumentParser(description='To analyze the relationship between PauseIndex, IntronNum, and GeneLen.',
                                 epilog='Author: Fengli Zhao\tzhaofl@sustech.edu.cn      20220426 20220425 20220424 20211005')

parser.add_argument('--pcg_len', type=str, required=True, help='The length file of PCG.')
parser.add_argument('--intronNum', type=str, help='The intron number file. [TAIR10.LongRNA.intron_num]', default='TAIR10.LongRNA.intron_num')
parser.add_argument('--pause', type=str, required=True, help='The pause index result.')
parser.add_argument('--prefix', type=str, required=True, help='The prefix of output files.')

args = parser.parse_args()

pcg_len = args.pcg_len
intronNum = args.intronNum
pause = args.pause
prefix = args.prefix

GeneLen = {}; ii = 0
INpcg = open(pcg_len, 'r')
while True:
    line = INpcg.readline()
    if not line:
        break
    array = line.strip().split()
    if ii > 0:
        deci1, int1 = math.modf(int(array[1]) / 1000)
        if float(format(deci1, '.3f')) <= 0.2:
            GeneLen[array[0]] = 'Len' + '0' * (2 - len(str(int(int1)))) + str(int(int1)) + 'K'
        elif float(format(deci1, '.3f')) >= 0.8:
            GeneLen[array[0]] = 'Len' + '0' * (2 - len(str(int(int1 + 1)))) + str(int(int1 + 1)) + 'K'
    ii += 1
INpcg.close()

Number = {}
INinf = open(intronNum, 'r')
while True:
    line = INinf.readline()
    if not line:
        break
    if line[:3] != 'LOC':
        array = line.strip().split()
        Number[array[0]] = array[1]
INinf.close()

Genes = []; Index = []; LogIdx = []; IntronNum =  []; iNumber = []; GLen = []; iGLen = []
INidx = open(pause, 'r')
while True:
    line = INidx.readline()
    if not line:
        break
    array = line.strip().split()
    if line[:4] != 'Gene':
        if array[0] in Number.keys() and array[0] in GeneLen.keys():
            Genes.append(array[0])
            Index.append(float(array[12]))
            LogIdx.append(math.log2(float(array[12])))
            GLen.append(GeneLen[array[0]])
            len0 = GeneLen[array[0]]
            len1 = len0.replace('Len', '')
            len2 = len1.replace('K', '')
            iGLen.append(int(len2))
            if int(Number[array[0]]) <= 9:
                IntronNum.append('Intron_0' + Number[array[0]])
                iNumber.append(int(Number[array[0]]))
            else:
                IntronNum.append('Intron_m10')
                iNumber.append(int(Number[array[0]]))
INidx.close()

Group = sorted(list(set(IntronNum)))
LenArr = sorted(list(set(GLen)))

df = pd.DataFrame({"Genes":Genes, "Index":Index, "LogIdx":LogIdx, "Intron_Num":IntronNum, "iNumber":iNumber, "GLen":GLen, "iGLen":iGLen})

df.to_csv(prefix + '.PauseIndex_IntronNum_GeneLen.data', sep='\t', index=False)

outstat = prefix + '.PauseIndex_IntronNum_GeneLen.statistics'
if os.path.exists(outstat):
    os.remove(outstat)

CORR = open(prefix + '.PauseIndex_Intron-number_GeneLen.corr', 'w')

for len0 in LenArr:
    df_l = df[df['GLen'] == len0]
    test_arr = list(set(df_l['Intron_Num']))
    if len(df_l) >= 10 and len(test_arr) >= 2:
        len_0 = len0.replace('Len', ''); len_0 = len_0.replace('K', '')
        breaks_arr = sorted(list(set(df_l['Intron_Num'])))
        labels_arr = [x.replace('Intron_0', '') for x in breaks_arr]
        labels_arr = [x.replace('Intron_', '') for x in labels_arr]

        Count = []; yPos = [] # 20220425
        for ii in range(len(breaks_arr)):
            yPos.append(min(df_l['LogIdx'])-0.5)
            df_g = df_l[df_l['Intron_Num'] == breaks_arr[ii]]
            Count.append(format(len(df_g['Index']), ','))

        df_t = pd.DataFrame({"xPos":breaks_arr, "yPos":yPos, "Text":Count})

        spearman_r0, p_val0 = stats.spearmanr(df_l['iNumber'], df_l['Index'])

        CORR.write('\nThe correlation between PauseIndex and Intron_num in ' + prefix + ', when gene_len is ' + len_0 + ' :\n\n')
        CORR.write('\tspearman_r : ' + str(format(spearman_r0, '.2f')) + ' ( ' + str(spearman_r0) + ' )\n\tP_value : ' \
                   + str(format(p_val0, '.2e')) + ' ( ' + str(p_val0) + ' )\n\n\n')

        xPos_r = []; yPos_r = []; Corr = [] # 20220426
        sign = '' # 20220426
        if p_val0 >= 0.05:
            sign = 'ns'
        elif p_val0 >= 0.01:
            sign = '*'
        elif p_val0 >= 0.001:
            sign = '**'
        else:
            sign = '***'
        Corr.append('spearman R = ' + format(spearman_r0, '.2f') + ' (' + sign + ')')
        xPos_r.append(breaks_arr[int(0.75 * len(breaks_arr)) - 1])
        yPos_r.append(0.9 * max(df['LogIdx']))
        df_r = pd.DataFrame({"xPos":xPos_r, "yPos":yPos_r, "Text":Corr}) # 20220426

        boxplot_fig = (ggplot(df_l, aes(x='Intron_Num', y='LogIdx'))
        + geom_boxplot(df_l, aes(fill='Intron_Num'), outlier_alpha=0, outlier_size=0, position=position_dodge(0.85))
        + ggtitle('The relationship between PauseIndex and Intron_num,\nwhen the gene_len is ' + str(int(len_0)) + ' +- 0.2 Kb')
        + labs(x="Intron number", y="Log2(Pause_index)")
        + scale_x_discrete(breaks=breaks_arr, labels=labels_arr)
        + geom_text(df_t, aes(x='xPos', y='yPos', label='Text'), family='arial', size=10) # 20220425 correct the bug: label='Count', even it doesn't affect the outcome
        + geom_text(df_r, aes(x='xPos', y='yPos', label='Text'), family='arial', size=10) # 20220426
        + theme(plot_title=element_text(family="arial", size=14), axis_title=element_text(family="arial", size=12), \
                axis_text=element_text(family="arial", size=10), legend_position="none"))
     
        boxplot_fig.save(prefix + '.Gene' + len0 + '.PauseIndex_IntronNum.boxplot.pdf', format='pdf')
    
        MultiComp.OneFactor_multiComparison(df_l, 'Intron_Num', 'Index', 'Fix_GeneLen_' + len0, outstat)
 
for num0 in Group:
    df_n = df[df['Intron_Num'] == num0]
    test_arr = list(set(df_n['GLen']))
    if len(df_n) >= 10 and len(test_arr) >= 2:
        num_0 = num0.replace('Intron_', '')
        if 'm' in num_0:
            num_0 = num_0.replace('m', '')
        breaks_arr = sorted(list(set(df_n['GLen'])))
        labels_arr = sorted(list(set(df_n['iGLen'])))

        Count = []; yPos = [] # 20220425
        for ii in range(len(breaks_arr)):
            yPos.append(min(df_n['LogIdx'])-0.5)
            df_g = df_n[df_n['GLen'] == breaks_arr[ii]]
            Count.append(format(len(df_g['Index']), ','))

        df_t = pd.DataFrame({"xPos":breaks_arr, "yPos":yPos, "Text":Count})

        spearman_r0, p_val0 = stats.spearmanr(df_n['iGLen'], df_n['Index'])

        CORR.write('The correlation between PauseIndex and GeneLen in ' + prefix + ', when IntronNum is ' + num_0 + ' :\n\n')
        CORR.write('\tspearman_r : ' + str(format(spearman_r0, '.2f')) + ' ( ' + str(spearman_r0) + ' )\n\tP_value : ' \
                   + str(format(p_val0, '.2e')) + ' ( ' + str(p_val0) + ' )\n\n\n')

        xPos_r = []; yPos_r = []; Corr = [] # 20220426
        sign = '' # 20220426
        if p_val0 >= 0.05:
            sign = 'ns'
        elif p_val0 >= 0.01:
            sign = '*'
        elif p_val0 >= 0.001:
            sign = '**'
        else:
            sign = '***'
        Corr.append('spearman R = ' + format(spearman_r0, '.2f') + ' (' + sign + ')')
        xPos_r.append(breaks_arr[int(0.75 * len(breaks_arr)) - 1])
        yPos_r.append(0.9 * max(df['LogIdx']))
        df_r = pd.DataFrame({"xPos":xPos_r, "yPos":yPos_r, "Text":Corr}) # 20220426

        boxplot_fig = (ggplot(df_n, aes(x='GLen', y='LogIdx'))
        + geom_boxplot(df_n, aes(fill='GLen'), outlier_alpha=0, outlier_size=0, position=position_dodge(0.85))
        + ggtitle('The relationship between PauseIndex and GeneLen,\nwhen the Intron_num is ' + str(int(num_0)))
        + labs(x="Gene length (Kb)", y="Log2(Pause_index)")
        + scale_x_discrete(breaks=breaks_arr, labels=labels_arr)
        + geom_text(df_t, aes(x='xPos', y='yPos', label='Text'), family='arial', size=10) # 20220425 correct the bug: label='Count', even it doesn't affect the outcome
        + geom_text(df_r, aes(x='xPos', y='yPos', label='Text'), family='arial', size=10) # 20220426
        + theme(plot_title=element_text(family="arial", size=14), axis_title=element_text(family="arial", size=12), \
                axis_text=element_text(family="arial", size=10), legend_position="none"))
     
        boxplot_fig.save(prefix + '.IntronNum' + num_0 + '.PauseIndex_GeneLen.boxplot.pdf', format='pdf')

        MultiComp.OneFactor_multiComparison(df_n, 'GLen', 'Index', 'Fix_IntronNum_' + num0, outstat)
     
CORR.close()

