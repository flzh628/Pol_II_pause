
import os
import numpy as np
import pandas as pd
import scipy.stats as stats
import MultipleComparisons as MultiComp # 20220328
from plotnine import *
import argparse

parser = argparse.ArgumentParser(description='To analyze the relationship between PauseIndex and IntronNum.',
                                 epilog='Author: Fengli Zhao\tzhaofl@sustech.edu.cn      20220427 20220328 20211226 20211111 20210908')

parser.add_argument('--pcg', type=str, required=True, help='The PCG list, or gene list to analyze.')
parser.add_argument('--intronNum', type=str, help='The intron number file. [TAIR10.LongRNA.intron_num]', default='TAIR10.LongRNA.intron_num')
parser.add_argument('--pause', type=str, required=True, help='The pause index result.')
parser.add_argument('--prefix', type=str, required=True, help='The prefix of output files.')

args = parser.parse_args()

pcg = args.pcg
intronNum = args.intronNum
pause = args.pause
prefix = args.prefix

PCG = []
INpcg = open(pcg, 'r')
while True:
    line = INpcg.readline()
    if not line:
        break
    PCG.append(line.strip())
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

Genes = []; Index = []; LogIdx = []; IntronNum =  []; iNumber = []
INidx = open(pause, 'r')
while True:
    line = INidx.readline()
    if not line:
        break
    array = line.strip().split()
    if line[:4] != 'Gene':
        if array[0] in Number.keys() and array[0] in PCG and float(array[12]) > 0:
            Genes.append(array[0])
            Index.append(float(array[12]))
            LogIdx.append(np.log2(float(array[12])))
            if int(Number[array[0]]) <= 9:
                IntronNum.append('Intron_0' + Number[array[0]])
                iNumber.append(int(Number[array[0]]))
            else:
                IntronNum.append('Intron_m10')
                iNumber.append(int(Number[array[0]]))
INidx.close()

Group = sorted(list(set(IntronNum)))
df = pd.DataFrame({"Genes":Genes, "Index":Index, "LogIdx":LogIdx, "Intron_Num":IntronNum, "iNumber":iNumber})

df.to_csv(prefix + '.PauseIndex_Intron-number.data', sep='\t', index=False)

outstat = prefix + '.PauseIndex_Intron-number.statistics'
if os.path.exists(outstat):
    os.remove(outstat)
MultiComp.OneFactor_multiComparison(df, 'Intron_Num', 'Index', 'conditions', outstat)

breaks_arr = sorted(list(set(IntronNum)))
labels_arr = [x.replace('Intron_0', '') for x in breaks_arr[:-1]]
labels_arr.append(breaks_arr[-1].replace('Intron_', ''))
Count = []; yPos = []
for ii in range(len(breaks_arr)):
    yPos.append(min(df['LogIdx'])-0.5)
    df_g = df[df['Intron_Num'] == breaks_arr[ii]]
    Count.append(format(len(df_g['Index']), ','))

df_t = pd.DataFrame({"xPos":breaks_arr, "yPos":yPos, "Text":Count})

OUTstat = open(outstat, 'a')
OUTstat.write('\n\n' + '=' * 100 + '\n\n')

spearman_r0, p_val0 = stats.spearmanr(df['iNumber'], df['Index'])

OUTstat.write('The correlation between PauseIndex and Intron_num in ' + prefix + ' :\n\n')
OUTstat.write('\tspearman_r : ' + str(format(spearman_r0, '.2f')) + ' ( ' + str(spearman_r0) + ' )\n\tP_val : ' \
              + str(format(p_val0, '.2e')) + ' ( ' + str(p_val0) + ' )\n\n\n')
OUTstat.close()

xPos_r = []; yPos_r = []; Corr = [] # 20220427
sign = '' # 20220427
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
df_r = pd.DataFrame({"xPos":xPos_r, "yPos":yPos_r, "Text":Corr}) # 20220427

boxplot_fig = (ggplot(df, aes(x='Intron_Num', y='LogIdx'))
+ geom_boxplot(df, aes(fill='Intron_Num'), outlier_alpha=0, outlier_size=0, position=position_dodge(0.85))
+ ggtitle('The relationship between PauseIndex and Intron_num\nin ' + prefix)
+ labs(x="Intron number", y="Log2(Pause_index)")
+ scale_x_discrete(breaks=breaks_arr, labels=labels_arr)
+ geom_text(df_t, aes(x='xPos', y='yPos', label='Text'), family='arial', size=10)
+ geom_text(df_r, aes(x='xPos', y='yPos', label='Text'), family='arial', size=10) # 20220427
+ theme(plot_title=element_text(family="arial", size=14), axis_title=element_text(family="arial", size=12), \
        axis_text=element_text(family="arial", size=10), legend_position="none"))

boxplot_fig.save(prefix + '.PauseIndex_Intron-number.boxplot.pdf', format='pdf')

