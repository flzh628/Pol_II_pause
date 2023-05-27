#!/public/home/Wuz/software/anaconda3/bin/python3

import os
import math
import pandas as pd
import scipy.stats as stats
import MultipleComparisons as MultiComp
from plotnine import *
import argparse

parser = argparse.ArgumentParser(description='To analyze the relationship between pauseIndex and IntronNum, at different SS-ratio levels.',
                                 epilog='Author: Fengli Zhao\tzhaofl@sustech.edu.cn      20220426 20220425 20211111 20210912 20210908')

parser.add_argument('--gpart', type=str, required=True, help='The SSratio-part of gene.')
parser.add_argument('--intronNum', type=str, help='The intron number file. [TAIR10.LongRNA.intron_num]', default='TAIR10.LongRNA.intron_num')
parser.add_argument('--pause', type=str, required=True, help='The pause index result.')
parser.add_argument('--prefix', type=str, required=True, help='The prefix of output files.')

args = parser.parse_args()

gpart = args.gpart
intronNum = args.intronNum
pause = args.pause
prefix = args.prefix

Parts = {}
INpart = open(gpart, 'r')
while True:
    line = INpart.readline()
    if not line:
        break
    array = line.strip().split('\t')
    Parts[array[0]] = array[1]
INpart.close()

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

Genes = []; Index = []; LogIdx = []; IntronNum =  []; iNumber = []; Part = []
INidx = open(pause, 'r')
while True:
    line = INidx.readline()
    if not line:
        break
    array = line.strip().split()
    if line[:4] != 'Gene':
        if array[0] in Number.keys() and array[0] in Parts.keys():
            Genes.append(array[0])
            Index.append(float(array[12]))
            LogIdx.append(math.log2(float(array[12])))
            Part.append(Parts[array[0]])
            if int(Number[array[0]]) <= 9:
                IntronNum.append('Intron_0' + Number[array[0]])
                iNumber.append(int(Number[array[0]]))
            else:
                IntronNum.append('Intron_m10')
                iNumber.append(int(Number[array[0]]))
INidx.close()

parts = sorted(list(set(Part)))
df = pd.DataFrame({"Genes":Genes, "Index":Index, "LogIdx":LogIdx, "Intron_Num":IntronNum, "iNumber":iNumber, "Part":Part})
df.to_csv(prefix + '.PauseIndex_IntronNum_SSratio-parts.data', sep='\t', index=False)

outstat = prefix + '.PauseIndex_IntronNum_SSratio-parts.statistics'
if os.path.exists(outstat):
    os.remove(outstat)

CORR = open(prefix + '.PauseIndex_IntronNum_SSratio-parts.corr', 'w')

breaks_arr = sorted(list(set(IntronNum)))
labels_arr = [x.replace('Intron_0', '') for x in breaks_arr[:-1]]
labels_arr.append(breaks_arr[-1].replace('Intron_', ''))

xPos = []; yPos = []; Count = []; cPart = [] # 20220425
xPos_r = []; yPos_r = []; Corr = []; rPart = [] # 20220426
for part in parts:
    df_p = df[df['Part'] == part]

    MultiComp.OneFactor_multiComparison(df_p, 'Intron_Num', 'Index', 'Fix_Gene-SSratio_' + part, outstat)

    spearman_r0, p_val0 = stats.spearmanr(df_p['iNumber'], df_p['Index'])

    CORR.write('The correlation between PauseIndex and Intron_num in ' + prefix + '.' + part + ' :\n\n')
    CORR.write('\tspearman_r : ' + format(spearman_r0, '.2f') + ' ( ' + str(spearman_r0) + ' )\n\tP_val : ' \
               + format(p_val0, '.2e') + ' ( ' + str(p_val0) + ' )\n\n\n')
    xPos_r.append(breaks_arr[int(0.75 * len(breaks_arr)) - 1])
    yPos_r.append(0.9 * max(df['LogIdx']))
    sign = ''
    if p_val0 >= 0.05:
        sign = 'ns'
    elif p_val0 >= 0.01:
        sign = '*'
    elif p_val0 >= 0.001:
        sign = '**'
    else:
        sign = '***'
    Corr.append('spearman R = ' + format(spearman_r0, '.2f') + ' (' + sign + ')')
    rPart.append(part)

    for str0 in breaks_arr: # 20220425
        xPos.append(str0)
        yPos.append(min(df_p['LogIdx'])-0.5)
        df_g = df_p[df_p['Intron_Num'] == str0]
        Count.append(format(len(df_g['Index']), ','))
        cPart.append(part)

CORR.close()

df_t = pd.DataFrame({"xPos":xPos, "yPos":yPos, "Text":Count, "Part":cPart}) # 20220425
df_r = pd.DataFrame({"xPos":xPos_r, "yPos":yPos_r, "Text":Corr, "Part":rPart}) # 20220426

boxplot_fig = (ggplot(df, aes(x='Intron_Num', y='LogIdx'))
+ geom_boxplot(df, aes(fill='Intron_Num'), outlier_alpha=0, outlier_size=0, position=position_dodge(0.85))
+ facet_grid('Part~.')
+ ggtitle('The relationship between PauseIndex and Intron_num\nin ' + prefix + ' at different SS-ratio levels')
+ labs(x="Intron number", y="Log2(Pause_index)")
+ scale_x_discrete(breaks=breaks_arr, labels=labels_arr)
+ geom_text(df_t, aes(x='xPos', y='yPos', label='Text'), family='arial', size=10) # 20220425
+ geom_text(df_r, aes(x='xPos', y='yPos', label='Text'), family='arial', size=10) # 20220426
+ theme_bw()
+ theme(plot_title=element_text(family="arial", size=14), axis_title=element_text(family="arial", size=12), axis_text=element_text(family="arial", size=10), \
        strip_text=element_text(family="arial", size=12), legend_position="none"))
#+ theme(plot_title=element_text(family="arial", size=14), axis_title=element_text(family="arial", size=12), axis_text=element_text(family="arial", size=10), \
#        legend_title= element_text(family="arial", size=12), legend_text= element_text(family="arial", size=10)))

boxplot_fig.save(prefix + '.PauseIndex_IntronNum_SSratio-parts.boxplot.pdf', width=4.8, height=9, format='pdf')
