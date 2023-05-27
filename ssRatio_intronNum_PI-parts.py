#!/public/home/Wuz/software/anaconda3/bin/python3

import os
import pandas as pd
import scipy.stats as stats
import MultipleComparisons as MultiComp
from plotnine import *
import argparse

parser = argparse.ArgumentParser(description='To analyze the relationship between SS-ratio and IntronNum, at different PI levels.',
                                 epilog='Author: Fengli Zhao\tzhaofl@sustech.edu.cn      20220704 20220426 20211111 20210912 20210908')

parser.add_argument('--gpart', type=str, required=True, help='The PI-part of gene.')
parser.add_argument('--intronNum', type=str, help='The intron number file. [TAIR10.LongRNA.intron_num]', default='TAIR10.LongRNA.intron_num')
parser.add_argument('--ss', type=str, required=True, help='The intron SS-ratio result.')
parser.add_argument('--ss_type', type=str, required=True, help='The SS3, or SS5.')
parser.add_argument('--prefix', type=str, required=True, help='The prefix of output files.')

args = parser.parse_args()

gpart = args.gpart
intronNum = args.intronNum
ss = args.ss
ss_type = args.ss_type
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

idx = 0
if ss_type == 'SS3': # one-sample result
    idx = -3
else:
    idx = -4
Genes = []; Ratio = []; IntronNum =  []; iNumber = []; Part = []
loc = ''; nn = 0; ss_arr = [] # 20220704
INidx = open(ss, 'r')
while True:
    line = INidx.readline()
    if not line:
        break
    array = line.strip().split()
    if line[:3] != 'LOC':
        if array[0] in Number.keys() and array[0] in Parts.keys():
            nn += 1
            if nn == 1:
                loc = array[0]
                if array[idx] != 'unknown' and float(array[idx]) <= 1:
                    ss_arr.append(float(array[idx]))
            else:
                if array[0] == loc:
                    if array[idx] != 'unknown' and float(array[idx]) <= 1:
                        ss_arr.append(float(array[idx]))
                else:
                    if len(ss_arr) > 0:
                        ss_mean = round(sum(ss_arr) / len(ss_arr), 4)
                        Genes.append(loc); Ratio.append(ss_mean); Part.append(Parts[loc])
                        if int(Number[loc]) <= 9:
                            IntronNum.append('Intron_0' + Number[loc])
                            iNumber.append(int(Number[loc]))
                        else:
                            IntronNum.append('Intron_m10')
                            iNumber.append(int(Number[loc]))
                    loc = array[0]; ss_arr = []
                    if array[idx] != 'unknown' and float(array[idx]) <= 1:
                        ss_arr.append(float(array[idx]))
INidx.close()

if len(ss_arr) > 0:
    ss_mean = round(sum(ss_arr) / len(ss_arr), 4)
    Genes.append(loc); Ratio.append(ss_mean); Part.append(Parts[loc])
    if int(Number[loc]) <= 9:
        IntronNum.append('Intron_0' + Number[loc])
        iNumber.append(int(Number[loc]))
    else:
        IntronNum.append('Intron_m10')
        iNumber.append(int(Number[loc]))

parts = sorted(list(set(Part)))
df = pd.DataFrame({"Genes":Genes, "Ratio":Ratio, "Intron_Num":IntronNum, "iNumber":iNumber, "Part":Part})
df.to_csv(prefix + '.' + ss_type + '-ratio_IntronNum_PI-parts.data', sep='\t', index=False)

outstat = prefix + '.' + ss_type + '-ratio_IntronNum_PI-parts.statistics'
if os.path.exists(outstat):
    os.remove(outstat)

CORR = open(prefix + '.' + ss_type + '-ratio_IntronNum_PI-parts.corr', 'w')

breaks_arr = sorted(list(set(IntronNum)))
labels_arr = [x.replace('Intron_0', '') for x in breaks_arr[:-1]]
labels_arr.append(breaks_arr[-1].replace('Intron_', ''))

xPos = []; yPos = []; Count = []; cPart = [] # 20220425
xPos_r = []; yPos_r = []; Corr = []; rPart = [] # 20220426
for part in parts:
    df_p = df[df['Part'] == part]

    MultiComp.OneFactor_multiComparison(df_p, 'Intron_Num', 'Ratio', 'Fix_GenePI_' + part, outstat)

    spearman_r0, p_val0 = stats.spearmanr(df_p['iNumber'], df_p['Ratio'])

    CORR.write('The correlation between ' + ss_type + '-ratio and Intron_num in ' + prefix + '.' + part + ' :\n\n')
    CORR.write('\tspearman_r : ' + str(format(spearman_r0, '.2f')) + ' ( ' + str(spearman_r0) + ' )\n\tP_val : ' \
               + str(format(p_val0, '.2e')) + ' ( ' + str(p_val0) + ' )\n\n\n')

    xPos_r.append(breaks_arr[int(0.75 * len(breaks_arr)) - 1])
    yPos_r.append(0.95 * max(df['Ratio']))
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
        yPos.append(min(df_p['Ratio'])-0.05)
        df_g = df_p[df_p['Intron_Num'] == str0]
        Count.append(format(len(df_g['Ratio']), ','))
        cPart.append(part)

CORR.close()

df_t = pd.DataFrame({"xPos":xPos, "yPos":yPos, "Text":Count, "Part":cPart}) # 20220425
df_r = pd.DataFrame({"xPos":xPos_r, "yPos":yPos_r, "Text":Corr, "Part":rPart}) # 20220426

boxplot_fig = (ggplot(df, aes(x='Intron_Num', y='Ratio'))
+ geom_boxplot(df, aes(fill='Intron_Num'), outlier_alpha=0, outlier_size=0, position=position_dodge(0.85))
+ facet_grid('Part~.')
+ ggtitle('The relationship between ' + ss_type + '-ratio and Intron_num\nin ' + prefix + ' at different PI levels')
+ labs(x="Intron number", y="SS-ratio")
+ scale_x_discrete(breaks=breaks_arr, labels=labels_arr)
+ geom_text(df_t, aes(x='xPos', y='yPos', label='Text'), family='arial', size=8) # 20220425
+ geom_text(df_r, aes(x='xPos', y='yPos', label='Text'), family='arial', size=8) # 20220426
+ theme_bw()
+ theme(plot_title=element_text(family="arial", size=14), axis_title=element_text(family="arial", size=12), axis_text=element_text(family="arial", size=10), \
        strip_text=element_text(family="arial", size=12), legend_position="none"))
#+ theme(plot_title=element_text(family="arial", size=14), axis_title=element_text(family="arial", size=12), axis_text=element_text(family="arial", size=10), \
#        legend_title= element_text(family="arial", size=12), legend_text= element_text(family="arial", size=10)))

boxplot_fig.save(prefix + '.' + ss_type + '-ratio_IntronNum_PI-parts.boxplot.pdf', width=4.8, height=9, format='pdf')
