import os
import glob
import math
import random
import argparse
import numpy as np
import seaborn as sns
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def just_plot_backup(num_list_depth, num_list_gc, ctg_len_list, ctg_color_list, bin_file_list, legend_elements, output_plot):

    num_array_depth = np.array(num_list_depth)
    num_array_gc = np.array(num_list_gc)
    num_array_len = np.array(ctg_len_list)
    fig, ax = plt.subplots()
    scatter = plt.scatter(num_array_depth, num_array_gc, c=ctg_color_list, s=num_array_len, alpha=0.3, linewidths=0)
    plt.title('Contig length (Kbp): %s - %s' % (min(ctg_len_list), max(ctg_len_list)), fontsize=9)
    plt.xlabel("Depth (X)", size=9)
    plt.ylabel("GC (%)", size=9)

    # add size legend
    # kw = dict(prop="sizes", num=6, color='grey', fmt="{x:.2f}", func=lambda s: np.sqrt(s/.3)/3)
    # handles, labels = scatter.legend_elements(**kw, alpha=0.5)
    # legend2 = ax.legend(handles, labels, loc="upper left", title="Size", prop={'size': 5}, bbox_to_anchor=(1, 1))
    # ax.add_artist(legend2)

    # add color legend
    plt.legend(handles=legend_elements, title='MAG', loc='upper left', prop={'size': 6}, bbox_to_anchor=(1, 0.5))

    # save and close
    plt.tight_layout()
    plt.savefig(output_plot)
    plt.close()


num_list_depth  = []
num_list_gc     = []
ctg_len_list    = []
ctg_color_list  = ['green']
bin_file_list   = ['bin1.fa']
output_plot     = '/Users/songweizhi/Desktop/888/plot.pdf'


pwd_each_mag    = '/Users/songweizhi/Desktop/888/4.fa'
depth_file      = '/Users/songweizhi/Desktop/888/filtered_stats_Maxbin_nonEukaryota.txt'


# read in depth info
ctg_depth_dict = dict()
for each_ctg in open(depth_file):
    each_ctg_split = each_ctg.strip().split('\t')
    ctg_depth_dict[each_ctg_split[0]] = float(each_ctg_split[1])

for each_ctg in SeqIO.parse(pwd_each_mag, 'fasta'):
    ctg_depth = ctg_depth_dict[each_ctg.id]
    ctg_seq = str(each_ctg.seq).upper()
    ctg_gc = ((ctg_seq.count('G')) + (ctg_seq.count('C'))) * 100 / len(ctg_seq)
    num_list_depth.append(ctg_depth)
    num_list_gc.append(ctg_gc)
    ctg_len_kbp = len(ctg_seq)/(1024)
    ctg_len_kbp = float("{0:.2f}".format(ctg_len_kbp))
    ctg_len_list.append(ctg_len_kbp)

legend_elements = []
legend_elements.append(Patch(facecolor='green', edgecolor='green', label='bin1.fa'))

just_plot(num_list_depth, num_list_gc, ctg_len_list, ctg_color_list, bin_file_list, legend_elements, output_plot)


