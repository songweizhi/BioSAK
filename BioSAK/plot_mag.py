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


Plot_MAG_parser_usage = '''
====================================== Plot_MAG example commands ======================================

# annotate protein sequences
BioSAK Plot_MAG -i refined_MAGs -x fasta -d contig_depth.txt -o refined_MAGs.pdf

# depth file: tab separated, no header

======================================================================================================
'''


def get_color_list(color_num):

    if color_num <= 8:
        color_list_combined = ['#3787c0', '#39399f', '#ffb939', '#399f39', '#9f399f', '#fb694a', '#9f9f39', '#959595']

    elif 8 < color_num <= 16:
        color_list_combined = ['#2b7bba', '#89bedc', '#2e2e99', '#8a8acc', '#ffa500', '#ffc55c', '#2e992e', '#8acc8a', '#992e99', '#cc8acc', '#d52221', '#fc8161', '#99992e', '#cccc8a', '#5c5c5c', '#adadad']

    else:
        color_num_each = math.ceil(color_num/8) + 2

        color_list_1 = sns.color_palette('Blues',  n_colors=color_num_each).as_hex()
        color_list_2 = sns.light_palette('navy',   n_colors=color_num_each).as_hex()
        color_list_3 = sns.light_palette('orange', n_colors=color_num_each).as_hex()
        color_list_4 = sns.light_palette('green',  n_colors=color_num_each).as_hex()
        color_list_5 = sns.light_palette('purple', n_colors=color_num_each).as_hex()
        color_list_6 = sns.color_palette('Reds',   n_colors=color_num_each).as_hex()
        color_list_7 = sns.light_palette('olive',  n_colors=color_num_each).as_hex()
        color_list_8 = sns.color_palette('Greys', n_colors=color_num_each).as_hex()

        color_list_combined = []
        for color_list in [color_list_1, color_list_2, color_list_3, color_list_4, color_list_5, color_list_6, color_list_7, color_list_8]:
            for color in color_list[2:][::-1]:
                color_list_combined.append(color)

    color_list_to_return = random.sample(color_list_combined, color_num)

    color_list_to_return_sorted = []
    for color_to_return in color_list_combined:
        if color_to_return in color_list_to_return:
            color_list_to_return_sorted.append(color_to_return)

    return color_list_to_return_sorted


def Plot_MAG(argument_list):

    bin_folder =  argument_list['i']
    bin_ext =     argument_list['x']
    depth_file =  argument_list['d']
    output_plot = argument_list['o']

    # read in depth info
    ctg_depth_dict = dict()
    for each_ctg in open(depth_file):
        each_ctg_split = each_ctg.strip().split('\t')
        if not each_ctg.startswith('contigName	contigLen	totalAvgDepth	sample1	sample1-var'):
            ctg_depth_dict[each_ctg_split[0]] = float(each_ctg_split[1])

    # get mag file list
    bin_file_re = '%s/*.%s' % (bin_folder, bin_ext)
    bin_file_list = [os.path.basename(file_name) for file_name in glob.glob(bin_file_re)]

    if len(bin_file_list) == 0:
        print('No MAG detected, program exited')
        exit()

    color_list = get_color_list(len(bin_file_list))

    legend_elements = []
    num_list_depth = []
    num_list_gc = []
    ctg_color_list = []
    ctg_len_list = []
    index = 0
    for each_mag in bin_file_list:
        pwd_each_mag = '%s/%s' % (bin_folder, each_mag)
        for each_ctg in SeqIO.parse(pwd_each_mag, 'fasta'):
            ctg_depth = ctg_depth_dict[each_ctg.id]
            ctg_seq = str(each_ctg.seq).upper()
            ctg_gc = ((ctg_seq.count('G')) + (ctg_seq.count('C')))*100/len(ctg_seq)
            num_list_depth.append(ctg_depth)
            num_list_gc.append(ctg_gc)
            ctg_color_list.append(color_list[index])
            ctg_len_list.append(len(ctg_seq)/2000)
        legend_elements.append(Patch(facecolor=color_list[index], edgecolor=color_list[index], label=each_mag))
        index += 1

    num_array_depth = np.array(num_list_depth)
    num_array_gc = np.array(num_list_gc)
    num_array_len = np.array(ctg_len_list)

    plt.scatter(num_array_depth, num_array_gc, c=ctg_color_list, s=num_array_len, alpha=0.3, linewidths=0)
    #plt.scatter(num_array_depth, num_array_gc, c=ctg_color_list, s=5, alpha=0.6, linewidths=0)
    plt.title('Number of genome: %s' % len(bin_file_list), fontsize=9)
    plt.xlabel("Depth (X)", size=9)
    plt.ylabel("GC (%)", size=9)
    plt.legend(handles=legend_elements, loc='center left', prop={'size': 6}, bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(output_plot)
    plt.close()


if __name__ == '__main__':

    Plot_MAG_parser = argparse.ArgumentParser()
    Plot_MAG_parser.add_argument('-i',   required=True,  help='MAG folder')
    Plot_MAG_parser.add_argument('-x',   required=True,  help='file extension')
    Plot_MAG_parser.add_argument('-d',   required=True,  help='contig depth')
    Plot_MAG_parser.add_argument('-o',   required=True,  help='output plot')
    args = vars(Plot_MAG_parser.parse_args())
    Plot_MAG(args)
