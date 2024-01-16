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
================= Plot_MAG example commands =================

BioSAK Plot_MAG -i MAG_dir -x fa -d ctg_depth.txt 
BioSAK Plot_MAG -i MAG_dir -x fa -d ctg_depth.txt -sep

# depth file: tab separated, no header
contig_1    30
contig_2    5

=============================================================
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


def just_plot(num_list_depth, num_list_gc, ctg_len_list, ctg_color_list, legend_elements, output_plot):

    ctg_len_list_adjusted = [l/2 for l in ctg_len_list]

    num_array_depth = np.array(num_list_depth)
    num_array_gc = np.array(num_list_gc)
    num_array_len = np.array(ctg_len_list_adjusted)

    plt.scatter(num_array_depth, num_array_gc, c=ctg_color_list, s=num_array_len, alpha=0.3, linewidths=0)
    #plt.scatter(num_array_depth, num_array_gc, c=ctg_color_list, s=5, alpha=0.6, linewidths=0)
    #plt.title('Number of genome: %s' % len(bin_file_list), fontsize=9)
    plt.title('Range of sequence length (Kbp): %s - %s' % (min(ctg_len_list), max(ctg_len_list)), fontsize=8)
    plt.xlabel("Depth (X)", size=9)
    plt.ylabel("GC (%)", size=9)
    plt.legend(handles=legend_elements, title='Genome', loc='upper left', prop={'size': 6}, bbox_to_anchor=(1, 1))
    #plt.legend(handles=legend_elements, loc='lower left', prop={'size': 6}, bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(output_plot)
    plt.close()


def just_plot_backup(num_list_depth, num_list_gc, ctg_len_list, ctg_color_list, legend_elements, output_plot):

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


def Plot_MAG(argument_list):

    bin_folder = argument_list['i']
    bin_ext    = argument_list['x']
    depth_file = argument_list['d']
    plot_sep   = argument_list['sep']

    # read in depth info
    ctg_depth_dict = dict()
    for each_ctg in open(depth_file):
        each_ctg_split = each_ctg.strip().split('\t')
        ctg_depth_dict[each_ctg_split[0]] = float(each_ctg_split[1])

    # get mag file list
    bin_file_re = '%s/*.%s' % (bin_folder, bin_ext)
    bin_file_list = [os.path.basename(file_name) for file_name in glob.glob(bin_file_re)]

    if len(bin_file_list) == 0:
        print('No MAG detected, program exited')
        exit()

    color_list = get_color_list(len(bin_file_list))

    output_txt  = '%s_plot_GC_vs_depth.txt'     % bin_folder
    output_plot = '%s_plot_GC_vs_depth.pdf'     % bin_folder
    plot_folder = '%s_plot_GC_vs_depth'         % bin_folder
    txt_folder  = '%s_plot_GC_vs_depth/data'    % bin_folder
    if plot_sep is True:
        os.mkdir(plot_folder)
        os.mkdir(txt_folder)

    output_txt_handle = open(output_txt, 'w')
    output_txt_handle.write('Genome\tContig\tDepth\tGC\tLength(bp)\n')
    legend_elements = []
    num_list_depth = []
    num_list_gc = []
    ctg_len_list = []
    ctg_color_list = []
    index = 0
    for each_mag in bin_file_list:

        mag_no_ext        = each_mag[:-(len(bin_ext) + 1)]
        pwd_each_mag      = '%s/%s'     % (bin_folder, each_mag)
        pwd_each_mag_plot = '%s/%s.pdf' % (plot_folder, mag_no_ext)
        pwd_each_mag_txt  = '%s/%s.txt' % (txt_folder, mag_no_ext)

        if plot_sep is True:
            pwd_each_mag_txt_handle = open(pwd_each_mag_txt, 'w')
            pwd_each_mag_txt_handle.write('Contig\tDepth\tGC\tLength(bp)\n')
        legend_elements_sep = []
        num_list_depth_sep = []
        num_list_gc_sep = []
        ctg_len_list_sep = []
        ctg_color_list_sep = []
        for each_ctg in SeqIO.parse(pwd_each_mag, 'fasta'):
            ctg_depth = ctg_depth_dict[each_ctg.id]
            ctg_seq = str(each_ctg.seq).upper()
            ctg_gc = ((ctg_seq.count('G')) + (ctg_seq.count('C')))*100/len(ctg_seq)
            ctg_len_kbp = len(ctg_seq) / (1024)
            ctg_len_kbp = float("{0:.2f}".format(ctg_len_kbp))

            num_list_depth.append(ctg_depth)
            num_list_depth_sep.append(ctg_depth)

            num_list_gc.append(ctg_gc)
            num_list_gc_sep.append(ctg_gc)

            ctg_len_list.append(ctg_len_kbp)
            ctg_len_list_sep.append(ctg_len_kbp)

            ctg_color_list.append(color_list[index])
            ctg_color_list_sep.append(color_list[index])

            output_txt_handle.write('%s\t%s\t%s\t%s\t%s\n' % (mag_no_ext, each_ctg.id, ctg_depth, float("{0:.2f}".format(ctg_gc)), len(ctg_seq)))
            if plot_sep is True:
                pwd_each_mag_txt_handle.write('%s\t%s\t%s\t%s\n' % (each_ctg.id, ctg_depth, float("{0:.2f}".format(ctg_gc)), len(ctg_seq)))

        legend_elements.append(Patch(facecolor=color_list[index], edgecolor=color_list[index], label=mag_no_ext))
        legend_elements_sep.append(Patch(facecolor=color_list[index], edgecolor=color_list[index], label=mag_no_ext))
        if plot_sep is True:
            just_plot(num_list_depth_sep, num_list_gc_sep, ctg_len_list_sep, ctg_color_list_sep, legend_elements_sep, pwd_each_mag_plot)
            pwd_each_mag_txt_handle.close()
        index += 1

    just_plot(num_list_depth, num_list_gc, ctg_len_list, ctg_color_list, legend_elements, output_plot)
    output_txt_handle.close()


if __name__ == '__main__':

    Plot_MAG_parser = argparse.ArgumentParser(usage=Plot_MAG_parser_usage)
    Plot_MAG_parser.add_argument('-i',   required=True,                         help='MAG folder')
    Plot_MAG_parser.add_argument('-x',   required=True,                         help='file extension')
    Plot_MAG_parser.add_argument('-d',   required=True,                         help='contig depth')
    Plot_MAG_parser.add_argument('-sep', required=False, action='store_true',   help='get plot for individual MAGs')
    args = vars(Plot_MAG_parser.parse_args())
    Plot_MAG(args)
