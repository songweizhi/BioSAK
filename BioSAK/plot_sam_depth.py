
import os
import argparse
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


def take_kmer_mean(num_list, k_mer):

    k_mer_average_list = []
    n = 0

    while n < len(num_list):

        if n <= len(num_list) - k_mer:

            to_add_value = list(range(1, k_mer + 1))

            pos_list = []
            for each in to_add_value:
                pos_list.append(num_list[n + each - 1])

            # get average
            current_average = sum(pos_list) / k_mer
            k_mer_average_list.append(float("{0:.2f}".format(current_average)))

            n += 1

        elif (len(num_list) - k_mer) < n < len(num_list):
            k_mer_average_list.append(num_list[n])

            n += 1

    return k_mer_average_list


def plot_seq_depth(depth_file, seq_to_plot, start_pos, end_pos, bps_to_marker, plot_filename):
    #print('Extracting absolute depth from input file')
    x = []
    y = []
    bp_num = 1
    current_pos = 0
    for each_base in open(depth_file):
        each_base_split = each_base.strip().split('\t')
        seq_id = each_base_split[0]
        pos = int(each_base_split[1])
        depth = int(each_base_split[2])

        if seq_id == seq_to_plot:

            if pos < start_pos:
                pass

            elif pos == start_pos:
                x.append(pos)
                y.append(depth)
                current_pos = pos

            else:
                start_0 = None
                end_0 = None
                to_add = []

                if (pos == current_pos + 1) and (pos <= end_pos):
                    x.append(pos)
                    y.append(depth)
                    current_pos = pos

                elif (pos > current_pos + 1) and (pos <= end_pos):

                    # add zero
                    start_0 = current_pos + 1
                    end_0 = pos - 1
                    to_add = list(range(start_0, end_0 + 1))
                    for each_0 in to_add:
                        x.append(each_0)
                        y.append(0)

                    x.append(pos)
                    y.append(depth)
                    current_pos = pos

                elif (pos > current_pos + 1) and (pos > end_pos):

                    # add zero
                    start_0 = current_pos + 1
                    end_0 = end_pos
                    to_add = list(range(start_0, end_0 + 1))
                    for each_0 in to_add:
                        x.append(each_0)
                        y.append(0)

                    current_pos = pos

    #print('Calculating k-mer means')
    y = take_kmer_mean(y, k_mer)

    fig = plt.figure(figsize=(plot_width, plot_height), dpi=100)

    # Change the color and its transparency
    plt.plot(x, y, color="skyblue", alpha=0.7, linewidth=0.7)

    # titles
    plt.title(plot_filename,fontsize=7)
    plt.xlabel('Position (bp)', fontsize=10)
    plt.ylabel('Depth (X)', fontsize=10)

    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)

    ymax = 0
    if max(y) <= 2000:
        ymax = 2000
    elif (2000 < max(y) <= 5000):
        ymax = 5000
    elif (5000 < max(y) <= 10000):
        ymax = 10000
    else:
        ymax = max(y)

    plt.ylim(ymin=0, ymax=ymax)

    # add lines to specified positions
    if bps_to_marker != None:
        bps_to_marker_list = bps_to_marker.strip().split(',')
        for each_line in bps_to_marker_list:
            plt.axvline(x=int(each_line), c='red', linewidth=0.3)

    # Get plot
    plt.savefig('%s.png' % plot_filename, dpi=300)
    plt.close()


def plot_sam_depth(args):

    sequence_file = args['REF']
    depth_file = args['DEPTH']
    seq_to_plot = args['SeqID']
    start_pos = args['START']
    end_pos = args['END']
    k_mer = args['Kmer']
    plot_filename = args['Out']
    bps_to_marker = args['Lines']
    plot_width = args['xlen']
    plot_height = args['ylen']


    # get sequence length dict
    seq_id_length_dict = {}
    for each_seq in SeqIO.parse(sequence_file, 'fasta'):
        seq_id_length_dict[each_seq.id] = len(each_seq.seq)


    if (seq_to_plot != None) and (seq_to_plot not in seq_id_length_dict):
        print('Reference sequence %s not found in bam file, program exited!' % seq_to_plot)
        exit()


    # get depth_file base name
    depth_file_basename, depth_file_extension = os.path.splitext(depth_file)


    if seq_to_plot != None:

        # get the start and end of the region to plot
        if start_pos == None:
            start_pos = 1
        if end_pos == None:
            end_pos = seq_id_length_dict[seq_to_plot]

        print('Processing %s' % seq_to_plot)
        if plot_filename == None:
            plot_filename = '%s__%s__%s-%sbp__%smer' % (depth_file_basename, seq_to_plot, start_pos, end_pos, k_mer)
            plot_seq_depth(depth_file, seq_to_plot, start_pos, end_pos, bps_to_marker, plot_filename)

    if seq_to_plot == None:

        for each_ctg in SeqIO.parse(sequence_file, 'fasta'):
            print('Processing %s' % each_ctg.id)
            plot_filename = '%s__%s__%s-%sbp__%smer' % (depth_file_basename, each_ctg.id, 1, seq_id_length_dict[each_ctg.id], k_mer)
            plot_seq_depth(depth_file, each_ctg.id, 1, seq_id_length_dict[each_ctg.id], bps_to_marker, plot_filename)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='', add_help=False)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument('-h', action='help', help='Show this help message and exit')
    required.add_argument('-r', dest='REF',   nargs='?', required=True,  type=str, help='reference sequence file')
    required.add_argument('-d', dest='DEPTH', nargs='?', required=True,  type=str, help='depth file')
    optional.add_argument('-i', dest='SeqID', nargs='?', required=False, type=str, default=None, help='id of sequence to plot')
    optional.add_argument('-s', dest='START', nargs='?', required=False, type=int, default=None, help='start position to plot')
    optional.add_argument('-e', dest='END',   nargs='?', required=False, type=int, default=None, help='end position to plot')
    optional.add_argument('-k', dest='Kmer',  nargs='?', required=False, type=int, default=100, help='k-mer mean depth')
    optional.add_argument('-o', dest='Out',   nargs='?', required=False, type=str, default=None, help='output plot name')
    optional.add_argument('-l', dest='Lines', nargs='?', required=False, type=str, default=None, help='output plot name')
    optional.add_argument('-x', dest='xlen',  nargs='?', required=False, type=int, default=8, help='plot width')
    optional.add_argument('-y', dest='ylen',  nargs='?', required=False, type=int, default=3, help='plot height')

    args = vars(parser.parse_args())

    plot_sam_depth(args)
