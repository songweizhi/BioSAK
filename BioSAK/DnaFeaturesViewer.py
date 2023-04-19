import os
import glob
import argparse
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


DnaFeaturesViewer_usage = '''
============================ DnaFeaturesViewer example commands ============================

BioSAK DnaFeaturesViewer -target shc -gbk gbk_dir -x gbk -o demo.pdf -c gene_color.txt -ft

============================================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def gbk_to_feature_list(input_gbk, target_locus_tag_list, flk_len_to_plot, forward_strand_target, gene_color_dict):

    # get contigs to plot region dict
    ctg_to_plot_region_dict = dict()
    ctg_to_plot_region_len_dict = dict()
    locus_tag_to_direction_dict = dict()
    for seq_record in SeqIO.parse(input_gbk, 'genbank'):
        contig_id  = seq_record.id
        contig_len = len(seq_record.seq)
        for seq_feature in seq_record.features:
            if seq_feature.type == "CDS":
                locus_tag = seq_feature.qualifiers['locus_tag'][0]
                if locus_tag in target_locus_tag_list:

                    pos_start = seq_feature.location.start
                    pos_end = seq_feature.location.end
                    feature_strand = seq_feature.strand
                    locus_tag_to_direction_dict[locus_tag] = feature_strand

                    to_plot_start = pos_start - flk_len_to_plot
                    if to_plot_start < 1:
                        to_plot_start = 1

                    to_plot_end = pos_end + flk_len_to_plot
                    if to_plot_end > contig_len:
                        to_plot_end = contig_len

                    if contig_id not in ctg_to_plot_region_dict:
                        ctg_to_plot_region_dict[contig_id] = dict()
                        ctg_to_plot_region_len_dict[contig_id] = dict()

                    ctg_to_plot_region_dict[contig_id][locus_tag] = [to_plot_start, to_plot_end]
                    ctg_to_plot_region_len_dict[contig_id][locus_tag] = to_plot_end - to_plot_start + 1

    to_write_feature_list = []
    for seq_record in SeqIO.parse(input_gbk, 'genbank'):
        gbk_path, gbk_basename, gbk_ext = sep_path_basename_ext(input_gbk)
        contig_id = seq_record.id
        if contig_id in ctg_to_plot_region_dict:
            current_ctg_to_plot_region_dict = ctg_to_plot_region_dict[contig_id]
            current_ctg_to_plot_region_len_dict = ctg_to_plot_region_len_dict[contig_id]
            for target_region in current_ctg_to_plot_region_dict:
                to_plot_region = current_ctg_to_plot_region_dict[target_region]
                to_plot_region_len = current_ctg_to_plot_region_len_dict[target_region]
                target_region_direction = locus_tag_to_direction_dict[target_region]

                reverse_feature = False
                if forward_strand_target is True:
                    if target_region_direction == -1:
                        reverse_feature = True
                else:
                    if target_region_direction == 1:
                        reverse_feature = True

                location_label = '%s___%s___%s' % (gbk_basename, contig_id, target_region)
                for seq_feature in seq_record.features:
                    if seq_feature.type == "CDS":
                        locus_tag      = seq_feature.qualifiers['locus_tag'][0]
                        gene_name      = seq_feature.qualifiers.get('gene', ['na'])[0]
                        pos_start      = seq_feature.location.start
                        pos_end        = seq_feature.location.end
                        feature_strand = seq_feature.strand
                        gene_color     = gene_color_dict.get(gene_name, '')

                        to_label = gene_name
                        if to_label == 'na':
                            to_label = locus_tag

                        strand_index = ''
                        strand_index_reverse = ''
                        if feature_strand == 1:
                            strand_index = 'f'
                            strand_index_reverse = 'r'
                        if feature_strand == -1:
                            strand_index = 'r'
                            strand_index_reverse = 'f'

                        if (pos_start >= to_plot_region[0]) and (pos_end <= to_plot_region[1]):
                            pos_start_new = pos_start - to_plot_region[0]
                            pos_end_new   = pos_end - to_plot_region[0]

                            if reverse_feature is False:
                                to_write_feature_list.append('%s,%s,%s,%s,%s,%s' % (location_label, to_label, strand_index, pos_start_new, pos_end_new, gene_color))
                            else:
                                pos_start_new_reversed = to_plot_region_len - pos_end_new + 1
                                pos_end_new_reversed   = to_plot_region_len - pos_start_new + 1
                                to_write_feature_list.append('%s,%s,%s,%s,%s,%s' % (location_label, to_label, strand_index_reverse, pos_start_new_reversed, pos_end_new_reversed, gene_color))

    return to_write_feature_list


def read_in_ctg_feature(ctg_feature_txt):

    ctg_pos_dict = dict()
    ctg_feature_dict = dict()
    for each_gene in open(ctg_feature_txt):
        gene_split = each_gene.strip().split(',')
        ctg_id = gene_split[0]
        if ctg_id not in ctg_feature_dict:
            ctg_feature_dict[ctg_id] = []
            ctg_pos_dict[ctg_id] = set()

        gene_id = gene_split[1]
        if gene_id == '':
            gene_id = ' '
        gene_start = int(gene_split[3])
        gene_end = int(gene_split[4])
        gene_color = gene_split[5]
        if gene_color == '':
            gene_color = "#FFFFCC"

        ctg_pos_dict[ctg_id].add(gene_start)
        ctg_pos_dict[ctg_id].add(gene_end)

        strand_index = ''
        if gene_split[2] == 'f':
            strand_index = +1
        elif gene_split[2] == 'r':
            strand_index = -1
        else:
            print('Strand value not valid for %s!, specify as "f" or "r".' % gene_id)
            exit()

        ctg_feature_dict[ctg_id].append(GraphicFeature(start=gene_start, end=gene_end, strand=strand_index, color=gene_color, label=gene_id))

    ctg_range_dict = {i: [min(ctg_pos_dict[i]), max(ctg_pos_dict[i])] for i in ctg_pos_dict}

    return ctg_feature_dict, ctg_range_dict


def plot_ctg_feature_txt(ctg_feature_txt, hide_title, hide_ruler, plot_height, plot_width, out_plot):

    ctg_feature_dict, ctg_range_dict = read_in_ctg_feature(ctg_feature_txt)
    ctg_list_sorted = sorted([i for i in ctg_feature_dict])

    if len(ctg_list_sorted) == 1:
        ctg_id = ctg_list_sorted[0]
        ctg_feature_list = ctg_feature_dict[ctg_id]
        ctg_range_list = ctg_range_dict[ctg_id]
        graph_record = GraphicRecord(sequence_length=(ctg_range_list[1]), features=ctg_feature_list)
        ax, _ = graph_record.plot(figure_width=5)
        ax.figure.savefig(out_plot)
    else:
        fig, ax = plt.subplots(nrows=len(ctg_feature_dict), ncols=1, figsize=(plot_width, len(ctg_feature_dict)*plot_height*1.2))
        i = 0
        for each_ctg in ctg_list_sorted:
            ctg_feature_list = ctg_feature_dict[each_ctg]
            ctg_range_list = ctg_range_dict[each_ctg]
            graph_record = GraphicRecord(sequence_length=(ctg_range_list[1]), features=ctg_feature_list)

            if hide_title is False:
                ax[i].set_xlabel(each_ctg)

            graph_record.plot(ax=ax[i], with_ruler=(not hide_ruler))
            i += 1

        #plt.subplots_adjust(hspace=plot_height*0.1)
        plt.savefig(out_plot)


def DnaFeaturesViewer(args):

    targeted_gene_name      = args['target']
    gbk_dir                 = args['gbk']
    gbk_file_ext            = args['x']
    gene_color_txt          = args['c']
    flk_len_to_plot         = args['l']
    plot_width              = args['width']
    plot_height             = args['height']
    no_title                = args['no_title']
    no_ruler                = args['no_ruler']
    forward_strand_target   = args['ft']
    output_plot             = args['o']

    gbk_file_re     = '%s/*.%s'   % (gbk_dir, gbk_file_ext)
    seq_feature_csv = '%s.%s.csv' % (output_plot, flk_len_to_plot)

    # read in gene color
    gene_color_dict = dict()
    for each_gene in open(gene_color_txt):
        each_gene_split = each_gene.strip().split('\t')
        gene_color_dict[each_gene_split[0]] = each_gene_split[1]

    gnm_to_gene_dict = dict()
    for each_gbk in glob.glob(gbk_file_re):
        gbk_path, gbk_basename, gbk_ext = sep_path_basename_ext(each_gbk)
        print('Processing %s' % gbk_basename)
        for seq_record in SeqIO.parse(each_gbk, 'genbank'):
            for seq_feature in seq_record.features:
                if seq_feature.type == "CDS":
                    locus_tag = seq_feature.qualifiers['locus_tag'][0]
                    gene_name = seq_feature.qualifiers.get('gene', ['na'])[0]
                    if gene_name == targeted_gene_name:
                        if gbk_basename not in gnm_to_gene_dict:
                            gnm_to_gene_dict[gbk_basename] = {locus_tag}
                        else:
                            gnm_to_gene_dict[gbk_basename].add(locus_tag)

    # write out features
    with open(seq_feature_csv, 'w') as seq_feature_csv_handle:
        for each_gnm in gnm_to_gene_dict:
            print('Extracting features to plot: %s' % each_gnm)
            pwd_gbk = '%s/%s.%s' % (gbk_dir, each_gnm, gbk_file_ext)
            locus_tag_to_plot = gnm_to_gene_dict[each_gnm]
            features_to_write_list = gbk_to_feature_list(pwd_gbk, locus_tag_to_plot, flk_len_to_plot, forward_strand_target, gene_color_dict)
            seq_feature_csv_handle.write('\n'.join(features_to_write_list) + '\n')

    # plot
    print('Getting plot')
    plot_ctg_feature_txt(seq_feature_csv, no_title, no_ruler, plot_height, plot_width, output_plot)

    print('Plot exported to %s' % output_plot)
    print('Done!')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-target',      required=True,                          help='gene id')
    parser.add_argument('-gbk',         required=True,                          help='gbk file')
    parser.add_argument('-x',           required=False, default='gbk',          help='file extension, default: gbk')
    parser.add_argument('-o',           required=True,                          help='output plot')
    parser.add_argument('-c',           required=False, default=None,           help='gene to color file, tab separated')
    parser.add_argument('-l',           required=False, type=int, default=7500, help='length (in bp) of flanking sequences to plot')
    parser.add_argument('-width',       required=False, type=int, default=22,   help='plot_width, default: 22')
    parser.add_argument('-height',      required=False, type=int, default=2,    help='subplot_height, default: 2')
    parser.add_argument('-no_title',    required=False, action='store_true',    help='no_title')
    parser.add_argument('-no_ruler',    required=False, action='store_true',    help='no_ruler')
    parser.add_argument('-ft',          required=False, action='store_true',    help='forward target strand')
    args = vars(parser.parse_args())
    DnaFeaturesViewer(args)

'''
cd /Users/songweizhi/Desktop/Brady/draw_shc_wd
python3 /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/DnaFeaturesViewer.py -target shc -gbk gbk_dir -x gbk -o demo.pdf -c gene_color.txt -ft
'''
