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
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


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

    input_gbk                 = args['gbk']
    gbk_file_ext            = args['x']
    target_to_plot          = args['tgt']
    plot_by_gene            = args['by_gene']
    plot_by_fun             = args['by_fun']
    gene_fun_txt            = args['gf']
    gene_or_fun_color_txt   = args['c']
    output_plot             = args['o']
    flk_len_to_plot         = args['l']
    plot_width              = args['width']
    plot_height             = args['height']
    no_title                = args['no_title']
    no_ruler                = args['no_ruler']
    forward_strand_target   = args['ft']

    seq_feature_csv = '%s.%s.csv' % (output_plot, flk_len_to_plot)

    # check input gbk file
    gbk_file_list = []
    if os.path.isdir(input_gbk) is True:
        gbk_file_re = '%s/*.%s' % (input_gbk, gbk_file_ext)
        gbk_file_list = glob.glob(gbk_file_re)
        if len(gbk_file_list) == 0:
            print('gbk file not found in %s. program exited!' % input_gbk)
            exit()
    elif os.path.isfile(input_gbk) is True:
        gbk_file_list = [input_gbk]
    else:
        print('gbk file not found. program exited!')
        exit()

    # read in target_to_plot
    target_gene_dict = dict()
    target_fun_dict = dict()
    if plot_by_gene is True:
        if os.path.isdir(input_gbk) is True:
            if os.path.isfile(target_to_plot) is True:
                for each_g in open(target_to_plot):
                    each_g_split = each_g.strip().split('\t')
                    gnm_id = each_g_split[0]
                    gnm_target_list = each_g_split[1].split(',')
                    pwd_gbk = '%s/%s.%s' % (input_gbk, gnm_id, gbk_file_ext)
                    target_gene_dict[pwd_gbk] = gnm_target_list
            else:
                print('Please provide to-plot target genes for all genomes in a single txt file, please separates genome (without extension) and target(s) with tab, separates targets from the same genome with comma.')
                exit()
        elif os.path.isfile(input_gbk) is True:
            gbk_path, gbk_base, gbk_ext = sep_path_basename_ext(input_gbk)
            if os.path.isfile(target_to_plot) is True:
                for each_g in open(target_to_plot):
                    each_g_split = each_g.strip().split('\t')
                    gnm_id = each_g_split[0]
                    gnm_target_list = each_g_split[1].split(',')
                    if gnm_id == gbk_base:
                        target_gene_dict[input_gbk] = gnm_target_list
            else:
                target_list = target_to_plot.split(',')
                target_gene_dict[input_gbk] = target_list
    elif plot_by_fun is True:
        target_list = target_to_plot.split(',')
        for each_gbk in gbk_file_list:
            target_fun_dict[each_gbk] = target_list
    else:
        print('Please specify either "-by_gene" or "-by_fun"')


    # read in gene function
    gene_fun_dict = dict()
    if gene_fun_txt is not None:
        if os.path.isfile(gene_fun_txt) is True:
            for each_gene in open(gene_fun_txt):
                each_gene_split = each_gene.strip().split('\t')
                gene_fun_dict[each_gene_split[0]] = each_gene_split[1]
        else:
            print('%s not found, program exited!' % gene_fun_txt)
            exit()

    # read in gene/function color
    fun_color_dict = dict()
    if gene_or_fun_color_txt is not None:
        if os.path.isfile(gene_or_fun_color_txt) is True:
            for each_fun in open(gene_or_fun_color_txt):
                each_fun_split = each_fun.strip().split('\t')
                fun_color_dict[each_fun_split[0]] = each_fun_split[1]
        else:
            print('%s not found, program exited!' % gene_or_fun_color_txt)
            exit()

    gnm_to_gene_dict = target_gene_dict
    if plot_by_fun is True:
        gnm_to_gene_dict = dict()
        for each_gbk in gbk_file_list:
            gnm_fun_to_plot = target_fun_dict.get(each_gbk)
            for seq_record in SeqIO.parse(each_gbk, 'genbank'):
                for seq_feature in seq_record.features:
                    if seq_feature.type == "CDS":
                        locus_tag = seq_feature.qualifiers['locus_tag'][0]
                        gene_fun = gene_fun_dict.get(locus_tag, 'na')
                        if plot_by_fun is True:
                            if gene_fun in gnm_fun_to_plot:
                                if each_gbk not in gnm_to_gene_dict:
                                    gnm_to_gene_dict[each_gbk] = {locus_tag}
                                else:
                                    gnm_to_gene_dict[each_gbk].add(locus_tag)

    print('target_gene_dict\t%s'    % target_gene_dict)
    print('target_fun_dict\t\t%s'   % target_fun_dict)
    print('gene_fun_dict\t\t%s'     % gene_fun_dict)
    print('fun_color_dict\t\t%s'    % fun_color_dict)
    print('gnm_to_gene_dict\t%s'  % gnm_to_gene_dict)

    # write out features
    print('Extracting features to plot')
    with open(seq_feature_csv, 'w') as seq_feature_csv_handle:
        for each_gnm in gnm_to_gene_dict:
            pwd_gbk = each_gnm
            locus_tag_to_plot = gnm_to_gene_dict[each_gnm]
            features_to_write_list = gbk_to_feature_list(pwd_gbk, locus_tag_to_plot, flk_len_to_plot, forward_strand_target, fun_color_dict)
            seq_feature_csv_handle.write('\n'.join(features_to_write_list) + '\n')

    # plot
    print('Getting plot')
    plot_ctg_feature_txt(seq_feature_csv, no_title, no_ruler, plot_height, plot_width, output_plot)

    print('Plot exported to %s' % output_plot)
    print('Done!')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-gbk',         required=True,                          help='gbk file')
    parser.add_argument('-x',           required=False,                         default='gbk', help='file extension, default: gbk')
    parser.add_argument('-tgt',         required=True,                          help='targets to plot')
    parser.add_argument('-by_gene',     required=False, action='store_true',    help='target plot region by gene id (locus_tag)')
    parser.add_argument('-by_fun',      required=False, action='store_true',    help='target plot region by gene function')
    parser.add_argument('-gf',          required=False, default=None,           help='gene to function file, tab separated')
    parser.add_argument('-c',           required=False, default=None,           help='gene/function to color file, tab separated')
    parser.add_argument('-o',           required=True,                          help='output plot')
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
python3 /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/DnaFeaturesViewer.py -gbk gbk_dir -x gbk -tgt  -by_gene -o demo.pdf

gene_to_fun.txt
fun_to_color.txt

# single gbk
python3 /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/DnaFeaturesViewer.py -gbk gbk_dir/Bradyrhizobium_oligotrophicum_SZCCHNR1093.gbk -by_gene -o demo_by_g.pdf -tgt SZCCHNR1093_01706,SZCCHNR1093_01729
python3 /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/DnaFeaturesViewer.py -gbk gbk_dir/Bradyrhizobium_oligotrophicum_SZCCHNR1093.gbk -by_gene -o demo_by_g.pdf -tgt targeted_gene.txt

# multiple gbk
python3 /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/DnaFeaturesViewer.py -gbk gbk_dir -x gbk -by_gene -o demo_by_g.pdf -tgt targeted_gene.txt





python3 /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/DnaFeaturesViewer.py -gbk gbk_dir -x gbk -by_fun -o demo_by_f.pdf -tgt shc,mntB -c fun_color.txt
python3 /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/DnaFeaturesViewer.py -gbk gbk_dir/Bradyrhizobium_oligotrophicum_SZCCHNR1091.gbk -by_fun -o demo_by_f.pdf -tgt shc,mntB -c fun_color.txt

single file/multi file

gene by str (for single file)
gene by file (for multi file)

fun by str
fun by file


'''
