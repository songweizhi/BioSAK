import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation
from reportlab.lib import colors
from reportlab.lib.units import cm


VisGeneFlk_usage = '''
====================== VisGeneFlk example commands ======================

BioSAK VisGeneFlk -gene gene_01 -gbk input.gbk -len 2000 -fmt pdf
BioSAK VisGeneFlk -gene gene_01 -gbk input.gbk -len 5000 -scale 300

=========================================================================
'''


def get_flanking_region(input_gbk_file, HGT_candidate, flanking_length):

    wd, gbk_file = os.path.split(input_gbk_file)
    new_gbk_file = '%s/%s_%sbp_temp.gbk' % (wd, HGT_candidate, flanking_length)
    new_gbk_final_file = '%s/%s_%sbp.gbk' % (wd, HGT_candidate, flanking_length)

    # get flanking range of candidate
    input_gbk = SeqIO.parse(input_gbk_file, "genbank")
    new_start = 0
    new_end = 0
    contig_length = 0
    for record in input_gbk:
        contig_length = len(record.seq)
        for gene in record.features:
            # get contig length
            if gene.type == 'source':
                pass
            # get new start and end points
            elif 'locus_tag' in gene.qualifiers:
                if gene.qualifiers['locus_tag'][0] == HGT_candidate:
                    # get new start point
                    new_start = gene.location.start - flanking_length
                    if new_start < 0:
                        new_start = 0
                    # get new end point
                    new_end = gene.location.end + flanking_length
                    if new_end > contig_length:
                        new_end = contig_length

    # get genes within flanking region
    keep_gene_list = []
    input_gbk = SeqIO.parse(input_gbk_file, "genbank")
    for record in input_gbk:
        for gene in record.features:
            if 'locus_tag' in gene.qualifiers:
                if (gene.location.start < new_start) and (gene.location.end >= new_start):
                    keep_gene_list.append(gene.qualifiers['locus_tag'][0])
                    new_start = gene.location.start
                elif (gene.location.start > new_start) and (gene.location.end < new_end):
                    keep_gene_list.append(gene.qualifiers['locus_tag'][0])
                elif (gene.location.start <= new_end) and (gene.location.end > new_end):
                    keep_gene_list.append(gene.qualifiers['locus_tag'][0])
                    new_end = gene.location.end

    # remove genes not in flanking region from gbk file
    input_gbk = SeqIO.parse(input_gbk_file, "genbank")
    new_gbk = open(new_gbk_file, 'w')
    for record in input_gbk:
        new_record_features = []
        for gene in record.features:
            if gene.type == 'source':
                new_record_features.append(gene)
            elif 'locus_tag' in gene.qualifiers:
                if gene.qualifiers['locus_tag'][0] in keep_gene_list:
                    new_record_features.append(gene)
        record.features = new_record_features
        SeqIO.write(record, new_gbk, 'genbank')
    new_gbk.close()

    # remove sequences not in flanking region
    new_gbk_full_length = SeqIO.parse(new_gbk_file, "genbank")
    new_gbk_final = open(new_gbk_final_file, 'w')
    for record in new_gbk_full_length:

        # get new sequence
        new_seq = record.seq[new_start:new_end]
        new_contig_length = len(new_seq)
        new_record = SeqRecord(new_seq, id=record.id, name=record.name,
                               description=record.description,
                               annotations=record.annotations)

        # get new location
        new_record_features_2 = []
        for gene in record.features:
            if gene.type == 'source':
                gene_location_new = ''
                if gene.location.strand == 1:
                    gene_location_new = FeatureLocation(0, new_contig_length, strand=+1)
                if gene.location.strand == -1:
                    gene_location_new = FeatureLocation(0, new_contig_length, strand=-1)
                gene.location = gene_location_new
                new_record_features_2.append(gene)
            elif 'locus_tag' in gene.qualifiers:

                strand_index = ''
                if gene.location.strand == 1:
                    strand_index = +1
                if gene.location.strand == -1:
                    strand_index = -1

                if gene.location.start - new_start < 0:
                    gene_location_new = FeatureLocation(0, gene.location.end - new_start, strand=strand_index)
                else:
                    gene_location_new = FeatureLocation(gene.location.start - new_start, gene.location.end - new_start, strand=strand_index)

                gene.location = gene_location_new
                new_record_features_2.append(gene)
        new_record.features = new_record_features_2
        SeqIO.write(new_record, new_gbk_final, 'genbank')

    new_gbk_final.close()
    os.system('rm %s' % new_gbk_file)


def set_contig_track_features(gene_contig, gene_to_highlight, feature_set, hide_label):

    if hide_label is False:
        show_label = True
    else:
        show_label = False

    # add features to feature set
    for feature in gene_contig.features:
        if feature.type == "CDS":

            # define label color
            if feature.qualifiers['locus_tag'][0] == gene_to_highlight:
                label_color = colors.red
                label_size = 10
            else:
                label_color = colors.black
                label_size = 10

            # strands
            color = None
            label_angle = 0
            if feature.location.strand == 1:
                label_angle = 45
                color = colors.lightblue
            elif feature.location.strand == -1:
                label_angle = -225
                color = colors.lightgreen

            # add feature
            feature_set.add_feature(feature, color=color, sigil='BIGARROW', arrowshaft_height=0.5, arrowhead_length=0.3,
                                    label=show_label, label_color=label_color, label_size=label_size, label_angle=label_angle, label_position="middle")


def VisGeneFlk(args):

    gene_id         = args['gene']
    input_gbk       = args['gbk']
    flanking_length = args['len']
    plot_scale      = args['scale']
    plot_fmt        = args['fmt']
    hide_label      = args['no_label']

    plot_wd                  = '%s_flk%s_wd'    % (gene_id, flanking_length)
    gbk_subset_located_seq   = '%s/%s.gbk'      % (plot_wd, gene_id)
    gbk_subset_flanking_gene = '%s/%s_%sbp.gbk' % (plot_wd, gene_id, flanking_length)
    plot_file                = '%s_flk%sbp.%s'  % (gene_id, flanking_length, plot_fmt)

    if os.path.isdir(plot_wd) is False:
        os.mkdir(plot_wd)
    else:
        os.system('rm -r %s' % plot_wd)
        os.mkdir(plot_wd)

    dict_value_list = []
    for seq_record in SeqIO.parse(input_gbk, 'genbank'):
        for gene_feature in seq_record.features:
            if 'locus_tag' in gene_feature.qualifiers:
                if gene_id in gene_feature.qualifiers["locus_tag"]:
                    dict_value_list.append([gene_id, int(gene_feature.location.start), int(gene_feature.location.end), gene_feature.location.strand, len(seq_record.seq)])
                    SeqIO.write(seq_record, gbk_subset_located_seq, 'genbank')
                    get_flanking_region(gbk_subset_located_seq, gene_id, flanking_length)

    # get the distance of the gene to contig ends
    gene_1_left_len = dict_value_list[0][1]
    gene_1_right_len = dict_value_list[0][4] - dict_value_list[0][2]

    # read in gbk file
    sequence_record = SeqIO.read(gbk_subset_flanking_gene, "genbank")

    # create an empty diagram
    diagram = GenomeDiagram.Diagram()
    plot_len_cm = len(sequence_record)/plot_scale

    # add tracks to diagram
    track_footnote = '%s (left %sbp, right %sbp)' % (sequence_record.name, gene_1_left_len, gene_1_right_len)
    track_footnote = sequence_record.name
    seq_track = diagram.new_track(1, name=track_footnote, greytrack=True,
                                  greytrack_labels=1, greytrack_font='Helvetica', greytrack_fontsize=12,
                                  height=0.35, start=0, end=len(sequence_record),
                                  scale=True, scale_fontsize=6, scale_ticks=1,
                                  scale_smalltick_interval=10000, scale_largetick_interval=10000)

    # create blank feature set and add gene features to it
    feature_set = seq_track.new_set(type='feature')
    set_contig_track_features(sequence_record, gene_id, feature_set, hide_label)

    # draw and export
    diagram.draw(format='linear', orientation='landscape', pagesize=(20*cm, plot_len_cm*cm), fragments=1, start=0, end=len(sequence_record))
    diagram.write(plot_file, plot_fmt)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-gene',    required=True,                          help='gene id')
    parser.add_argument('-gbk',     required=True,                          help='gbk file')
    parser.add_argument('-len',     required=True, type=int,                help='length (in bp) of flanking sequences to plot')
    parser.add_argument('-scale',   required=False, type=int, default=200,  help='scale for plotting, default: 200bp per cm')
    parser.add_argument('-fmt',     required=False, default='svg',          help='output format (svg or pdf), default: svg')
    parser.add_argument('-no_label',required=False, action='store_true',    help='output format (svg or pdf), default: svg')
    args = vars(parser.parse_args())
    VisGeneFlk(args)
