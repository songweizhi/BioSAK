#!/usr/bin/env python3

import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation
from reportlab.lib import colors
from reportlab.lib.units import cm

VisGeneFlk_usage = '''
=========================== VisGeneFlk example commands ===========================

BioSAK VisGeneFlk -gene NorthSea_bin068_00045 -gbk NorthSea_bin068.gbk -len 5000

===================================================================================
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
        new_record = SeqRecord(new_seq,
                               id=record.id,
                               name=record.name,
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
                gene_location_new = ''
                if gene.location.strand == 1:
                    if gene.location.start - new_start < 0:
                        gene_location_new = FeatureLocation(0,
                                                            gene.location.end - new_start,
                                                            strand=+1)
                    else:
                        gene_location_new = FeatureLocation(gene.location.start - new_start,
                                                            gene.location.end - new_start,
                                                            strand=+1)
                if gene.location.strand == -1:
                    if gene.location.start - new_start < 0:
                        gene_location_new = FeatureLocation(0,
                                                            gene.location.end - new_start,
                                                            strand=-1)
                    else:
                        gene_location_new = FeatureLocation(gene.location.start - new_start,
                                                            gene.location.end - new_start,
                                                            strand=-1)
                gene.location = gene_location_new
                new_record_features_2.append(gene)
        new_record.features = new_record_features_2
        SeqIO.write(new_record, new_gbk_final, 'genbank')

    new_gbk_final.close()
    os.system('rm %s' % new_gbk_file)


def set_contig_track_features(gene_contig, candidate_list, feature_set):

    # add features to feature set
    for feature in gene_contig.features:
        if feature.type == "CDS":

            # define label color
            if feature.qualifiers['locus_tag'][0] in candidate_list:
                label_color = colors.blue
                label_size = 16
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
            feature_set.add_feature(feature, color=color, label=True, sigil='ARROW', arrowshaft_height=0.5, arrowhead_length=0.4, label_color=label_color, label_size=label_size, label_angle=label_angle, label_position="middle")


def VisGeneFlk(args):

    gene_1              = args['gene']
    pwd_genome_1_gbk    = args['gbk']
    flanking_length     = args['len']

    plot_wd = '%s_flk%s_wd' % (gene_1, flanking_length)

    if os.path.isdir(plot_wd) is False:
        os.mkdir(plot_wd)
    else:
        os.system('rm -r %s' % plot_wd)
        os.mkdir(plot_wd)

    dict_value_list = []
    # Extract gbk and fasta files for gene 1
    for genome_1_record in SeqIO.parse(pwd_genome_1_gbk, 'genbank'):
        for gene_1_f in genome_1_record.features:
            if 'locus_tag' in gene_1_f.qualifiers:
                if gene_1 in gene_1_f.qualifiers["locus_tag"]:
                    dict_value_list.append([gene_1, int(gene_1_f.location.start), int(gene_1_f.location.end), gene_1_f.location.strand, len(genome_1_record.seq)])
                    pwd_gene_1_gbk_file = '%s/%s.gbk' % (plot_wd, gene_1)
                    pwd_gene_1_fasta_file = '%s/%s.fasta' % (plot_wd, gene_1)
                    SeqIO.write(genome_1_record, pwd_gene_1_gbk_file, 'genbank')
                    SeqIO.write(genome_1_record, pwd_gene_1_fasta_file, 'fasta')
                    get_flanking_region(pwd_gene_1_gbk_file, gene_1, flanking_length)


    ############################## prepare for flanking plot ##############################

    # read in gbk files
    matche_pair_list = []
    path_to_gbk_file = '%s/%s_%sbp.gbk' % (plot_wd, gene_1, flanking_length)
    gene_contig = SeqIO.read(path_to_gbk_file, "genbank")
    matche_pair_list.append(gene_contig)

    # get the distance of the gene to contig ends
    gene_1_left_len = dict_value_list[0][1]
    gene_1_right_len = dict_value_list[0][4] - dict_value_list[0][2]

    # create an empty diagram
    diagram = GenomeDiagram.Diagram()

    # add tracks to diagram
    for gene1_contig in matche_pair_list:

        # add gene content track for gene1_contig
        contig_1_gene_content_track = diagram.new_track(1,
                                                        name='%s (left %sbp, right %sbp)' % (gene1_contig.name, gene_1_left_len, gene_1_right_len),
                                                        greytrack=True,
                                                        greytrack_labels=1,
                                                        greytrack_font='Helvetica',
                                                        greytrack_fontsize=12,
                                                        height=0.35,
                                                        start=0,
                                                        end=len(gene1_contig),
                                                        scale=True,
                                                        scale_fontsize=6,
                                                        scale_ticks=1,
                                                        scale_smalltick_interval=10000,
                                                        scale_largetick_interval=10000)


        # add blank feature/graph sets to each track
        feature_sets_1 = contig_1_gene_content_track.new_set(type='feature')

        # add gene features to 2 blank feature sets
        set_contig_track_features(gene1_contig, [gene_1], feature_sets_1)

        ############################################### Draw and Export ################################################

        diagram.draw(format="linear",
                     orientation="landscape",
                     pagesize=(75 * cm, 25 * cm),
                     fragments=1,
                     start=0,
                     end=len(gene1_contig))

        diagram.write('%s_flk%sbp.svg' % (gene_1, flanking_length), "svg")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-gene',    required=True,           help='gene id')
    parser.add_argument('-gbk',     required=True,           help='gbk file')
    parser.add_argument('-len',     required=True, type=int, help='length of flanking region to plot')
    args = vars(parser.parse_args())
    VisGeneFlk(args)
