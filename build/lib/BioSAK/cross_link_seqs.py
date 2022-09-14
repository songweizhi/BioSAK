import os
import argparse
from Bio import SeqIO
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink

cross_link_seqs_usage = '''
==================== cross_link_seqs example commands ====================

module load blast+
BioSAK cross_link_seqs -p Demo -1 ctg1.fa -2 ctg2.fa -l 500 -i 99.9 
BioSAK cross_link_seqs -p Demo -1 ctg1.fa -2 ctg2.fa -l 500 -i 95 -f SVG

==========================================================================
'''


def cross_link_seqs(args):

    prefix          = args['p']
    query_c         = args['1']
    subject_c       = args['2']
    iden_cutoff     = args['i']
    aln_len_cutoff  = args['l']
    plot_format     = args['f']

    # check plot_format
    if plot_format not in ['PDF', 'SVG', 'EPS', 'PNG']:
        print('please choose plot format from PDF, SVG, EPS and PNG')
        exit()

    blast_result_txt          = '%s_blastn.txt'                 % prefix
    blast_result_txt_filtered = '%s_blastn_iden%s_aln%sbp.txt'  % (prefix, iden_cutoff, aln_len_cutoff)
    output_plot               = '%s_iden%s_aln%sbp.%s'          % (prefix, iden_cutoff, aln_len_cutoff, plot_format)

    ############################## prepare for flanking plot ##############################

    # Run Blast
    command_blast = 'blastn -query %s -subject %s -out %s -evalue 1e-5 -outfmt 6 -task blastn' % (query_c, subject_c, blast_result_txt)
    os.system(command_blast)

    # filter blast result
    aln_len_set = set()
    iden_set = set()
    blast_result_txt_filtered_handle = open(blast_result_txt_filtered, 'w')
    for blast_hit in open(blast_result_txt):
        blast_hit_split = blast_hit.strip().split('\t')
        iden = float(blast_hit_split[2])
        align_len = int(blast_hit_split[3])
        if (iden >= iden_cutoff) and (align_len >= aln_len_cutoff):
            blast_result_txt_filtered_handle.write(blast_hit)
            iden_set.add(iden)
            aln_len_set.add(align_len)
    blast_result_txt_filtered_handle.close()

    if len(aln_len_set) == 0:
        print('No matched region was found with specified identity and alignment length cut-offs, program exited!')
        exit()

    iden_min    = min([i for i in iden_set])
    iden_max    = max([i for i in iden_set])
    aln_len_min = min([i for i in aln_len_set])
    aln_len_max = max([i for i in aln_len_set])

    # read in sequences
    ctg_1_record = SeqIO.read(query_c, "fasta")
    ctg_2_record = SeqIO.read(subject_c, "fasta")

    # create an empty diagram and set track length
    diagram = GenomeDiagram.Diagram("TITLE HERE")
    max_len = max(len(ctg_1_record), len(ctg_2_record))

    # add gene content track for gene1_contig
    ctg_1_name_str = '%s (%sbp)    Identity (%s-%s%s)    Alignment length (%s-%sbp)' % (ctg_1_record.id, len(ctg_1_record), iden_min, iden_max, '%', aln_len_min, aln_len_max)
    contig_1_gene_content_track = diagram.new_track(1, name=ctg_1_name_str, greytrack=True, greytrack_labels=1,
                                                    greytrack_font='Helvetica', greytrack_fontsize=18,
                                                    height=0.35, start=0, end=len(ctg_1_record),
                                                    scale=True, scale_fontsize=12, scale_ticks=1,
                                                    scale_smalltick_interval=10000, scale_largetick_interval=10000)

    # add gene content track for gene2_contig
    ctg_2_name_str = '%s (%sbp)    Identity (%s-%s%s)    Alignment length (%s-%sbp)' % (ctg_2_record.id, len(ctg_2_record), iden_min, iden_max, '%', aln_len_min, aln_len_max)
    ctg_2_name_str = '%s (%sbp)' % (ctg_2_record.id, len(ctg_2_record))
    contig_2_gene_content_track = diagram.new_track(1, name=ctg_2_name_str, greytrack=True, greytrack_labels=1,
                                                    greytrack_font='Helvetica', greytrack_fontsize=18,
                                                    height=0.35, start=0, end=len(ctg_2_record),
                                                    scale=True, scale_fontsize=12, scale_ticks=1,
                                                    scale_smalltick_interval=10000, scale_largetick_interval=10000)

    ####################################### add crosslink from blast results #######################################

    # parse blast results
    for each_line in open(blast_result_txt_filtered):
        each_line_split = each_line.split('\t')
        query        = each_line_split[0]
        identity     = float(each_line_split[2])
        query_start  = int(each_line_split[6])
        query_end    = int(each_line_split[7])
        target_start = int(each_line_split[8])
        target_end   = int(each_line_split[9])

        # use color to reflect identity
        color = colors.linearlyInterpolatedColor(colors.white, colors.red, iden_min, iden_max, identity)

        # if query is contig_1
        if query == ctg_1_record.id:
            link = CrossLink((contig_1_gene_content_track, query_start, query_end),
                             (contig_2_gene_content_track, target_start, target_end),
                             color=color, border=color, flip=False)
            diagram.cross_track_links.append(link)

        # if query is contig_2
        elif query == ctg_2_record.id:
            link = CrossLink((contig_2_gene_content_track, query_start, query_end),
                             (contig_1_gene_content_track, target_start, target_end),
                             color=color, border=color, flip=False)
            diagram.cross_track_links.append(link)

    # Draw and Export
    diagram.draw(format="linear", orientation="landscape", pagesize=(50*cm, 25*cm), fragments=1, start=0, end=max_len)
    diagram.write(output_plot, plot_format)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', required=True,                 help='output prefix')
    parser.add_argument('-1', required=True,                 help='sequence 1, in fasta format')
    parser.add_argument('-2', required=True,                 help='sequence 2, in fasta format')
    parser.add_argument('-i', required=True,  type=float,    help='identity cutoff, 0-100')
    parser.add_argument('-l', required=True,  type=int,      help='alignment length cutoff, in bp')
    parser.add_argument('-f', required=False, default='PDF', help='plot format, choose from PDF, SVG, EPS and PNG, default: PDF')
    args = vars(parser.parse_args())
    cross_link_seqs(args)
