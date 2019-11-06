import os
import argparse
from Bio import SeqIO
from datetime import datetime
from BioSAK.global_functions import time_format
from BioSAK.global_functions import sep_path_basename_ext


get_gene_depth_parser_usage = '''
========================== get_gene_depth example commands ==========================

# Get gene depth by contig depth, together with gbk or gff (prefered) file
BioSAK get_gene_depth -gff test.gff -ctg_depth contig_depth.txt
BioSAK get_gene_depth -gff test.gff -ctg_depth contig_depth.txt -skip_header
BioSAK get_gene_depth -gbk test.gbk -ctg_depth contig_depth.txt -skip_header

# Contig depth file format (one contig per line, tab separated)
contig_1    30.16
contig_2    26
contig_3    0.726

=====================================================================================
'''

def get_gene_depth(args):

    gbk_file =                  args['gbk']
    gff_file =                  args['gff']
    ctg_depth_file =            args['ctg_depth']
    id_column =                 args['id_column']
    depth_column =              args['depth_column']
    skip_depth_file_header =    args['skip_header']

    ################################################# check input file #################################################

    annotation_file = None
    if (gbk_file is None) and (gff_file is None):
        print(datetime.now().strftime(time_format) + 'Please provide either a gbk file or a gff file, program exited!')
        exit()

    if (gbk_file is not None) and (gff_file is not None):
        print(datetime.now().strftime(time_format) + 'Both gbk and gff file were provided, will parse annotation results from gff file')
        annotation_file = gff_file

    if (gbk_file is not None) and (gff_file is None):
        annotation_file = gbk_file

    if (gbk_file is None) and (gff_file is not None):
        annotation_file = gff_file

    # define output gene deptp file name
    annotation_file_path, annotation_file_basename, annotation_file_extension = sep_path_basename_ext(annotation_file)
    pwd_gene_depth_file = '%s/%s.depth' % (annotation_file_path, annotation_file_basename)
    if os.path.isfile(pwd_gene_depth_file) is True:
        pwd_gene_depth_file = '%s/%s.depth.txt' % (annotation_file_path, annotation_file_basename)

    ################################################ read in ctg depth #################################################

    ctg_depth_dict = {}
    line = 0
    for ctg in open(ctg_depth_file):
        ctg_split = ctg.strip().split('\t')

        if skip_depth_file_header is True:
            if line > 0:
                ctg_depth_dict[ctg_split[id_column - 1]] = float(ctg_split[depth_column - 1])
        else:
            ctg_depth_dict[ctg_split[id_column - 1]] = float(ctg_split[depth_column - 1])

        line += 1

    ########################################### get gene depth with gbk file ###########################################

    if annotation_file == gbk_file:
        gene_depth_file_handle = open(pwd_gene_depth_file, 'w')
        for seq_record in SeqIO.parse(gbk_file, 'genbank'):
            seq_id = seq_record.id
            seq_depth = ctg_depth_dict[seq_id]
            for feature in seq_record.features:
                if (feature.type != 'source') and (feature.type != 'assembly_gap'):
                    gene_id = feature.qualifiers['locus_tag'][0]
                    for_out = '%s\t%s\n' % (gene_id, seq_depth)
                    gene_depth_file_handle.write(for_out)
        gene_depth_file_handle.close()

    ########################################### get gene depth with gff file ###########################################

    if annotation_file == gff_file:
        gene_depth_file_handle = open(pwd_gene_depth_file, 'w')
        for each_line in open(gff_file):
            if not each_line.startswith('#'):
                each_line_split = each_line.strip().split('\t')
                if len(each_line_split) > 1:
                    seq_id = each_line_split[0]
                    seq_depth = ctg_depth_dict[seq_id]
                    gene_id = each_line_split[8].split(';')[0].split('ID=')[-1]
                    for_out = '%s\t%s\n' % (gene_id, seq_depth)
                    gene_depth_file_handle.write(for_out)
        gene_depth_file_handle.close()

    ###################################################### report ######################################################

    print(datetime.now().strftime(time_format) + 'Gene depth exported to %s' % pwd_gene_depth_file)
    print(datetime.now().strftime(time_format) + 'Done!')


if __name__ == '__main__':

    get_gene_depth_parser = argparse.ArgumentParser()

    # arguments for get_gene_depth_parser
    get_gene_depth_parser.add_argument('-gbk',          required=False, default=None,        help='gbk file')
    get_gene_depth_parser.add_argument('-gff',          required=False, default=None,        help='gff file')
    get_gene_depth_parser.add_argument('-ctg_depth',    required=True,                       help='contig depth file')
    get_gene_depth_parser.add_argument('-id_column',    required=False, default=1, type=int, help='contig id column, default is 1')
    get_gene_depth_parser.add_argument('-depth_column', required=False, default=2, type=int, help='contig depth column, default is 2')
    get_gene_depth_parser.add_argument('-skip_header',  required=False, action='store_true', help='skip the 1st line in contig depth file')

    args = vars(get_gene_depth_parser.parse_args())

    get_gene_depth(args)
