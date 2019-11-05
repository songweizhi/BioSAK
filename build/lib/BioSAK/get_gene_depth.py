import argparse
from Bio import SeqIO
from MyBioTools.global_functions import sep_path_basename_ext


get_gene_depth_parser_usage = '''
========================== get_gene_depth example commands ==========================

# Get gene depth by contig depth (together with gbk file)
MyBioTools get_gene_depth -gbk test.gbk -ctg_depth contig_depth.txt
MyBioTools get_gene_depth -gbk test.gbk -ctg_depth contig_depth.txt -skip_header

# Contig depth file format (one contig per line, tab separated)
contig_1   30.16
contig_2   26
contig_3   0.726

=====================================================================================
'''

def get_gene_depth(args):

    gbk_file =                  args['gbk']
    ctg_depth_file =            args['ctg_depth']
    id_column =                 args['id_column']
    depth_column =              args['depth_column']
    skip_depth_file_header =    args['skip_header']

    gbk_file_path, gbk_file_basename, gbk_file_extension = sep_path_basename_ext(gbk_file)
    pwd_depth_file = '%s/%s.depth' % (gbk_file_path, gbk_file_basename)

    # read in depth
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

    # get gene depth
    gene_depth_file_handle = open(pwd_depth_file, 'w')
    for seq_record in SeqIO.parse(gbk_file, 'genbank'):

        seq_id = seq_record.id
        seq_depth = ctg_depth_dict[seq_id]

        for feature in seq_record.features:
            if feature.type == 'CDS':
                gene_id = feature.qualifiers['locus_tag'][0]
                for_out = '%s\t%s\n' % (gene_id, seq_depth)
                gene_depth_file_handle.write(for_out)

    gene_depth_file_handle.close()


if __name__ == '__main__':

    get_gene_depth_parser = argparse.ArgumentParser()

    # arguments for get_gene_depth_parser
    get_gene_depth_parser.add_argument('-gbk',          required=True,                       help='gbk file')
    get_gene_depth_parser.add_argument('-ctg_depth',    required=True,                       help='contig depth file')
    get_gene_depth_parser.add_argument('-id_column',    required=False, default=1, type=int, help='contig id column, default is 1')
    get_gene_depth_parser.add_argument('-depth_column', required=False, default=2, type=int, help='contig depth column, default is 2')
    get_gene_depth_parser.add_argument('-skip_header',  required=False, action='store_true', help='skip the 1st line in contig depth file')

    args = vars(get_gene_depth_parser.parse_args())

    get_gene_depth(args)

