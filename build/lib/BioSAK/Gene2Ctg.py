import os
from Bio import SeqIO


# input files
input_gbk_file = '/Users/songweizhi/Desktop/MT_contigs_lt500bp_prodigal_gbk/MT9_contig_lt500.gbk'
ctg_depth_file = '/Users/songweizhi/Desktop/MT_contigs_lt500bp_prodigal_gbk/MT9_contig_all_depth.txt'


# prepare file name
gbk_file_path, gbk_file_name = os.path.split(input_gbk_file)
gbk_file_name_no_extension, gbk_file_extension = os.path.splitext(gbk_file_name)
output_file = '%s/%s_gene_depth.txt' % (gbk_file_path, gbk_file_name_no_extension)


# get ctg_depth_dict
ctg_depth_dict = {}
for ctg in open(ctg_depth_file):
    if not ctg.startswith('contigName'):
        ctg_split = ctg.strip().split('\t')
        ctg_id = ctg_split[0]
        ctg_depth = ctg_split[2]
        ctg_depth_dict[ctg_id] = ctg_depth


output_file_handle = open(output_file, 'w')
output_file_handle.write('Gene\tContig\tCtg_length\tCtg_depth\n')
for seq_record in SeqIO.parse(input_gbk_file, 'genbank'):

    seq_id = seq_record.id
    seq_len = len(seq_record.seq)
    seq_depth = 'NA'

    if seq_id in ctg_depth_dict:
        seq_depth = ctg_depth_dict[seq_id]

    for feature in seq_record.features:
        feature_id = feature.qualifiers['locus_tag'][0]
        output_file_handle.write('%s\t%s\t%s\t%s\n' % (feature_id, seq_id, seq_len, seq_depth))

output_file_handle.close()








