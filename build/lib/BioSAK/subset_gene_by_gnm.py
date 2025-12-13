import os
import glob
import argparse
from Bio import SeqIO


subset_gene_by_gnm_usage = '''
================ subset_gene_by_gnm example commands ================

BioSAK subset_gene_by_gnm -i dir_in -x da -g gnm.txt -o dir_out -f 

=====================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def subset_gene_by_gnm(args):

    file_dir        = args['i']
    file_ext        = args['x']
    gnm_txt         = args['g']
    op_dir          = args['o']
    force_overwrite = args['f']

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    gnm_set = set()
    for each in open(gnm_txt):
        gnm_set.add(each.strip())

    file_re = '%s/*.%s' % (file_dir, file_ext)
    file_list = glob.glob(file_re)

    for each in file_list:
        f_name, f_path, f_base, f_ext = sep_path_basename_ext(each)
        op_fa = '%s/%s' % (op_dir, f_name)
        op_fa_handle = open(op_fa, 'w')
        for each_seq in SeqIO.parse(each, 'fasta'):
            seq_id = each_seq.id
            gnm_id = '_'.join(each_seq.id.split('_')[:-1])
            if gnm_id in gnm_set:
                op_fa_handle.write('>%s\n' % seq_id)
                op_fa_handle.write('%s\n' % str(each_seq.seq))
            else:
                print(seq_id)
        op_fa_handle.close()


if __name__ == '__main__':

    subset_gene_by_gnm_parser = argparse.ArgumentParser(usage=subset_gene_by_gnm_usage)
    subset_gene_by_gnm_parser.add_argument('-i',  required=True,                          help='input folder')
    subset_gene_by_gnm_parser.add_argument('-x',  required=False, default=None,           help='file extension')
    subset_gene_by_gnm_parser.add_argument('-g',  required=False, default=None,           help='genome id')
    subset_gene_by_gnm_parser.add_argument('-o',  required=True,                          help='output folder')
    subset_gene_by_gnm_parser.add_argument('-f',  required=False, action="store_true",    help='force overwrite')
    args = vars(subset_gene_by_gnm_parser.parse_args())
    subset_gene_by_gnm(args)
