import os
import glob
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp


gbk2ffn_usage = '''
========= gbk2ffn example commands =========

BioSAK gbk2ffn -i MAG.gbk -o MAG.ffn
BioSAK gbk2ffn -i gbk_dir -o ffn_dir -t 6

============================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def gbk2ffn_single(arg_list):

    gbk_in  = arg_list[0]
    ffn_out = arg_list[1]

    ffn_out_handle = open(ffn_out, 'w')
    for seq_record in SeqIO.parse(gbk_in, 'genbank'):
        sequence_str = str(seq_record.seq)
        for feature in seq_record.features:
            if feature.type == 'CDS':
                feature_loci = feature.location
                feature_loci_start = feature_loci.start
                feature_loci_end = feature_loci.end
                feature_loci_strand = feature_loci.strand
                coding_seq = sequence_str[feature_loci_start:feature_loci_end]
                feature_locus_tag = feature.qualifiers['locus_tag'][0]
                coding_seq_5_3 = coding_seq
                if feature_loci_strand == -1:
                    coding_seq_5_3 = str(Seq(coding_seq).reverse_complement())
                ffn_out_handle.write('>%s\n' % feature_locus_tag)
                ffn_out_handle.write('%s\n' % coding_seq_5_3)
    ffn_out_handle.close()


def gbk2ffn(args):

    gbk_in      = args['i']
    in_ext      = args['ix']
    ffn_out     = args['o']
    out_ext     = args['ox']
    num_threads = args['t']

    if os.path.isfile(gbk_in) is True:
        gbk2ffn_single([gbk_in, ffn_out])

    elif os.path.isdir(gbk_in) is True:

        gbk_file_re   = '%s/*.%s' % (gbk_in, in_ext)
        gbk_file_list = glob.glob(gbk_file_re)

        if len(gbk_file_list) == 0:
            print('No input gbk found, program exited!')
            exit()

        if os.path.isdir(ffn_out) is True:
            print('output folder detected, program exited!')
            exit()
        os.system('mkdir %s' % ffn_out)

        mp_arg_lol = []
        for each_gbk in gbk_file_list:
            gbk_path, gbk_base, gbk_ext = sep_path_basename_ext(each_gbk)
            pwd_ffn_out = '%s/%s.%s' % (ffn_out, gbk_base, out_ext)
            mp_arg_lol.append([each_gbk, pwd_ffn_out])

        pool = mp.Pool(processes=num_threads)
        pool.map(gbk2ffn_single, mp_arg_lol)
        pool.close()
        pool.join()
    else:
        print('No input gbk found, program exited!')
        exit()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',  required=True,                       help='input gbk')
    parser.add_argument('-ix', required=False, default='gbk',       help='input file extension, default: gbk')
    parser.add_argument('-o',  required=True,                       help='output ffn')
    parser.add_argument('-ox', required=False, default='ffn',       help='output file extension, default: ffn')
    parser.add_argument('-t',  required=False, type=int, default=1, help='number of threads')
    args = vars(parser.parse_args())
    gbk2ffn(args)
