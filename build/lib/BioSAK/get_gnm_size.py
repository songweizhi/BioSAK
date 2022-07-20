import os
import glob
import argparse
from Bio import SeqIO


get_gnm_size_parser_usage = '''
============= get_gnm_size example commands ============= 

BioSAK get_gnm_size -i genome_1.fa
BioSAK get_gnm_size -i refined_MAGs -x fasta 
BioSAK get_gnm_size -i refined_MAGs -x fasta -total

=========================================================
'''


def get_gnm_size(argument_list):

    mag_folder =        argument_list['i']
    mag_extension =     argument_list['x']
    get_total_size =    argument_list['total']


    if os.path.isfile(mag_folder) is True:

        total_len = 0
        for each_seq in SeqIO.parse(mag_folder, 'fasta'):
            total_len += len(each_seq.seq)

        # in Kbp
        total_len_Kbp = total_len/1024
        total_len_Kbp = float("{0:.2f}".format(total_len_Kbp))

        # in Mbp
        total_len_Mbp = total_len/(1024*1024)
        total_len_Mbp = float("{0:.2f}".format(total_len_Mbp))

        print('%s\t%s (bp)\t%s (Kbp)\t%s (Mbp)' % (mag_folder, total_len, total_len_Kbp, total_len_Mbp))


    elif os.path.isdir(mag_folder) is True:

        mag_file_re = '%s/*.%s' % (mag_folder, mag_extension)
        mag_file_list = [os.path.basename(file_name) for file_name in glob.glob(mag_file_re)]

        if len(mag_file_list) == 0:
            print('No file detected, program exited')
            exit()

        print('Genome\tlength(bp)\tlength(Kbp)\tlength(Mbp)')
        all_mag_size = 0
        for each_mag in sorted(mag_file_list):
            pwd_mag = '%s/%s' % (mag_folder, each_mag)
            total_len = 0
            for each_seq in SeqIO.parse(pwd_mag, 'fasta'):
                total_len += len(each_seq.seq)
            all_mag_size += total_len

            # in Kbp
            total_len_Kbp = total_len/1024
            total_len_Kbp = float("{0:.2f}".format(total_len_Kbp))

            # in Mbp
            total_len_Mbp = total_len/(1024*1024)
            total_len_Mbp = float("{0:.2f}".format(total_len_Mbp))

            print('%s\t%s\t%s\t%s' % (each_mag, total_len, total_len_Kbp, total_len_Mbp))

        if get_total_size is True:
            get_total_size_Kbp = all_mag_size/1024
            get_total_size_Mbp = all_mag_size/(1024*1024)
            get_total_size_Kbp = float("{0:.2f}".format(get_total_size_Kbp))
            get_total_size_Mbp = float("{0:.2f}".format(get_total_size_Mbp))
            print('Folder: %s, Number of genome: %s, total size: %s (bp)\t%s (Kbp)\t%s (Mbp)' % (mag_folder, len(mag_file_list), all_mag_size, get_total_size_Kbp, get_total_size_Mbp))


if __name__ == '__main__':

    get_gnm_size_parser = argparse.ArgumentParser()
    get_gnm_size_parser.add_argument('-i',      required=True,                          help='MAG file/folder')
    get_gnm_size_parser.add_argument('-x',      required=False,                         help='file extension')
    get_gnm_size_parser.add_argument('-total',  required=False,  action='store_true',   help='get total size')
    args = vars(get_gnm_size_parser.parse_args())
    get_gnm_size(args)
