import os
import glob
import argparse
import subprocess
import multiprocessing as mp


blca_usage = '''
============================== blca example commands ==============================

BioSAK blca -f -t 12 -i otu_fa -x fna -o demo -r ssu_r220.fasta -c ssu_r220.tax

# GTDB SSU (BLCA compatible)
-r /Users/songweizhi/DB/BLCA/GTDB_SSU/ssu_all_r220.blca.fa
-c /Users/songweizhi/DB/BLCA/GTDB_SSU/ssu_all_r220.blca.tax

# SILVA SSU (BLCA compatible)
-r /Users/songweizhi/DB/BLCA/SILVA_SSU/SILVA_138.2_SSURef_NR99_tax_silva.blca.fa
-c /Users/songweizhi/DB/BLCA/SILVA_SSU/SILVA_138.2_SSURef_NR99_tax_silva.blca.tax

===================================================================================
'''


def check_executables(program_list):

    not_detected_programs = []
    for needed_program in program_list:

        if subprocess.call(['which', needed_program], stdout=open(os.devnull, 'wb')) != 0:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not detected, program exited!' % ','.join(not_detected_programs))
        exit()


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def parse_blca_op(blca_output):

    blca_op_name, blca_op_path, blca_op_base, blca_op_ext = sep_path_basename_ext(blca_output)
    output_file_1 = '%s/%s.1.txt' % (blca_op_path, blca_op_base)
    output_file_2 = '%s/%s.2.txt' % (blca_op_path, blca_op_base)

    # read in input file
    s16_taxon_blca_dict = {}
    for each_16s_taxon in open(blca_output):
        each_16s_taxon_split = each_16s_taxon.strip().split('\t')
        s16_taxon_blca_dict[each_16s_taxon_split[0]] = each_16s_taxon_split[1]

    taxon_dict_formatted_with_num = {}
    taxon_dict_formatted_no_num = {}
    for each_16s in s16_taxon_blca_dict:
        taxon_blca_raw = s16_taxon_blca_dict[each_16s]
        formatted_taxon_str_with_num = 'Unclassified'
        formatted_taxon_str_no_num = 'Unclassified'
        if taxon_blca_raw != 'Unclassified':
            split2 = taxon_blca_raw.strip().split(';')[:-1]
            formatted_taxon_list_with_num2 = []
            formatted_taxon_list_no_num2 = []
            element_index = 0
            while element_index < len(split2):
                taxon_blca = split2[element_index]
                confidence_value = split2[element_index + 1]
                confidence_value = float("{0:.2f}".format(float(confidence_value)))
                taxon_rank_blca = taxon_blca.split(':')[0]

                taxon_rank_gtdb = ''
                if taxon_rank_blca == 'superkingdom':
                    taxon_rank_gtdb = 'd'
                elif taxon_rank_blca == 'phylum':
                    taxon_rank_gtdb = 'p'
                elif taxon_rank_blca == 'class':
                    taxon_rank_gtdb = 'c'
                elif taxon_rank_blca == 'order':
                    taxon_rank_gtdb = 'o'
                elif taxon_rank_blca == 'family':
                    taxon_rank_gtdb = 'f'
                elif taxon_rank_blca == 'genus':
                    taxon_rank_gtdb = 'g'
                elif taxon_rank_blca == 'species':
                    taxon_rank_gtdb = 's'

                taxon_name_blca  = ':'.join(taxon_blca.split(':')[1:]).replace(':', '')
                current_rank_with_num = '%s__%s(%s)' % (taxon_rank_gtdb, taxon_name_blca, confidence_value)
                current_rank_no_num   = '%s__%s' % (taxon_rank_gtdb, taxon_name_blca)
                formatted_taxon_list_with_num2.append(current_rank_with_num)
                formatted_taxon_list_no_num2.append(current_rank_no_num)
                element_index += 2

            formatted_taxon_str_with_num = ';'.join(formatted_taxon_list_with_num2)
            formatted_taxon_str_no_num = ';'.join(formatted_taxon_list_no_num2)

        formatted_taxon_str_with_numno_space = '_'.join(formatted_taxon_str_with_num.split(' '))
        formatted_taxon_str_no_num_no_space = '_'.join(formatted_taxon_str_no_num.split(' '))

        taxon_dict_formatted_with_num[each_16s] = formatted_taxon_str_with_numno_space
        taxon_dict_formatted_no_num[each_16s] = formatted_taxon_str_no_num_no_space

    output_file_1_handle = open(output_file_1, 'w')
    output_file_2_handle = open(output_file_2, 'w')
    for each_seq in taxon_dict_formatted_with_num:
        output_file_1_handle.write('%s\t%s\n' % (each_seq, taxon_dict_formatted_with_num[each_seq]))
        output_file_2_handle.write('%s\t%s\n' % (each_seq, taxon_dict_formatted_no_num[each_seq]))
    output_file_1_handle.close()
    output_file_2_handle.close()


def blca(args):

    fa_dir              = args['i']
    fa_ext              = args['x']
    output_folder       = args['o']
    ref_seq             = args['r']
    ref_tax             = args['c']
    num_threads         = args['t']
    force_overwrite     = args['f']

    check_executables(['blastn', 'blastdbcmd', 'clustalo', 'muscle'])

    pwd_current_file    = os.path.realpath(__file__)
    current_file_path   = '/'.join(pwd_current_file.split('/')[:-1])
    blca_main_py        = '%s/blca_main.py' % (current_file_path)

    seq_file_re = '%s/*.%s' % (fa_dir, fa_ext)
    seq_file_list = glob.glob(seq_file_re)

    if len(seq_file_list) == 0:
        print('Sequence file not found, program exited')
        exit()

    # create output folder
    if os.path.isdir(output_folder) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % output_folder)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % output_folder)

    # prepare blca commands
    blca_cmd_list = []
    for seq_file in seq_file_list:
        seq_name, seq_path, seq_base, seq_ext = sep_path_basename_ext(seq_file)
        blca_cmd = 'python3 %s -q %s -r %s -p %s -i %s -o %s' % (blca_main_py, ref_seq, ref_tax, 1, seq_file, output_folder)
        print(blca_cmd)
        #os.system(blca_cmd)
        blca_cmd_list.append(blca_cmd)

    # run blca with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, blca_cmd_list)
    pool.close()
    pool.join()

    # parse blca output
    blca_op_re = '%s/*.blca.out' % output_folder
    blca_op_list = glob.glob(blca_op_re)
    for blca_op in blca_op_list:
        parse_blca_op(blca_op)

    print('Done!')


if __name__ == '__main__':

    blca_parser = argparse.ArgumentParser(usage=blca_usage)
    blca_parser.add_argument('-i',   required=True,                              help='path to input sequences (in multi-fasta format)')
    blca_parser.add_argument('-x',   required=False,                             help='file extension')
    blca_parser.add_argument('-o',   required=True,                              help='output directory')
    blca_parser.add_argument('-r',   required=True,                              help='reference sequences')
    blca_parser.add_argument('-c',   required=True,                              help='reference taxonomy')
    blca_parser.add_argument('-t',   required=False, type=int, default=1,        help='number of threads')
    blca_parser.add_argument('-f',   required=False, action="store_true",        help='force overwrite')
    args = vars(blca_parser.parse_args())
    blca(args)
