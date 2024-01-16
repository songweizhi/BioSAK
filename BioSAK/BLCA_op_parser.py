#!/usr/bin/env python3

import argparse
from BioSAK.global_functions import sep_path_basename_ext


BLCA_op_parser_usage = '''
========== BLCA_op_parser example commands ==========

BioSAK BLCA_op_parser -in OTUs.fasta.blca.out

=====================================================
'''


def BLCA_op_parser(args):

    blca_output = args['in']

    file_in_path, file_in_basename, file_in_ext = sep_path_basename_ext(blca_output)
    output_file_1 = '%s/%s_reformatted_1.txt' % (file_in_path, file_in_basename)
    output_file_2 = '%s/%s_reformatted_2.txt' % (file_in_path, file_in_basename)

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
            taxon_blca_raw_split_1 = taxon_blca_raw.strip().split(':')[1:]
            formatted_taxon_list_with_num = []
            formatted_taxon_list_no_num = []
            for each_str in taxon_blca_raw_split_1:
                each_str_split = each_str.split(';')

                # determine_current_rank
                current_rank = ''
                if each_str_split[-1] == 'phylum':
                    current_rank = 'd'
                elif each_str_split[-1] == 'class':
                    current_rank = 'p'
                elif each_str_split[-1] == 'order':
                    current_rank = 'c'
                elif each_str_split[-1] == 'family':
                    current_rank = 'o'
                elif each_str_split[-1] == 'genus':
                    current_rank = 'f'
                elif each_str_split[-1] == 'species':
                    current_rank = 'g'
                elif each_str_split[-1] == '':
                    current_rank = 's'

                taxon_with_confidence = '%s(%s)' % (each_str_split[0], each_str_split[1][:5])
                taxon_without_confidence = '%s__%s' % (current_rank, each_str_split[0])

                formatted_taxon_list_with_num.append(taxon_with_confidence)
                formatted_taxon_list_no_num.append(taxon_without_confidence)

            formatted_taxon_str_with_num = ';'.join(formatted_taxon_list_with_num)
            formatted_taxon_str_no_num = ';'.join(formatted_taxon_list_no_num)

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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage=BLCA_op_parser_usage)
    parser.add_argument('-in', required=True, help='BLCA output')
    args = vars(parser.parse_args())
    BLCA_op_parser(args)
