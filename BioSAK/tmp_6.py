import os


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


for n in range(1, 11):
    blca_op = '/Users/songweizhi/Desktop/SMP/00_fa_files_Usearch_BLCA_Silva/s10_AllSamples_BLCA_classifications/s06_AllSamples_unoise_nc_%s.blca.out' % n
    parse_blca_op(blca_op)
    n += 1

