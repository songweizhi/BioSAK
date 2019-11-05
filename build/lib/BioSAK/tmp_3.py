


def annotation_normalization(file_in, skip_file_in_header, file_in_value_column, Divisor_value, file_out, file_out_header):

    file_out_handle = open(file_out, 'w')
    file_out_handle.write(file_out_header)
    line_num = 0
    for each_line in open(file_in):

        each_line_split = each_line.strip().split('\t')
        value_str = each_line_split[file_in_value_column - 1]

        if (skip_file_in_header is True and line_num > 0) or (skip_file_in_header is False):
            value_pct = float(value_str)*100/Divisor_value
            each_line_split[file_in_value_column - 1] = str(float("{0:.2f}".format(value_pct)))
            file_out_handle.write('%s\n' % '\t'.join(each_line_split))

        line_num += 1

    file_out_handle.close()



file_in = '/Users/songweizhi/Desktop/BH_ER_050417_Refined_11_ko_stats_A_GeneNumber.txt'
file_out = '/Users/songweizhi/Desktop/BH_ER_050417_Refined_11_ko_stats_A_GeneNumber_pct.txt'
file_out_header = 'KO\tGeneNumber_pct\tDescription\n'
Divisor_value = 10000
skip_file_in_header = True
file_in_value_column = 2



annotation_normalization(file_in, skip_file_in_header, file_in_value_column, Divisor_value, file_out, file_out_header)