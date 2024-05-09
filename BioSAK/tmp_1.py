
file_in = '/Users/songweizhi/Desktop/to_be_p/tmp/GCA_013203245.1_data_report.tsv'
file_out = '/Users/songweizhi/Desktop/to_be_p/tmp/GCA_013203245.1_data_report_fmt.tsv'


def get_nonredundant_table(file_in, file_out, ignore_na):

    col_dict = dict()
    col_list = []
    col_index_dict = dict()
    line_num_index = 0
    for each_line in open(file_in):
        line_num_index += 1
        line_split = each_line.strip().split('\t')
        if line_num_index == 1:
            col_list = line_split
            col_index_dict = {key: i for i, key in enumerate(line_split)}
        else:
            for col in col_list:
                if col not in col_dict:
                    col_dict[col] = [line_split[col_index_dict[col]]]
                else:
                    col_dict[col].append(line_split[col_index_dict[col]])

    # write out
    file_out_handle = open(file_out, 'w')
    for col in col_list:

        value_list_uniq = list(set(col_dict[col]))
        if value_list_uniq == ['']:
            value_list_uniq = ['na']

        if col not in ['Assembly BioSample Attribute Name', 'Assembly BioSample Attribute Value']:
            if ignore_na is True:
                if value_list_uniq != ['na']:
                    file_out_handle.write('%s:\t%s\n' % (col, ','.join(value_list_uniq)))
            else:
                file_out_handle.write('%s:\t%s\n' % (col, ','.join(value_list_uniq)))
        else:
            biosample_attribute_name_list  = col_dict.get('Assembly BioSample Attribute Name', [])
            biosample_attribute_value_list = col_dict.get('Assembly BioSample Attribute Value', [])
            for (name, value) in zip(biosample_attribute_name_list, biosample_attribute_value_list):
                file_out_handle.write('Assembly BioSample Attribute %s:\t%s\n' % (name, value))
    file_out_handle.close()


get_nonredundant_table(file_in, file_out, ignore_na=True)

