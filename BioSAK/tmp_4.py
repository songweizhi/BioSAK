import os
import glob


def summarise_metadata(file_dir, file_ext, summary_txt):


    file_re = '%s/*.%s' % (file_dir, file_ext)
    file_list = glob.glob(file_re)

    all_col_name_set  = set()
    sample_metadata_dod = dict()
    for each_file in file_list:
        file_name = os.path.basename(each_file)
        sample_id = file_name[:-(len(file_ext) + 1)]
        current_sample_metadata_dict = dict()
        col_index = dict()
        line_num_index = 0
        for each_line in open(each_file):
            line_num_index += 1
            line_split = each_line.split('\t')
            if line_num_index == 1:
                col_index = {key: i for i, key in enumerate(line_split)}
            else:
                bioSample_attribute_name  = line_split[col_index['Assembly BioSample Attribute Name']]
                bioSample_attribute_value = line_split[col_index['Assembly BioSample Attribute Value']]
                for each_col in col_index:
                    if each_col not in ['Assembly BioSample Attribute Name', 'Assembly BioSample Attribute Value']:
                        if each_col not in current_sample_metadata_dict:
                            current_sample_metadata_dict[each_col] = set()
                        current_sample_metadata_dict[each_col].add(line_split[col_index[each_col]].strip())
                        all_col_name_set.add(each_col.strip())
                    else:
                        if each_col == 'Assembly BioSample Attribute Name':
                            new_col_name = 'Assembly BioSample Attribute Name - %s' % bioSample_attribute_name
                            if bioSample_attribute_name not in current_sample_metadata_dict:
                                current_sample_metadata_dict[new_col_name] = set()
                            current_sample_metadata_dict[new_col_name].add(bioSample_attribute_value)
                            all_col_name_set.add(new_col_name.strip())
        sample_metadata_dod[sample_id] = current_sample_metadata_dict


    summary_txt_handle = open(summary_txt, 'w')
    summary_txt_handle.write('ID\t%s\n' % '\t'.join(sorted(list(all_col_name_set))))
    for each_sample in sample_metadata_dod:
        metadata_dict = sample_metadata_dod[each_sample]
        value_list = []
        value_list.append(each_sample)
        for each_col in sorted(list(all_col_name_set)):
            col_value = metadata_dict.get(each_col, set())
            col_value_str = ','.join(sorted(list(col_value)))
            if col_value_str == '':
                col_value_str = 'na'
            value_list.append(col_value_str)
        str_to_write = '\t'.join(value_list)
        summary_txt_handle.write(str_to_write + '\n')
    summary_txt_handle.close()


file_dir    = '/Users/songweizhi/Desktop/Sponge_r226/02_AOA_genomes/GCA_GCF_1167_metadata/tmp2'
file_ext    = 'txt'
summary_txt = '/Users/songweizhi/Desktop/Sponge_r226/02_AOA_genomes/GCA_GCF_1167_metadata/tmp2_summary.txt'

