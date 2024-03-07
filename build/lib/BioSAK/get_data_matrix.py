import os
import glob

###############################################################

input_dir    = '/Users/songweizhi/Desktop/Examples'
file_ext     = 'txt'
note_txt     = '/Users/songweizhi/Desktop/0.Basic.infor'
op_file      = '/Users/songweizhi/Desktop/op_dm.txt'
skip_1st_row = True  # specify as True or False

###############################################################


def sep_path_basename_ext(file_in):
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)
    return f_path, f_base, f_ext


note_dict = dict()
for each_line in open(note_txt):
    each_line_split = each_line.strip().split('\t')
    note_dict[each_line_split[0]] = each_line_split[1]

id_set = set()
annotation_dod = dict()
for each_file in glob.glob('%s/*.%s' % (input_dir, file_ext)):
    f_path, f_base, f_ext = sep_path_basename_ext(each_file)
    line_num_index = 1
    current_annotation_dict = dict()
    for each_line in open(each_file):
        each_line_split = each_line.strip().split('\t')
        if line_num_index == 1:
            if skip_1st_row is False:
                current_annotation_dict[each_line_split[0]] = each_line_split[1]
                id_set.add(each_line_split[0])
        else:
            current_annotation_dict[each_line_split[0]] = each_line_split[1]
            id_set.add(each_line_split[0])
        line_num_index += 1
    annotation_dod[f_base] = current_annotation_dict

id_list_sorted   = sorted([i for i in id_set])
desc_list_sorted = [note_dict.get(i, 'na') for i in id_list_sorted]
file_list_sorted = sorted(list(annotation_dod.keys()))

op_file_handle = open(op_file, 'w')
op_file_handle.write('Class\tNote\t' + '\t'.join(file_list_sorted) + '\n')
for each_id in id_list_sorted:
    id_desc = note_dict.get(each_id, 'na')
    current_id_value_list = [each_id, id_desc]
    for each_sample in file_list_sorted:
        current_dict = annotation_dod[each_sample]
        current_id_value_list.append(current_dict.get(each_id, 'na'))
    op_file_handle.write('\t'.join(current_id_value_list) + '\n')
op_file_handle.close()
