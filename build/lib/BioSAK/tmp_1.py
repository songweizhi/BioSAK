import os
import glob


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def get_df(overview_file_list, cazy_fam_activities_txt, op_df_txt, op_df_txt_with_desc):
    # read in CAZyDB.fam-activities.txt
    cazy_desc_dict = {}
    for each_fam in open(cazy_fam_activities_txt):
        each_fam_split = each_fam.strip().split('	  ')
        if len(each_fam_split) == 2:
            fam_id = each_fam_split[0]
            fam_activities = each_fam_split[1]
            cazy_desc_dict[fam_id] = fam_activities

    stats_dod = dict()
    all_detected_cazy_set = set()
    for each_overview_file in overview_file_list:
        _, _, overview_file_base, _ = sep_path_basename_ext(each_overview_file)
        gnm_id = overview_file_base.split('_overview')[0]
        current_stats_dict = dict()
        for each_line in open(each_overview_file):
            if '#ofTools' not in each_line:
                each_line_split = each_line.strip().split('\t')
                cazy_id_by_hmmer = each_line_split[2]
                cazy_id_by_dbcan_sub = each_line_split[3]
                cazy_id_by_diamond = each_line_split[4]
                cazy_id_str = ''
                if cazy_id_by_diamond != '-':
                    cazy_id_str = cazy_id_by_diamond
                elif cazy_id_by_hmmer != '-':
                    cazy_id_str = cazy_id_by_hmmer.split('(')[0]
                elif cazy_id_by_dbcan_sub != '-':
                    if '_e' in cazy_id_by_dbcan_sub:
                        cazy_id_str = cazy_id_by_dbcan_sub.split('_e')[0]
                    else:
                        cazy_id_str = cazy_id_by_dbcan_sub

                if '+' in cazy_id_str:
                    cazy_id_list = cazy_id_str.split('+')
                else:
                    cazy_id_list = [cazy_id_str]

                for each_cazy_id in cazy_id_list:
                    if '.' not in each_cazy_id:
                        if each_cazy_id not in current_stats_dict:
                            current_stats_dict[each_cazy_id] = 0
                        current_stats_dict[each_cazy_id] += 1
                        all_detected_cazy_set.add(each_cazy_id)
        stats_dod[gnm_id] = current_stats_dict

    all_detected_cazy_list_sorted = sorted(list(all_detected_cazy_set))
    all_detected_cazy_list_sorted_with_desc = []
    for i in all_detected_cazy_list_sorted:
        i_desc = cazy_desc_dict.get(i, '')
        if i_desc == '':
            all_detected_cazy_list_sorted_with_desc.append(i)
        else:
            i_desc = i_desc.replace('\t', '_')
            i_desc = i_desc.replace(' ', '_')
            i_desc = i_desc.replace('_/_', '/')
            i_desc = i_desc.replace(';_', ';')
            i_desc = i_desc.replace('_(', '(')
            i_desc = i_desc.replace(');;', ');')
            i_desc = i_desc.replace('._', '.')
            i_desc = i_desc.replace('EC_', 'EC')
            i_desc = i_desc.replace(':_', ':')
            i_desc = i_desc.replace(';_', ';')
            i_desc = i_desc.replace(',_', ',')
            i_desc = i_desc.replace(':', '_')
            i_desc = i_desc.replace('_;', ';')
            i_desc = i_desc.replace('_[', '[')
            i_desc = i_desc.strip()
            if (i_desc.endswith('.')) or (i_desc.endswith(';')):
                i_desc = i_desc[:-1]
            all_detected_cazy_list_sorted_with_desc.append('%s__%s' % (i, i_desc))

    op_df_txt_handle = open(op_df_txt, 'w')
    op_df_txt_handle.write('\t%s\n' % '\t'.join(all_detected_cazy_list_sorted))
    op_df_txt_with_desc_handle = open(op_df_txt_with_desc, 'w')
    op_df_txt_with_desc_handle.write('\t%s\n' % '\t'.join(all_detected_cazy_list_sorted_with_desc))
    for each_gnm in sorted(list(stats_dod.keys())):
        current_gnm_stats_dict = stats_dod[each_gnm]
        value_list = [each_gnm]
        for each_cazy in all_detected_cazy_list_sorted:
            value_list.append(str(current_gnm_stats_dict.get(each_cazy, 0)))
        value_str = '\t'.join(value_list)
        op_df_txt_handle.write(value_str + '\n')
        op_df_txt_with_desc_handle.write(value_str + '\n')
    op_df_txt_handle.close()
    op_df_txt_with_desc_handle.close()



file_dir = '/Users/songweizhi/Desktop/overview_files'
file_ext = 'txt'

file_re = '%s/*.%s' % (file_dir, file_ext)

overview_file_list = glob.glob(file_re)



cazy_fam_activities_txt = '/Users/songweizhi/DB/dbcan3/CAZyDB.fam-activities.txt'
op_df_txt               = '/Users/songweizhi/Desktop/annotation_dbCAN3_df.txt'
op_df_txt_with_desc     = '/Users/songweizhi/Desktop/annotation_dbCAN3_df_desc.txt'

get_df(overview_file_list, cazy_fam_activities_txt, op_df_txt, op_df_txt_with_desc)
