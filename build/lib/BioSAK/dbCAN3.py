import os
import glob
import argparse
import multiprocessing as mp
from datetime import datetime


dbCAN3_parser_usage = '''
============================================ dbCAN3 example commands ============================================

BioSAK dbCAN3 -f -t 6 -db /project/spongeholobiont/DB/dbcan3 -i recipient.faa -o op_dir
BioSAK dbCAN3 -f -t 6 -db /project/spongeholobiont/DB/dbcan3 -i faa_files -x faa -o op_dir

# Note
1. Only accept amino acid sequences as inputs.
2. Prepare database as described in https://dbcan.readthedocs.io/en/latest/user_guide/database_preparation.html
3. The CAZyDB.fam-activities.txt file need to be in the db folder.

=================================================================================================================
'''

time_format = '[%Y-%m-%d %H:%M:%S] '


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def dbCAN3_worker(argument_list):

    pwd_input_faa   = argument_list[0]
    op_prefix       = argument_list[1]
    output_folder   = argument_list[2]
    db_folder       = argument_list[3]
    dt_value        = argument_list[4]

    if os.path.isdir(output_folder) is False:
        os.system('mkdir %s' % output_folder)

    dbcan3_cmd = 'run_dbcan %s protein --out_dir %s --db_dir %s -dt %s --out_pre %s_' % (pwd_input_faa, output_folder, db_folder, dt_value, op_prefix)
    os.system(dbcan3_cmd)


def get_df(overview_file_list, cazy_desc_dict, op_df_txt):

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
    op_df_txt_handle.write('\t%s\n' % '\t'.join(all_detected_cazy_list_sorted_with_desc))
    for each_gnm in sorted(list(stats_dod.keys())):
        current_gnm_stats_dict = stats_dod[each_gnm]
        value_list = [each_gnm]
        for each_cazy in all_detected_cazy_list_sorted:
            value_list.append(str(current_gnm_stats_dict.get(each_cazy, 0)))
        value_str = '\t'.join(value_list)
        op_df_txt_handle.write(value_str + '\n')
    op_df_txt_handle.close()


def dbCAN3(args):

    file_in             = args['i']
    file_extension      = args['x']
    DB_dir              = args['db']
    num_threads         = args['t']
    op_dir              = args['o']
    force_create_op_dir = args['f']

    cazy_fam_activities_txt = '%s/CAZyDB.fam-activities.txt' % DB_dir
    if os.path.isfile(cazy_fam_activities_txt) is False:
        print('%s not found, program exited!' % cazy_fam_activities_txt)
        exit()

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output directory detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    ################################################## if input is file ################################################

    # if input is file
    if os.path.isfile(file_in) is True:
        print(datetime.now().strftime(time_format) + 'Running dbCAN3 for 1 file with %s cores' % (num_threads))
        _, _, file_in_base, _ = sep_path_basename_ext(file_in)
        dbCAN3_worker([file_in, file_in_base, op_dir, DB_dir, num_threads])

    ################################################ if input is folder ################################################

    # if input is folder
    else:
        # check whether input folder exist
        if os.path.isdir(file_in) is False:
            print(datetime.now().strftime(time_format) + 'input folder not found, program exited!')
            exit()
        else:
            # check whether input genome exist
            input_file_re = '%s/*.%s' % (file_in, file_extension)
            input_file_name_list = [os.path.basename(file_name) for file_name in glob.glob(input_file_re)]

            if len(input_file_name_list) == 0:
                print(datetime.now().strftime(time_format) + 'input file not found, program exited!')
                exit()

            ################################################### define file name ###################################################

            tmp_dir   = '%s/tmp'           % op_dir
            op_df_txt = '%s/dbCAN3_df.txt' % op_dir
            os.system('mkdir %s' % tmp_dir)

            ######################################################### main #########################################################

            print(datetime.now().strftime(time_format) + 'Running dbCAN3 for %s input files with %s cores' % (len(input_file_name_list), num_threads))

            overview_file_list = []
            list_for_multiple_arguments_dbCAN3 = []
            for input_file in input_file_name_list:
                input_file_basename = '.'.join(input_file.split('.')[:-1])
                pwd_input_faa       = '%s/%s'                   % (file_in, input_file)
                current_op_dir      = '%s/%s'                   % (tmp_dir, input_file_basename)
                overview_file       = '%s/%s/%s_overview.txt'   % (tmp_dir, input_file_basename, input_file_basename)
                list_for_multiple_arguments_dbCAN3.append([pwd_input_faa, input_file_basename, current_op_dir, DB_dir, 1])
                overview_file_list.append(overview_file)

            # run COG annotaion files with multiprocessing
            pool = mp.Pool(processes=num_threads)
            pool.map(dbCAN3_worker, list_for_multiple_arguments_dbCAN3)
            pool.close()
            pool.join()

            ######################################################### get dataframe #########################################################

            # read in CAZyDB.fam-activities.txt
            fam_to_activities_dict = {}
            for each_fam in open(cazy_fam_activities_txt):
                each_fam_split = each_fam.strip().split('	  ')
                if len(each_fam_split) == 2:
                    fam_id = each_fam_split[0]
                    fam_activities = each_fam_split[1]
                    fam_to_activities_dict[fam_id] = fam_activities

            # get dataframe
            get_df(overview_file_list, fam_to_activities_dict, op_df_txt)

    print(datetime.now().strftime(time_format) + 'Done!')


if __name__ == '__main__':

    dbCAN3_parser = argparse.ArgumentParser(usage=dbCAN3_parser_usage)
    dbCAN3_parser.add_argument('-i',         required=True,                          help='path to input sequences (in multi-fasta format)')
    dbCAN3_parser.add_argument('-o',         required=True,                          help='output directory')
    dbCAN3_parser.add_argument('-x',         required=False,                         help='file extension')
    dbCAN3_parser.add_argument('-db',        required=True,                          help='db folder')
    dbCAN3_parser.add_argument('-t',         required=False, type=int, default=1,    help='number of threads')
    dbCAN3_parser.add_argument('-f',         required=False, action="store_true",    help='force overwrite')
    args = vars(dbCAN3_parser.parse_args())
    dbCAN3(args)
