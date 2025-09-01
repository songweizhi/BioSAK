import os
import glob
import argparse
import subprocess
import multiprocessing as mp


metaAssembly_usage = '''
==================== metaAssembly example commands ====================

BioSAK metaAssembly -i GCA_947846245.1 -o output_dir -f
BioSAK metaAssembly -i assembly_id.txt -o output_dir -f -t 12

# format of assembly_id.txt (one id per line)
GCA_947846245.1
GCF_026914265.1
GCF_002022765.2

Dependencies: datasets and dataformat
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

=======================================================================
'''


def check_executables(program_list):

    not_detected_programs = []
    for needed_program in program_list:

        if subprocess.call(['which', needed_program], stdout=open(os.devnull, 'wb')) != 0:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not detected, program exited!' % ','.join(not_detected_programs))
        exit()


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
            biosample_attribute_name_list = col_dict.get('Assembly BioSample Attribute Name', [])
            biosample_attribute_value_list = col_dict.get('Assembly BioSample Attribute Value', [])
            for (name, value) in zip(biosample_attribute_name_list, biosample_attribute_value_list):
                file_out_handle.write('Assembly BioSample Attribute %s:\t%s\n' % (name, value))
    file_out_handle.close()


def summarise_metadata(file_dir, file_ext, summary_txt):

    file_re = '%s/*.%s' % (file_dir, file_ext)
    file_list = glob.glob(file_re)

    all_col_name_set = set()
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
                bioSample_attribute_name = line_split[col_index['Assembly BioSample Attribute Name']]
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


def metaAssembly(args):

    assembly_id_txt = args['i']
    op_dir          = args['o']
    num_threads     = args['t']
    force_overwrite = args['f']

    check_executables(['datasets', 'dataformat'])

    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()

    tmp_dir             = '%s/tmp'          % op_dir
    tmp_dir_reformatted = '%s/report'       % op_dir
    op_txt              = '%s/summary.txt'  % op_dir
    cmd_txt             = '%s/cmds.txt'     % op_dir

    os.mkdir(op_dir)
    os.mkdir(tmp_dir)

    assembly_id_set = set()
    if os.path.isfile(assembly_id_txt):
        for assembly_id in open(assembly_id_txt):
            assembly_id_set.add(assembly_id.strip().split()[0])
        if len(assembly_id_set) == 0:
            print('No id found in %s, program exited!' % assembly_id_txt)
            exit()
    else:
        assembly_id_set.add(assembly_id_txt)

    # prepare cmds for datasets and dataformat
    cmd_list = []
    cmd_txt_handle = open(cmd_txt, 'w')
    for assembly_id in assembly_id_set:
        op_metadata_tsv   = '%s/%s.txt' % (tmp_dir, assembly_id)
        ncbi_datasets_cmd = 'datasets summary genome accession %s --as-json-lines | dataformat tsv genome > %s' % (assembly_id, op_metadata_tsv)
        cmd_list.append(ncbi_datasets_cmd)
        cmd_txt_handle.write(ncbi_datasets_cmd + '\n')
    cmd_txt_handle.close()

    # run datasets and dataformat with mp
    print('Running %s commands with %s cores' % (len(cmd_list), num_threads))
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, cmd_list)
    pool.close()
    pool.join()

    # parse outputs from datasets and dataformat
    with open(op_txt, 'w') as op_txt_handle:
        op_txt_handle.write('Assembly\tBioSample\tDescription\tModel\tIsolation_Source\tMetagenome_Source\tHost\tLocation\tComment\n')

    #os.mkdir(tmp_dir_reformatted)

    print('Parsing outputs from datasets and dataformat')
    missing_output_id_set = set()
    for assembly_id in sorted(list(assembly_id_set)):
        op_metadata_tsv             = '%s/%s.txt' % (tmp_dir, assembly_id)
        op_metadata_tsv_reformatted = '%s/%s.txt' % (tmp_dir_reformatted, assembly_id)
        if os.path.isfile(op_metadata_tsv) is False:
            missing_output_id_set.add(assembly_id)

    summarise_metadata(tmp_dir, 'txt', op_txt)

    if len(missing_output_id_set) > 0:
        print('datasets and dataformat failed on the following IDs:')
        print(','.join(missing_output_id_set))

    print('Done!')


if __name__ == '__main__':

    metaAssembly_parser = argparse.ArgumentParser(usage=metaAssembly_usage)
    metaAssembly_parser.add_argument('-i', required=True,                         help='file contains assembly ids')
    metaAssembly_parser.add_argument('-o', required=True,                         help='output directory')
    metaAssembly_parser.add_argument('-f', required=False, action="store_true",   help='force overwrite')
    metaAssembly_parser.add_argument('-t', required=False, type=int, default=1,   help='number of cores, default: 1')
    args = vars(metaAssembly_parser.parse_args())
    metaAssembly(args)
