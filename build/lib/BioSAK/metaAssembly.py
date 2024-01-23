import os
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

    tmp_dir = '%s/tmp'        % op_dir
    op_txt  = '%s/output.txt' % op_dir
    cmd_txt = '%s/cmds.txt'   % op_dir

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
        op_metadata_tsv   = '%s/%s_data_report.tsv' % (tmp_dir, assembly_id)
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

    print('Parsing outputs from datasets and dataformat')
    missing_output_id_set = set()
    for assembly_id in assembly_id_set:
        op_metadata_tsv = '%s/%s_data_report.tsv' % (tmp_dir, assembly_id)
        if os.path.isfile(op_metadata_tsv) is False:
            missing_output_id_set.add(assembly_id)
        else:
            host = 'na'
            geo_loc_name = 'na'
            biosample_model = 'na'
            isolation_source = 'na'
            metagenome_source = 'na'
            biosample_accession = 'na'
            biosample_description_title = 'na'
            biosample_description_comment = 'na'
            col_index = dict()
            line_num_index = 0
            for each_line in open(op_metadata_tsv):
                line_num_index += 1
                line_split = each_line.strip().split('\t')
                if line_num_index == 1:
                    col_index = {key: i for i, key in enumerate(line_split)}
                else:
                    if len(line_split) == len(col_index):
                        if line_split[col_index['Assembly BioSample Accession']] != '':
                            biosample_accession = line_split[col_index['Assembly BioSample Accession']]
                        if line_split[col_index['Assembly BioSample Description Title']] != '':
                            biosample_description_title = line_split[col_index['Assembly BioSample Description Title']]
                        if line_split[col_index['Assembly BioSample Description Comment']] != '':
                            biosample_description_comment = line_split[col_index['Assembly BioSample Description Comment']]
                        if line_split[col_index['Assembly BioSample Models']] != '':
                            biosample_model = line_split[col_index['Assembly BioSample Models']]

                        biosample_attribute_name = line_split[col_index['Assembly BioSample Attribute Name']]
                        biosample_attribute_value = line_split[col_index['Assembly BioSample Attribute Value']]
                        if biosample_attribute_name == 'isolation_source':
                            isolation_source = biosample_attribute_value
                        elif biosample_attribute_name == 'metagenome_source':
                            metagenome_source = biosample_attribute_value
                        elif biosample_attribute_name == 'host':
                            host = biosample_attribute_value
                        elif biosample_attribute_name == 'geo_loc_name':
                            geo_loc_name = biosample_attribute_value
                    else:
                        missing_output_id_set.add(assembly_id)

            str_to_write = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (assembly_id, biosample_accession, biosample_description_title, biosample_model, isolation_source, metagenome_source, host, geo_loc_name, biosample_description_comment)

            with open(op_txt, 'a') as op_txt_handle:
                op_txt_handle.write(str_to_write + '\n')

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
