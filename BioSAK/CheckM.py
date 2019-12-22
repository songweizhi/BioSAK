#!/usr/bin/env python3
import os
import glob
import shutil
import argparse


CheckM_Runner_usage = '''
=================================== CheckM_Runner example commands ===================================

# Example command
BioSAK CheckM_Runner -i bin_folder -x fasta -qsub
BioSAK CheckM_Runner -i bin_folder -x fasta -js_header checkm_js_header.sh -qsub (to be added)

======================================================================================================
'''


CheckM_output_parser_usage = '''
================================ CheckM_output_parser example commands ===============================

# only parse quality file
BioSAK CheckM_op_parser -i combined_qualities.txt -o bin_qualities.txt

# get the quality of qualified bins 
BioSAK CheckM_op_parser -i combined_qualities.txt -complete 99 -o bin_qualities_complete99.txt
BioSAK CheckM_op_parser -i combined_qualities.txt -complete 99 -contain 5 -o bin_qualities_complete99_contain5.txt

# get the quality of qualified bins and copy them into a separate folder
BioSAK CheckM_op_parser -i combined_qualities.txt -bin bin_folder -x fasta -complete 99 -o bin_qualities_complete99.txt
BioSAK CheckM_op_parser -i combined_qualities.txt -bin bin_folder -x fasta -complete 99 -contain 5 -o bin_qualities_complete99_contain5.txt

======================================================================================================
'''


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


def CheckM_Runner(args):

    input_bin_folder =      args['i']
    bin_file_extension =    args['x']
    nodes_number =          args['nodes']
    ppn_number =            args['ppn']
    memory =                args['memory']
    walltime_needed =       args['walltime']
    email =                 args['e']

    modules_needed = [args['python'],
                      args['hmmer'],
                      args['pplacer'],
                      args['prodigal']]

    if input_bin_folder[-1] == '/':
        input_bin_folder = input_bin_folder[:-1]

    input_bin_folder_no_path = input_bin_folder
    if '/' in input_bin_folder:
        input_bin_folder_no_path = input_bin_folder.split('/')[-1]

    ########################################################################################################################

    # define folder/file name
    wd = os.getcwd()
    pwd_checkm_wd = '%s/%s_CheckM_wd' % (wd, input_bin_folder_no_path)
    force_create_folder(pwd_checkm_wd)

    # prepare qsub file header and module lines
    header = ''
    module_lines = ''
    if args['qsub'] is True:

        # prepare qsub file header
        line_1 = '#!/bin/bash\n'
        line_2 = '#PBS -l nodes=' + str(nodes_number) + ':ppn=' + str(ppn_number) + '\n'
        line_3 = '#PBS -l mem=' + str(memory) + 'gb\n'
        line_4 = '#PBS -l walltime=' + walltime_needed + '\n'
        line_5 = '#PBS -j oe\n'
        line_6 = ''
        if args['e'] is not None:
            line_6 = '#PBS -M ' + email + '\n'
        line_7 = '#PBS -m ae\n'
        line_8 = 'cd $PBS_O_WORDIR\n'
        # line_8 = 'cd %s\n' % pwd_qsub_files_folder
        header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7 + line_8

        # Prepare module lines
        module_lines = ''
        for module in modules_needed:
            module_lines += 'module load ' + module + '\n'

    # get bin name list
    bin_files = '%s/%s/*.%s' % (wd, input_bin_folder, bin_file_extension)
    bin_list = [os.path.basename(file_name) for file_name in glob.glob(bin_files)]

    if len(bin_list) == 0:
        print('No input bin detected from %s, please double-check.' % input_bin_folder)
        exit()

    # get qsub file and submit it
    pwd_qsub_files_folder = '%s/job_scripts' % (pwd_checkm_wd)
    if args['qsub'] == True:
        os.mkdir(pwd_qsub_files_folder)

    pwd_checkm_cmds_file = '%s/CheckM_cmd_%s.txt' % (pwd_checkm_wd, input_bin_folder)

    checkm_cmds_file = open(pwd_checkm_cmds_file, 'w')
    for each_bin in bin_list:

        bin_name = each_bin[:-(len(bin_file_extension) + 1)]

        # create a folder for current bin
        current_bin_folder = '%s/%s' % (pwd_checkm_wd, bin_name)
        current_bin_checkm_wd = '%s/%s_CheckM_wd' % (pwd_checkm_wd, bin_name)
        os.mkdir(current_bin_folder)
        pwd_bin = '%s/%s/%s' % (wd, input_bin_folder, each_bin)
        os.system('cp %s %s' % (pwd_bin, current_bin_folder))

        checkm_cmd = 'checkm lineage_wf %s %s -f %s/%s_quality.txt -x %s -t %s\n' % (current_bin_folder,
                                                                                     current_bin_checkm_wd,
                                                                                     pwd_checkm_wd,
                                                                                     bin_name,
                                                                                     bin_file_extension,
                                                                                     ppn_number)
        checkm_cmds_file.write(checkm_cmd)

        if args['qsub'] is True:

            pwd_qsub_file = '%s/%s.sh' % (pwd_qsub_files_folder, bin_name)
            qsub_file_handle = open(pwd_qsub_file, 'w')
            qsub_file_handle.write('%s\n%s' % (header, module_lines))
            qsub_file_handle.write(checkm_cmd)
            qsub_file_handle.close()

            # submit job script
            os.chdir(pwd_qsub_files_folder)
            os.system('qsub %s.sh' % bin_name)
            os.chdir(wd)

    checkm_cmds_file.close()


def CheckM_output_parser(args):

    combined_checkm_quality_file =  args['i']
    input_bin_folder =              args['bin']
    bin_file_extension =            args['x']
    completeness_cutoff =           args['complete']
    contamination_cutoff =          args['contain']
    output_file =                   args['o']

    qualitfied_bin_folder = None
    if input_bin_folder is not None:

        if input_bin_folder[-1] == '/':
            input_bin_folder = input_bin_folder[:-1]

        input_bin_folder_no_path = input_bin_folder
        if '/' in input_bin_folder:
            input_bin_folder_no_path = input_bin_folder.split('/')[-1]

        if (completeness_cutoff is not None) and (contamination_cutoff is None):
            qualitfied_bin_folder = '%s_complete%s' % (input_bin_folder_no_path, completeness_cutoff)

        if (completeness_cutoff is None) and (contamination_cutoff is not None):
            qualitfied_bin_folder = '%s_contain%s' % (input_bin_folder_no_path, contamination_cutoff)

        if (completeness_cutoff is not None) and (contamination_cutoff is not None):
            qualitfied_bin_folder = '%s_complete%s_contain%s' % (input_bin_folder_no_path, completeness_cutoff, contamination_cutoff)

    if (input_bin_folder is not None) and ((completeness_cutoff is not None) or (contamination_cutoff is not None)):
        force_create_folder(qualitfied_bin_folder)

    qualified_bin_num = 0
    output_file_handle = open(output_file, 'w')
    output_file_handle.write('Genome,Completeness,Contamination\n')
    for each_line in open(combined_checkm_quality_file):

        if (not each_line.startswith('-')) and (not each_line.startswith('  Bin')):

            quality_split = each_line.strip().split(' ')
            quality_split_new = []
            for each_line in quality_split:
                if each_line != '':
                    quality_split_new.append(each_line)

            genome_id = quality_split_new[0]
            completeness = float(quality_split_new[12])
            contamination = float(quality_split_new[13])

            if (completeness_cutoff is None) and (contamination_cutoff is None):
                output_file_handle.write('%s,%s,%s\n' % (quality_split_new[0], completeness, contamination))
                qualified_bin_num += 1

            if (completeness_cutoff is not None) and (contamination_cutoff is None):
                if completeness >= completeness_cutoff:
                    output_file_handle.write('%s,%s,%s\n' % (quality_split_new[0], completeness, contamination))
                    qualified_bin_num += 1

                    if (input_bin_folder is not None):
                        os.system('cp %s/%s.%s %s/' % (input_bin_folder, genome_id, bin_file_extension, qualitfied_bin_folder))

            if (completeness_cutoff is None) and (contamination_cutoff is not None):
                if contamination <= contamination_cutoff:
                    output_file_handle.write('%s,%s,%s\n' % (quality_split_new[0], completeness, contamination))
                    qualified_bin_num += 1

                    if (input_bin_folder is not None):
                        os.system('cp %s/%s.%s %s/' % (input_bin_folder, genome_id, bin_file_extension, qualitfied_bin_folder))

            if (completeness_cutoff is not None) and (contamination_cutoff is not None):
                if (completeness >= completeness_cutoff) and (contamination <= contamination_cutoff):
                    output_file_handle.write('%s,%s,%s\n' % (quality_split_new[0], completeness, contamination))
                    qualified_bin_num += 1

                    if (input_bin_folder is not None):
                        os.system('cp %s/%s.%s %s/' % (input_bin_folder, genome_id, bin_file_extension, qualitfied_bin_folder))

    output_file_handle.close()

    # report
    if (completeness_cutoff is not None) or (contamination_cutoff is not None):
        print('The number of qualified genomes: %s' % qualified_bin_num)

