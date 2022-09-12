#!/usr/bin/env python3
import os
import glob
import shutil
import argparse


CheckM_output_parser_usage = '''
================================ CheckM_output_parser example commands ===============================

# only parse quality file
BioSAK CheckM_op_parser -i combined_qualities.txt -o bin_qualities.txt

# get the quality of qualified bins 
BioSAK CheckM_op_parser -i combined_qualities.txt -cpl 99 -o bin_qualities_cpl99.txt
BioSAK CheckM_op_parser -i combined_qualities.txt -cpl 99 -cont 5 -o bin_qualities_cpl99_cont5.txt

# get the quality of qualified bins and copy them into a separate folder
BioSAK CheckM_op_parser -i combined_qualities.txt -bin bin_folder -x fasta -cpl 99 -o bin_qualities_cpl99.txt
BioSAK CheckM_op_parser -i combined_qualities.txt -bin bin_folder -x fasta -cpl 99 -cont 5 -o bin_qualities_cpl99_cont5.txt

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


def CheckM_op_parser(args):

    combined_checkm_quality_file = args['i']
    input_bin_folder =             args['bin']
    bin_file_extension =           args['x']
    completeness_cutoff =          args['cpl']
    contamination_cutoff =         args['cont']
    output_file =                  args['o']

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
    output_file_handle.write('Genome\tCompleteness\tContamination\n')
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
                output_file_handle.write('%s\t%s\t%s\n' % (quality_split_new[0], completeness, contamination))
                qualified_bin_num += 1

            if (completeness_cutoff is not None) and (contamination_cutoff is None):
                if completeness >= completeness_cutoff:
                    output_file_handle.write('%s\t%s\t%s\n' % (quality_split_new[0], completeness, contamination))
                    qualified_bin_num += 1
                    if input_bin_folder is not None:
                        os.system('cp %s/%s.%s %s/' % (input_bin_folder, genome_id, bin_file_extension, qualitfied_bin_folder))

            if (completeness_cutoff is None) and (contamination_cutoff is not None):
                if contamination <= contamination_cutoff:
                    output_file_handle.write('%s\t%s\t%s\n' % (quality_split_new[0], completeness, contamination))
                    qualified_bin_num += 1
                    if input_bin_folder is not None:
                        os.system('cp %s/%s.%s %s/' % (input_bin_folder, genome_id, bin_file_extension, qualitfied_bin_folder))

            if (completeness_cutoff is not None) and (contamination_cutoff is not None):
                if (completeness >= completeness_cutoff) and (contamination <= contamination_cutoff):
                    output_file_handle.write('%s\t%s\t%s\n' % (quality_split_new[0], completeness, contamination))
                    qualified_bin_num += 1
                    if input_bin_folder is not None:
                        os.system('cp %s/%s.%s %s/' % (input_bin_folder, genome_id, bin_file_extension, qualitfied_bin_folder))
    output_file_handle.close()

    # report
    if (completeness_cutoff is not None) or (contamination_cutoff is not None):
        print('The number of qualified genomes: %s' % qualified_bin_num)
