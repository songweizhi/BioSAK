#!/usr/bin/env python3
import os
import glob
import shutil
import argparse


CheckM_usage = '''
====================================== CheckM example commands ======================================

# only parse quality file
BioSAK CheckM -i checkm_op.txt -o MAG_qualities.txt
BioSAK CheckM -i checkm_op.txt -o MAG_qualities.txt -r

# get the quality of qualified bins 
BioSAK CheckM -i checkm_op.txt -cpl 99 -o MAG_qualities_cpl99.txt
BioSAK CheckM -i checkm_op.txt -cpl 99 -ctm 5 -o MAG_qualities_cpl99_ctm5.txt

# get the quality of qualified MAG and copy them into a separate folder
BioSAK CheckM -i checkm_op.txt -g MAG_folder -x fa -cpl 99 -o MAG_qualities_cpl99.txt
BioSAK CheckM -i checkm_op.txt -g MAG_folder -x fa -cpl 99 -ctm 5 -o MAG_qualities_cpl99_ctm5.txt

# Note
adjusted_contamination = contamination x (100 - heterogeneity) / 100

=====================================================================================================
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


def CheckM(args):

    input_quality_file = args['i']
    output_file =        args['o']
    input_bin_folder =   args['g']
    bin_file_extension = args['x']
    cpl_cutoff =         args['cpl']
    ctm_cutoff =         args['ctm']
    recalculate_ctm =    args['r']

    qualitfied_bin_folder = None
    if input_bin_folder is not None:

        if input_bin_folder[-1] == '/':
            input_bin_folder = input_bin_folder[:-1]

        input_bin_folder_no_path = input_bin_folder
        if '/' in input_bin_folder:
            input_bin_folder_no_path = input_bin_folder.split('/')[-1]

        if (cpl_cutoff is not None) and (ctm_cutoff is None):
            qualitfied_bin_folder = '%s_cpl%s' % (input_bin_folder_no_path, cpl_cutoff)

        if (cpl_cutoff is None) and (ctm_cutoff is not None):
            qualitfied_bin_folder = '%s_ctm%s' % (input_bin_folder_no_path, ctm_cutoff)

        if (cpl_cutoff is not None) and (ctm_cutoff is not None):
            qualitfied_bin_folder = '%s_cpl%s_ctm%s' % (input_bin_folder_no_path, cpl_cutoff, ctm_cutoff)

    # create output folder
    if (input_bin_folder is not None) and ((cpl_cutoff is not None) or (ctm_cutoff is not None)):
        if os.path.isdir(qualitfied_bin_folder) is True:
            print('Output MAG folder detected, program exited!')
            print(qualitfied_bin_folder)
            exit()
        os.mkdir(qualitfied_bin_folder)

    qualified_bin_num = 0
    output_file_handle = open(output_file, 'w')
    if recalculate_ctm is False:
        output_file_handle.write('Genome\tCompleteness\tContamination\tHeterogeneity\n')
    else:
        output_file_handle.write('Genome\tCompleteness\tContamination\tHeterogeneity\tAdjusted_contamination\n')
    for each_line in open(input_quality_file):

        if (not each_line.startswith('-')) and (not each_line.startswith('  Bin')):
            quality_split = each_line.strip().split(' ')
            quality_split_new = []
            for each_line in quality_split:
                if each_line != '':
                    quality_split_new.append(each_line)

            genome_id     = quality_split_new[0]
            completeness  = float(quality_split_new[12])
            contamination = float(quality_split_new[13])
            heterogeneity = float(quality_split_new[14])

            if cpl_cutoff is None:
                cpl_cutoff = 0

            if ctm_cutoff is None:
                ctm_cutoff = 9999

            if recalculate_ctm is False:
                if (completeness >= cpl_cutoff) and (contamination <= ctm_cutoff):
                    output_file_handle.write('%s\t%s\t%s\t%s\n' % (genome_id, completeness, contamination, heterogeneity))
                    qualified_bin_num += 1
                    if input_bin_folder is not None:
                        os.system('cp %s/%s.%s %s/' % (input_bin_folder, genome_id, bin_file_extension, qualitfied_bin_folder))
            else:
                ctm_adjusted = float("{0:.2f}".format(contamination*(100-heterogeneity)/100))
                if (completeness >= cpl_cutoff) and (ctm_adjusted <= ctm_cutoff):
                    output_file_handle.write('%s\t%s\t%s\t%s\t%s\n' % (genome_id, completeness, contamination, heterogeneity, ctm_adjusted))
                    qualified_bin_num += 1
                    if input_bin_folder is not None:
                        os.system('cp %s/%s.%s %s/' % (input_bin_folder, genome_id, bin_file_extension, qualitfied_bin_folder))
    output_file_handle.close()

    # report
    if (cpl_cutoff is not 0) or (ctm_cutoff is not 999):
        print('The number of qualified genomes: %s' % qualified_bin_num)


if __name__ == "__main__":

    CheckM_parser = argparse.ArgumentParser(usage=CheckM_usage)
    CheckM_parser.add_argument('-i',    required=True,                       help='CheckM produced quality file')
    CheckM_parser.add_argument('-o',    required=True,                       help='output quality file (reformatted)')
    CheckM_parser.add_argument('-g',    required=False,                      help='MAG folder')
    CheckM_parser.add_argument('-x',    required=False, default='fasta',     help='bin file extension')
    CheckM_parser.add_argument('-cpl',  required=False, type=float,          help='completeness cutoff (0-100)')
    CheckM_parser.add_argument('-ctm',  required=False, type=float,          help='contamination cutoff (0-100)')
    CheckM_parser.add_argument('-r',    required=False, action="store_true", help='recalculate contamination')
    args = vars(CheckM_parser.parse_args())
    CheckM(args)
