#!/usr/bin/env python3

import os
import glob
import shutil
import argparse
import multiprocessing as mp
from distutils.spawn import find_executable


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


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def list1_in_list2(list1, list2):

    all_list1_elements_in_list2 = True
    for each_element in list1:
        if each_element not in list2:
            all_list1_elements_in_list2 = False

    return all_list1_elements_in_list2


def prodigal_worker(argument_list):
    prodigal_cmd = argument_list[0]
    os.system(prodigal_cmd)


def hmmsearch_worker(argument_list):
    hmmsearch_cmd = argument_list[0]
    os.system(hmmsearch_cmd)


def gapseq_worker(argument_list):
    gapseq_cmd = argument_list[0]
    os.system(gapseq_cmd)


def detectCFP(args):

    ###################################################### file in/out #####################################################

    genome_folder                    = args['g']
    genome_ext                       = args['x']
    faa_folder                       = args['faa']
    hmm_profiles                     = args['hmm']
    pwy_to_key_enzyme_file           = args['k']
    pwy_completeness_cutoff          = args['c']
    num_threads                      = args['t']
    gapseq_exe                       = args['gapseq']


    # define output file name
    genome_folder_basename = genome_folder
    if '/' in genome_folder:
        genome_folder_basename = genome_folder.split('/')[-1]

    detectCFP_wd            = '%s_detectCFP_wd'             % genome_folder_basename
    prodigal_output_folder  = '%s/prodigal_output'          % detectCFP_wd
    hmmsearch_output_folder = '%s/hmmsearch_output'         % detectCFP_wd
    gapseq_output_folder    = '%s/gapseq_output'            % detectCFP_wd
    output_df               = '%s/detected_CFP.txt'         % detectCFP_wd

    cmd_file_prodigal       = '%s/cmds_prodigal.txt'        % detectCFP_wd
    cmd_file_hmmsearch      = '%s/cmds_hmmsearch.txt'       % detectCFP_wd
    cmd_file_gapseq         = '%s/cmds_gapseq.txt'          % detectCFP_wd

    # create folder
    force_create_folder(detectCFP_wd)


    ################################################ check dependencies ################################################

    # check whether executables exist
    program_list = ['prodigal', 'hmmsearch', 'gapseq']
    not_detected_programs = []
    for needed_program in program_list:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not detected, program exited!' % ','.join(not_detected_programs))
        exit()


    ################################################# check input files ################################################

    genome_file_re = '%s/*.%s' % (genome_folder, genome_ext)
    genome_file_list = [os.path.basename(file_name) for file_name in glob.glob(genome_file_re)]
    genome_file_list_no_extension = ['.'.join(i.split('.')[:-1]) for i in genome_file_list]

    if faa_folder is None:

        # create folder
        os.mkdir(prodigal_output_folder)

        print('Running prodigal for %s genomes with %s cores' % (len(genome_file_list_no_extension), num_threads))

        # prepare argument list for prodigal worker
        argument_list_for_prodigal_worker = []
        cmd_file_prodigal_handle = open(cmd_file_prodigal, 'w')
        cmd_file_prodigal_handle.write('cd %s\n' % os.getcwd())
        for genome_file in genome_file_list_no_extension:
            pwd_genome_file = '%s/%s.%s' % (genome_folder, genome_file, genome_ext)
            prodigal_cmd = 'prodigal -i %s -o %s/%s.genes -a %s/%s.faa -p meta -q' % (pwd_genome_file, prodigal_output_folder, genome_file, prodigal_output_folder, genome_file)
            cmd_file_prodigal_handle.write('%s\n' % prodigal_cmd)
            argument_list_for_prodigal_worker.append([prodigal_cmd])
        cmd_file_prodigal_handle.close()

        # run prodigal with multiprocessing
        pool = mp.Pool(processes=num_threads)
        pool.map(prodigal_worker, argument_list_for_prodigal_worker)
        pool.close()
        pool.join()

        print('Commands for running prodigal exported to %s' % cmd_file_prodigal)

        intersect_file_list = genome_file_list_no_extension
        faa_folder = prodigal_output_folder

    else:
        # get faa file list
        faa_file_re = '%s/*.faa' % faa_folder
        faa_file_list = [os.path.basename(file_name) for file_name in glob.glob(faa_file_re)]
        faa_file_list_no_extension = ['.'.join(i.split('.')[:-1]) for i in faa_file_list]

        # get file intersection
        intersect_file_list = set(genome_file_list_no_extension).intersection(faa_file_list_no_extension)
        if len(genome_file_list_no_extension) != len(intersect_file_list):
            print('Found %s genomes with provided faa file, will detect pathways on these genomes.' % len(intersect_file_list))

    # read in pathway to key enzyme file
    pwy_to_key_enzyme_dict = {}
    for each_pathway in open(pwy_to_key_enzyme_file):
        if not each_pathway.startswith('PWY	HMM'):
            each_pathway_split = each_pathway.strip().split('\t')
            pwy_id = each_pathway_split[0]
            key_enzyme_list = each_pathway_split[1].split(',')
            if pwy_id not in pwy_to_key_enzyme_dict:
                pwy_to_key_enzyme_dict[pwy_id] = [key_enzyme_list]
            else:
                pwy_to_key_enzyme_dict[pwy_id].append(key_enzyme_list)

    pwys_to_detect = sorted([i for i in pwy_to_key_enzyme_dict])


    ##################################################### hmmsearch ####################################################

    print('Running hmmsearch for %s genomes with %s cores' % (len(genome_file_list_no_extension), num_threads))

    os.mkdir(hmmsearch_output_folder)

    # prepare argument list for hmmsearch worker
    cmd_file_hmmsearch_handle = open(cmd_file_hmmsearch, 'w')
    cmd_file_hmmsearch_handle.write('cd %s\n' % os.getcwd())
    argument_list_for_hmmsearch_worker = []
    for each_faa in intersect_file_list:
        pwd_faa = '%s/%s.faa' % (faa_folder, each_faa)
        hmmsearch_cmd = 'hmmsearch -o /dev/null --domtblout %s/%s.tbl -E 1e-99 %s %s' % (hmmsearch_output_folder, each_faa, hmm_profiles, pwd_faa)
        cmd_file_hmmsearch_handle.write('%s\n' % hmmsearch_cmd)
        argument_list_for_hmmsearch_worker.append([hmmsearch_cmd])
    cmd_file_hmmsearch_handle.close()


    # run hmmsearch with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(hmmsearch_worker, argument_list_for_hmmsearch_worker)
    pool.close()
    pool.join()

    print('Commands for running hmmsearch exported to %s' % cmd_file_hmmsearch)

    # parse hmmsearch outputs
    detected_key_enzymes_dict = {}
    for each_hmm_out in intersect_file_list:
        pwd_hmm_out = '%s/%s.tbl' % (hmmsearch_output_folder, each_hmm_out)
        detected_hmms = set()
        for each_line in open(pwd_hmm_out):
            if not each_line.startswith('#'):
                each_line_split = each_line.strip().split(' ')
                each_line_split_no_empty_element = []
                for eachelement in each_line_split:
                    if eachelement != '':
                        each_line_split_no_empty_element.append(eachelement)
                hmm_id = each_line_split_no_empty_element[3]
                detected_hmms.add(hmm_id)
        detected_key_enzymes_dict[each_hmm_out] = detected_hmms


    ###################################################### Gapseq ######################################################

    os.mkdir(gapseq_output_folder)

    # get the id of genomes for running Gapseq
    genomes_for_running_gapseq = set()
    for each_genome in intersect_file_list:
        detected_key_enzymes = detected_key_enzymes_dict[each_genome]
        at_least_one_detected = False
        for pwy in pwy_to_key_enzyme_dict:
            for key_enzyme_group in pwy_to_key_enzyme_dict[pwy]:
                if list1_in_list2(key_enzyme_group, detected_key_enzymes) is True:
                    at_least_one_detected = True

        if at_least_one_detected is True:
            genomes_for_running_gapseq.add(each_genome)

    print('Running Gapseq for %s genomes with %s cores' % (len(genomes_for_running_gapseq), num_threads))

    # prepare argument list for GapSeq worker
    cmd_file_gapseq_handle = open(cmd_file_gapseq, 'w')
    cmd_file_gapseq_handle.write('cd %s\n' % gapseq_output_folder)
    argument_list_for_gapseq_worker = []
    for qualified_genome in genomes_for_running_gapseq:
        pwd_genome_file = '%s/%s.%s' % (genome_folder, qualified_genome, genome_ext)
        if not pwd_genome_file.startswith('/'):
            pwd_genome_file = '../../%s' % pwd_genome_file

        gapseq_cmd = '%s find -p all -l custom -b 200 %s > %s_gapseq_stdout.txt' % (gapseq_exe, pwd_genome_file, qualified_genome)
        cmd_file_gapseq_handle.write('%s\n' % gapseq_cmd)
        argument_list_for_gapseq_worker.append([gapseq_cmd])
    cmd_file_gapseq_handle.close()

    # run GapSeq with multiprocessing
    os.chdir(gapseq_output_folder)
    pool = mp.Pool(processes=num_threads)
    pool.map(gapseq_worker, argument_list_for_gapseq_worker)
    pool.close()
    pool.join()
    os.chdir('../../')

    print('Commands for running Gapseq exported to %s' % cmd_file_gapseq)

    # parse GapSeq outputs
    pwy_completeness_dict = {}
    for each_gapseq_output in genomes_for_running_gapseq:
        pwd_gapseq_output = '%s/%s-all-Pathways.tbl' % (gapseq_output_folder, each_gapseq_output)
        current_genome_pwy_completeness_dict = {}
        for each_pwy in open(pwd_gapseq_output):
            if not each_pwy.startswith('ID'):
                each_pwy_split = each_pwy.strip().split('\t')
                pwy_id = each_pwy_split[0]
                completeness = each_pwy_split[3]
                current_genome_pwy_completeness_dict[pwy_id] = float(completeness)
        pwy_completeness_dict[each_gapseq_output] = current_genome_pwy_completeness_dict


    ##################################################### write out ####################################################

    print('Writing out detections to file')

    # write out as data matrix
    output_df_handle = open(output_df, 'w')
    header_printed = False
    for genome in sorted(intersect_file_list):
        current_genome_header_list = []
        current_genome_value_list = []
        for pathway in pwys_to_detect:
            current_pathway_key_enzymes = pwy_to_key_enzyme_dict[pathway]
            group_header_top = []
            group_value_top = []
            for key_enzyme_g in current_pathway_key_enzymes:
                group_header = '_n_'.join(key_enzyme_g)
                group_value = []
                for key_enzyme in key_enzyme_g:
                    if key_enzyme in detected_key_enzymes_dict[genome]:
                        group_value.append('1')
                    else:
                        group_value.append('0')

                group_value_top.append('_n_'.join(group_value))
                group_header_top.append(group_header)

            current_genome_header_list.append('%s__%s' % (pathway, '_v_'.join(group_header_top)))
            current_genome_value_list.append('_v_'.join(group_value_top))

            current_genome_header_list.append('%s_completeness' % pathway)
            current_genome_header_list.append('%s_found' % pathway)

            if pathway in pwy_completeness_dict[genome]:
                current_genome_value_list.append(str(pwy_completeness_dict[genome][pathway]))
            else:
                current_genome_value_list.append('0')

            key_enzymes_detected = False
            for enzyme_group in current_pathway_key_enzymes:
                if list1_in_list2(enzyme_group, detected_key_enzymes_dict[genome]) is True:
                    key_enzymes_detected = True

            if (key_enzymes_detected is True) and (pwy_completeness_dict[genome][pathway] >= pwy_completeness_cutoff):
                current_genome_value_list.append('y')
            else:
                current_genome_value_list.append('n')

        if header_printed is False:
            output_df_handle.write('Genome\t%s\n' % '\t'.join(current_genome_header_list))
            header_printed = True
        output_df_handle.write('%s\t%s\n' % (genome, '\t'.join(current_genome_value_list)))

    output_df_handle.close()


    print('Done!')


######################################################### main #########################################################

if __name__ == '__main__':

    detectCFP_parser = argparse.ArgumentParser()

    # arguments for detectCFP
    detectCFP_parser.add_argument('-g',        required=True,                             help='genome folder')
    detectCFP_parser.add_argument('-x',        required=False, default='fasta',           help='genome file extension, default: fasta')
    detectCFP_parser.add_argument('-faa',      required=False, default=None,              help='faa files, requires Prodigal if not provided')
    detectCFP_parser.add_argument('-hmm',      required=True,                             help='hmm profiles')
    detectCFP_parser.add_argument('-k',        required=True,                             help='pathway to hmm file')
    detectCFP_parser.add_argument('-c',        required=False, type=float, default=80,    help='pathway completeness cutoff, default: 80')
    detectCFP_parser.add_argument('-gapseq',   required=False, default='gapseq',          help='path to GapSeq executable file, default: gapseq')
    detectCFP_parser.add_argument('-t',        required=False, type=int,   default=1,     help='number of threads, default: 1')

    args = vars(detectCFP_parser.parse_args())

    detectCFP(args)


'''

module load python/3.7.3
module load perl/5.28.0
module load blast+/2.9.0
module load hmmer/3.3
module load prodigal/2.6.3
module load git/2.22.0
module load bedtools/2.27.1
module load glpk/4.65
module load barrnap/0.9
module load gcc/7.3.0
module load exonerate/2.2.0
module load parallel/20190522
module add R/3.6.1
module load cplex/12.9.0-academic  
cd /srv/scratch/z5039045/detectCFP_wd
# python3 detectCFP.py -g mag_files -x fna -faa faa_files -hmm combined.HMM -k pathwaysXhmmfiles.txt -c 80 -t 6 -gapseq /srv/scratch/z5039045/Softwares/gapseq/gapseq
python3 /srv/scratch/z5287533/scripts/detectCFP.py -g mag_files -x fna -faa faa_files -hmm combined.HMM -k pathwaysXhmmfiles.txt -c 80 -t 6 -gapseq /path/to/gapseq


# for help:
python3 /srv/scratch/z5287533/scripts/detectCFP.py -h

Note:
1. HMM id in the HMM file and the path2hmm file need to be consistent. 
2. detectCFP.py will produce the faa files with Prodigal if there are not provided. if you already have them, you can specify with "-faa"
2. To avoid deleting GapSeq's default pathway database, put and only put your own patways in gapseq's dat/custom_pwy.tbl file and provide "-l custom" to gapseq.
3. Be aware of the format of path2hmm file.
    1) file header need to be "PWY	HMM", separated by tab.
    2) pathway id and hmm id(s) are separated by tab. if there are multiple key enzymes in a pathway, hmm ids need to be separated by comma.
    3) example: /srv/scratch/z5287533/Willis/pathwaysXhmmfiles.txt
4. you can customize the cutoff for pathway completeness with "-c", default is 80
5. final output:
    1) each pathway in the final output has three columns: PWY_HMM, PWY_completeness and PWY_found
    2) PWY_HMM: "_n_" refers to "and " and "_v_" refers to "or ".
    3) PWY_completeness: Gapseq provided completeness
    4) PWY_found: y for detected and n for not.

'''
