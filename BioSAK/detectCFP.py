#!/usr/bin/env python3

import os
import glob
import argparse
import multiprocessing as mp
from distutils.spawn import find_executable


detectCFP_usage = '''
============================================== detectCFP example commands ==============================================

# modules needed
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
module unload R
module load R/3.6.1
module load cplex/12.9.0-academic  

# example commands
python3 detectCFP.py -p DeepSea -g MAG_files -x fna -hmm keyEnzymes.hmm -k path2hmm.txt -t 6 -gapseq /path/to/gapseq
python3 detectCFP.py -p DeepSea -g MAG_files -x fna -hmm keyEnzymes.hmm -k path2hmm.txt -t 6 -gapseq /path/to/gapseq -faa faa_files

========================================================================================================================
'''


def sep_path_basename_ext(file_in):
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'
    file_basename, file_extension = os.path.splitext(file_name)
    return file_path, file_basename, file_extension


def get_pwys_with_qualified_key_enzymes(mag_to_key_enzymes_dict, pwy_to_key_enzyme_dict, key_enzyme_pct_cutoff):

        mag_to_qualified_pwys_dict = {}
        for mag in mag_to_key_enzymes_dict:
            current_mag_detected_key_enzymes = mag_to_key_enzymes_dict[mag]
            for pathway in pwy_to_key_enzyme_dict:
                current_pathway_key_enzyme = pwy_to_key_enzyme_dict[pathway]
                current_pathway_key_enzyme_qualified = []
                for key_enzyme_group in current_pathway_key_enzyme:
                    key_enzyme_group_detected = set()
                    for key_enzyme in key_enzyme_group:
                        if key_enzyme in current_mag_detected_key_enzymes:
                            key_enzyme_group_detected.add(key_enzyme)

                    if (len(key_enzyme_group_detected) * 100 / len(key_enzyme_group)) >= key_enzyme_pct_cutoff:
                        current_pathway_key_enzyme_qualified.append('yes')
                    else:
                        current_pathway_key_enzyme_qualified.append('no')

                if 'yes' in current_pathway_key_enzyme_qualified:
                    if mag not in mag_to_qualified_pwys_dict:
                        mag_to_qualified_pwys_dict[mag] = [pathway]
                    else:
                        mag_to_qualified_pwys_dict[mag].append(pathway)

        return mag_to_qualified_pwys_dict


def detectCFP(args):

    ###################################################### file in/out #####################################################

    output_prefix                    = args['p']
    genome_folder                    = args['g']
    genome_ext                       = args['x']
    faa_folder                       = args['faa']
    hmm_profiles                     = args['hmm']
    hmm_evalue                       = args['e']
    pwy_to_key_enzyme_file           = args['k']
    key_enzyme_percentage_cutoff     = args['kc']
    pwy_completeness_cutoff          = args['pc']
    num_threads                      = args['t']
    gapseq_exe                       = args['gapseq']

    # define output file name
    genome_folder_basename = genome_folder
    if '/' in genome_folder:
        genome_folder_basename = genome_folder.split('/')[-1]
    if output_prefix == None:
        output_prefix = genome_folder_basename

    detectCFP_wd                                = '%s_detectCFP_wd'                % output_prefix
    prodigal_output_folder                      = '%s/%s_prodigal_output'          % (detectCFP_wd, output_prefix)
    hmmsearch_output_folder                     = '%s/%s_hmmsearch_output'         % (detectCFP_wd, output_prefix)
    gapseq_output_folder                        = '%s/%s_gapseq_output'            % (detectCFP_wd, output_prefix)
    gapseq_output_folder_combined_Pathways      = '%s/combined_Pathways'           % gapseq_output_folder
    cmd_file_prodigal                           = '%s/%s_cmds_prodigal.txt'        % (detectCFP_wd, output_prefix)
    cmd_file_hmmsearch                          = '%s/%s_cmds_hmmsearch.txt'       % (detectCFP_wd, output_prefix)
    cmd_file_gapseq                             = '%s/%s_cmds_gapseq.txt'          % (detectCFP_wd, output_prefix)
    output_df_1                                 = '%s/%s_detected_CFPs.1.txt'      % (detectCFP_wd, output_prefix)
    output_df_2                                 = '%s/%s_detected_CFPs.2.txt'      % (detectCFP_wd, output_prefix)

    # create folder
    if os.path.isdir(detectCFP_wd):
        print('Output folder detected, program exited: %s' % detectCFP_wd)
        exit()
    else:
        os.mkdir(detectCFP_wd)


    ################################################ check dependencies ################################################

    # check whether executables exist
    program_list = ['hmmsearch', gapseq_exe]
    if faa_folder is None:
        program_list.append('prodigal')

    not_detected_programs = []
    for needed_program in program_list:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not detected, program exited!' % ', '.join(not_detected_programs))
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
            argument_list_for_prodigal_worker.append(prodigal_cmd)
        cmd_file_prodigal_handle.close()

        # run prodigal with multiprocessing
        pool = mp.Pool(processes=num_threads)
        pool.map(os.system, argument_list_for_prodigal_worker)
        pool.close()
        pool.join()

        print('Commands for running prodigal exported to %s' % os.path.basename(cmd_file_prodigal))

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
        hmmsearch_cmd = 'hmmsearch -o /dev/null --domtblout %s/%s.tbl -E %s %s %s' % (hmmsearch_output_folder, each_faa, hmm_evalue, hmm_profiles, pwd_faa)
        cmd_file_hmmsearch_handle.write('%s\n' % hmmsearch_cmd)
        argument_list_for_hmmsearch_worker.append(hmmsearch_cmd)
    cmd_file_hmmsearch_handle.close()

    # run hmmsearch with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, argument_list_for_hmmsearch_worker)
    pool.close()
    pool.join()

    print('Commands for running hmmsearch exported to %s' % os.path.basename(cmd_file_hmmsearch))

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

    mag_to_paths_with_qualified_key_enzyme_pct = get_pwys_with_qualified_key_enzymes(detected_key_enzymes_dict, pwy_to_key_enzyme_dict, key_enzyme_percentage_cutoff)

    os.mkdir(gapseq_output_folder)

    print('Running Gapseq for %s genomes with %s cores' % (len(mag_to_paths_with_qualified_key_enzyme_pct), num_threads))

    # prepare argument list for GapSeq worker
    cmd_file_gapseq_handle = open(cmd_file_gapseq, 'w')
    cmd_file_gapseq_handle.write('cd %s\n' % gapseq_output_folder)
    argument_list_for_gapseq_worker = []
    for qualified_genome in mag_to_paths_with_qualified_key_enzyme_pct:

        qualified_pwys = mag_to_paths_with_qualified_key_enzyme_pct[qualified_genome]
        pwd_genome_file = '%s/%s.%s' % (genome_folder, qualified_genome, genome_ext)
        if not pwd_genome_file.startswith('/'):
            pwd_genome_file = '../../%s' % pwd_genome_file

        for qualified_pwy in qualified_pwys:
            gapseq_cmd = '%s find -p %s -b 200 %s > %s-%s-stdout.tbl' % (gapseq_exe, qualified_pwy[1:-1], pwd_genome_file, qualified_genome, qualified_pwy[1:-1],)
            cmd_file_gapseq_handle.write('%s\n' % gapseq_cmd)
            argument_list_for_gapseq_worker.append(gapseq_cmd)
    cmd_file_gapseq_handle.close()

    # run GapSeq with multiprocessing
    os.chdir(gapseq_output_folder)
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, argument_list_for_gapseq_worker)
    pool.close()
    pool.join()
    os.chdir('../../')

    print('Commands for running Gapseq exported to %s' % os.path.basename(cmd_file_gapseq))

    # combine Pathways files
    os.mkdir(gapseq_output_folder_combined_Pathways)
    for each_mag in mag_to_paths_with_qualified_key_enzyme_pct:
        combine_cmd = 'cat %s/%s-*-Pathways.tbl > %s/%s-combined-Pathways.tbl' % (gapseq_output_folder, each_mag, gapseq_output_folder_combined_Pathways, each_mag)
        os.system(combine_cmd)

    # parse GapSeq outputs
    pwy_completeness_dict = {}
    for each_gapseq_output in mag_to_paths_with_qualified_key_enzyme_pct:
        pwd_gapseq_output = '%s/%s-combined-Pathways.tbl' % (gapseq_output_folder_combined_Pathways, each_gapseq_output)
        current_genome_pwy_completeness_dict = {}
        for each_pwy in open(pwd_gapseq_output):
            if not each_pwy.startswith('ID'):
                each_pwy_split = each_pwy.strip().split('\t')
                pwy_id = each_pwy_split[0]
                completeness = each_pwy_split[3]
                current_genome_pwy_completeness_dict[pwy_id] = float(completeness)
        pwy_completeness_dict[each_gapseq_output] = current_genome_pwy_completeness_dict


    ##################################################### write out ####################################################

    print('Preparing data matrix')

    # write out as data matrix
    output_df_1_handle = open(output_df_1, 'w')
    output_df_2_handle = open(output_df_2, 'w')
    header_printed = False
    for genome in sorted(intersect_file_list):
        current_genome_header_list = []
        current_genome_value_list = []
        current_genome_header_list_2 = []
        current_genome_value_list_2 = []
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
            current_genome_header_list_2.append(pathway[1:-1])

            if pathway in pwy_completeness_dict[genome]:
                current_genome_value_list.append(str(pwy_completeness_dict[genome][pathway]))
            else:
                current_genome_value_list.append('na')

            if (pathway in mag_to_paths_with_qualified_key_enzyme_pct[genome]) and (pwy_completeness_dict[genome][pathway] >= pwy_completeness_cutoff):
                current_genome_value_list.append('1')
                current_genome_value_list_2.append('1')
            else:
                current_genome_value_list.append('0')
                current_genome_value_list_2.append('0')

        if header_printed is False:
            output_df_1_handle.write('Genome\t%s\n' % '\t'.join(current_genome_header_list))
            output_df_2_handle.write('Genome\t%s\n' % '\t'.join(current_genome_header_list_2))
            header_printed = True
        output_df_1_handle.write('%s\t%s\n' % (genome, '\t'.join(current_genome_value_list)))
        output_df_2_handle.write('%s\t%s\n' % (genome, '\t'.join(current_genome_value_list_2)))

    output_df_1_handle.close()
    output_df_2_handle.close()


    print('Detections exported to %s and %s' % (os.path.basename(output_df_1), os.path.basename(output_df_2)))

    print('Done!')


######################################################### main #########################################################

if __name__ == '__main__':

    detectCFP_parser = argparse.ArgumentParser(usage=detectCFP_usage)

    # arguments for detectCFP
    detectCFP_parser.add_argument('-p',       required=False, default=None,            help='output prefix')
    detectCFP_parser.add_argument('-g',       required=True,                           help='genome folder')
    detectCFP_parser.add_argument('-x',       required=False, default='fasta',         help='genome file extension, default: fasta')
    detectCFP_parser.add_argument('-faa',     required=False, default=None,            help='faa files, requires Prodigal if not provided')
    detectCFP_parser.add_argument('-hmm',     required=True,                           help='hmm profiles')
    detectCFP_parser.add_argument('-k',       required=True,                           help='pathway to hmm profiles')
    detectCFP_parser.add_argument('-e',       required=False, default='1e-99',         help='evalue cutoff for hmmsearch, default: 1e-99')
    detectCFP_parser.add_argument('-kc',      required=False, type=float, default=50,  help='minimum percentage of key enzymes in a pathway, default: 50')
    detectCFP_parser.add_argument('-pc',      required=False, type=float, default=80,  help='pathway completeness cutoff, default: 80')
    detectCFP_parser.add_argument('-gapseq',  required=False, default='gapseq',        help='path to GapSeq executable file, default: gapseq')
    detectCFP_parser.add_argument('-t',       required=False, type=int,   default=1,   help='number of threads, default: 1')

    args = vars(detectCFP_parser.parse_args())

    detectCFP(args)

'''

Note:
1. HMM id in the HMM file and the path2hmm file need to be consistent. 
2. detectCFP.py will produce the faa files with Prodigal if they are not provided. if you already have them, you can specify with "-faa"
2. You need to add your own pathways to gapseq's dat/custom_pwy.tbl file.
3. The format of path2hmm file:
    1) file header need to be "PWY	HMM", separated by tab.
    2) pathway id and hmm id(s) are separated by tab. if there are multiple key enzymes in a pathway, hmm ids need to be separated by comma.
4. final output 1:
    1) each pathway in the final output has three columns: PWY_HMM, PWY_completeness and PWY_found
    2) PWY_HMM: "_n_" refers to "and " and "_v_" refers to "or ".
    3) PWY_completeness: Gapseq provided completeness
    4) PWY_found: "1" for detected and "0" for not.
5. final output 2: presence/absence matrices of interested pathways among MAGs.


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
module unload R
module load R/3.6.1
module load cplex/12.9.0-academic  
cd /srv/scratch/z5039045/detectCFP_wd
python3 detectCFP.py -p Test -g mag_files -x fna -hmm combined.HMM -k pathwaysXhmmfiles.txt -t 6 -gapseq /srv/scratch/z5039045/Softwares/gapseq/gapseq -faa faa_files

'''
