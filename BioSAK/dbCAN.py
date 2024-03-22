import os
import glob
import shutil
import argparse
from Bio import SeqIO
import multiprocessing as mp
from datetime import datetime
from Bio.SeqRecord import SeqRecord


dbCAN_parser_usage = '''
================================================ dbCAN example commands ================================================

# Dependency: hmmer

BioSAK dbCAN -m P -t 6 -db_dir dbCAN_db_v12 -i recipient.faa 
BioSAK dbCAN -m P -t 6 -db_dir dbCAN_db_v12 -i faa_files -x faa

# Prepare DB files (v12):
cd dbCAN_db_v12
wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/hmmscan-parser.sh
wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.08062022.fam-activities.txt -O CAZyDB.fam-activities.txt
wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN-HMMdb-V12.txt -O dbCAN-fam-HMMs.txt
hmmpress dbCAN-fam-HMMs.txt

# How it works:
1. http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/readme.txt
2. The TotalDepth of a CAZy family is calculated by summing up the depth of all genes assigned to it.
3. The percentage of GeneNumber/TotalDepth of genes assigned to a CAZy family is calculated by dividing it 
   by the summation of GeneNumber/TotalDepth of all identified CAZy families. 

# Note!!!
If you run dbCAN for multiple files in a batch manner and want to have their depth info incorporated into the results, 
you need to provide a folder containing individual depth files for each of your input sequence file.
Name of the depth file needs to be exactly the same as its corresponding sequence file, except the extension which is ".depth".

# Depth file format (one gene per line, tab separated)
gene_1	30
gene_2	10.58

========================================================================================================================
'''

time_format = '[%Y-%m-%d %H:%M:%S] '


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


def sep_path_basename_ext(file_in):

    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def AnnotateNorm(file_in, skip_header, value_column, Divisor_value, file_out, file_out_header):

    file_out_handle = open(file_out, 'w')
    file_out_handle.write(file_out_header)
    line_num = 0
    for each_line in open(file_in):

        each_line_split = each_line.strip().split('\t')
        value_str = each_line_split[value_column - 1]

        if (skip_header is True and line_num > 0) or (skip_header is False):

            value_pct = 0
            if Divisor_value != 0:
                value_pct = float(value_str) * 100 / Divisor_value
            each_line_split[value_column - 1] = str(float("{0:.2f}".format(value_pct)))
            file_out_handle.write('%s\n' % '\t'.join(each_line_split))

        line_num += 1

    file_out_handle.close()


def dna2aa(dna_file, aa_file):
    query_aa_handle = open(aa_file, 'w')
    for each in SeqIO.parse(dna_file, 'fasta'):
        each_aa = each.seq.translate()
        each_aa_record = SeqRecord(each_aa)
        each_aa_record.id = each.id
        each_aa_record.description = each.description
        SeqIO.write(each_aa_record, query_aa_handle, 'fasta')
    query_aa_handle.close()


def dbCAN_worker(argument_list):

    pwd_input_file =         argument_list[0]
    pwd_hmmscan_parser =     argument_list[1]
    pwd_dbCAN_fam_HMMs =     argument_list[2]
    sequence_type =          argument_list[3]
    output_folder =          argument_list[4]
    fam_to_activities_dict = argument_list[5]
    depth_file =             argument_list[6]

    input_seq_path, input_seq_no_ext, input_seq_ext = sep_path_basename_ext(pwd_input_file)
    current_output_folder = '%s/%s_dbCAN_wd' % (output_folder, input_seq_no_ext)

    force_create_folder(current_output_folder)

    input_seq_aa = ''
    if (sequence_type == 'N') or (sequence_type == 'n'):
        input_seq_aa = '%s/%s_aa.fasta' % (current_output_folder, input_seq_no_ext)
        dna2aa(pwd_input_file, input_seq_aa)
    elif (sequence_type == 'P') or (sequence_type == 'p'):
        input_seq_aa = pwd_input_file
    else:
        print('Specified input sequence type unrecognizable, program exited!')
        exit()

    hmmscan_cmd         = "hmmscan --domtblout %s/%s.out.dm %s %s > %s/%s.out"                          % (current_output_folder, input_seq_no_ext, pwd_dbCAN_fam_HMMs, input_seq_aa, current_output_folder, input_seq_no_ext)
    hmmscan_parser_cmd  = "sh %s %s/%s.out.dm > %s/%s.out.dm.ps"                                        % (pwd_hmmscan_parser, current_output_folder, input_seq_no_ext, current_output_folder, input_seq_no_ext)
    final_cat_cmd       = "cat %s/%s.out.dm.ps | awk '$5<1e-18&&$10>0.35' > %s/%s.out.dm.ps.stringent"  % (current_output_folder, input_seq_no_ext, current_output_folder, input_seq_no_ext)
    os.system(hmmscan_cmd)
    os.system(hmmscan_parser_cmd)
    os.system(final_cat_cmd)

    ################################### get functional descriptions for query genes ####################################

    pwd_annotation_results =                      '%s/%s_dbCAN.txt'                      % (current_output_folder, input_seq_no_ext)
    pwd_annotation_results_stats_GeneNumber =     '%s/%s_dbCAN_stats_GeneNumber.txt'     % (current_output_folder, input_seq_no_ext)
    pwd_annotation_results_stats_GeneNumber_pct = '%s/%s_dbCAN_stats_GeneNumber_pct.txt' % (current_output_folder, input_seq_no_ext)
    pwd_annotation_results_stats_TotalDepth =     '%s/%s_dbCAN_stats_TotalDepth.txt'     % (current_output_folder, input_seq_no_ext)
    pwd_annotation_results_stats_TotalDepth_pct = '%s/%s_dbCAN_stats_TotalDepth_pct.txt' % (current_output_folder, input_seq_no_ext)

    # read in depth info
    gene_depth_dict = {}
    if depth_file is not None:
        for each_depth in open(depth_file):
            each_depth_split = each_depth.strip().split('\t')
            gene_depth_dict[each_depth_split[0]] = float(each_depth_split[1])

    # get all sequences in input seq file
    query_seq_list = []
    for query_seq in SeqIO.parse(pwd_input_file, 'fasta'):
        query_seq_list.append(query_seq.id)

    # get total number and depth of all genes in one file
    total_depth_for_all_query_genes = 0
    if depth_file is not None:
        for gene in query_seq_list:
            gene_depth = gene_depth_dict[gene]
            total_depth_for_all_query_genes += gene_depth

    # parse hmmscan results
    pwd_annotation_results_handle = open(pwd_annotation_results, 'w')
    pwd_annotation_results_handle.write('Query\tFamily\tActivities\n')
    hmm_to_gene_member_dict = {}
    for hmm_hit in open('%s/%s.out.dm.ps.stringent' % (current_output_folder, input_seq_no_ext)):
        hmm_hit_split = hmm_hit.strip().split('\t')
        query_id = hmm_hit_split[2]
        matched_hmm = hmm_hit_split[0]
        matched_hmm_id = matched_hmm.split('.hmm')[0]

        # get activities
        matched_hmm_activities = 'NA'
        matched_hmm_id_no_underscore = matched_hmm_id
        if '_' in matched_hmm_id_no_underscore:
            matched_hmm_id_no_underscore = matched_hmm_id_no_underscore.split('_')[0]
        if matched_hmm_id_no_underscore in fam_to_activities_dict:
            matched_hmm_activities = fam_to_activities_dict[matched_hmm_id_no_underscore]

        # get hmm_to_num_dict
        if matched_hmm_id not in hmm_to_gene_member_dict:
            hmm_to_gene_member_dict[matched_hmm_id] = [query_id]
        else:
            hmm_to_gene_member_dict[matched_hmm_id].append(query_id)

        # write out
        pwd_annotation_results_handle.write('%s\t%s\t%s\n' % (query_id, matched_hmm, matched_hmm_activities))
    pwd_annotation_results_handle.close()

    #################### get summary of annotation results GeneNumber ####################

    pwd_annotation_results_stats_GeneNumber_handle = open(pwd_annotation_results_stats_GeneNumber, 'w')
    pwd_annotation_results_stats_GeneNumber_handle.write('Family\tGeneNumber\tActivities\n')
    total_GeneNumber_identified = 0
    for each_hmm in hmm_to_gene_member_dict:
        each_hmm_id = each_hmm.split('.hmm')[0]
        each_hmm_GeneNumber = len(hmm_to_gene_member_dict[each_hmm_id])
        each_hmm_activities = 'NA'
        matched_hmm_id_no_underscore = each_hmm_id
        if '_' in matched_hmm_id_no_underscore:
            matched_hmm_id_no_underscore = matched_hmm_id_no_underscore.split('_')[0]
        if matched_hmm_id_no_underscore in fam_to_activities_dict:
            each_hmm_activities = fam_to_activities_dict[matched_hmm_id_no_underscore]
        pwd_annotation_results_stats_GeneNumber_handle.write('%s\t%s\t%s\n' % (each_hmm_id, each_hmm_GeneNumber, each_hmm_activities))
        total_GeneNumber_identified += each_hmm_GeneNumber
    pwd_annotation_results_stats_GeneNumber_handle.close()

    #################### get summary of annotation results GeneNumber pct ####################

    AnnotateNorm(file_in=pwd_annotation_results_stats_GeneNumber, skip_header=True, value_column=2, Divisor_value=total_GeneNumber_identified, file_out=pwd_annotation_results_stats_GeneNumber_pct, file_out_header='Family\tGeneNumber_pct\tActivities\n')

    #################### get summary of annotation results TotalDepth ####################

    if depth_file is not None:
        pwd_annotation_results_stats_TotalDepth_handle = open(pwd_annotation_results_stats_TotalDepth, 'w')
        pwd_annotation_results_stats_TotalDepth_handle.write('Family\tTotalDepth\tActivities\n')
        total_depth_identified = 0
        for each_hmm in hmm_to_gene_member_dict:
            each_hmm_id = each_hmm.split('.hmm')[0]
            each_hmm_TotalDepth = 0
            for each_gene in hmm_to_gene_member_dict[each_hmm_id]:
                each_gene_depth = gene_depth_dict[each_gene]
                each_hmm_TotalDepth += each_gene_depth
            each_hmm_TotalDepth = float("{0:.2f}".format(each_hmm_TotalDepth))

            each_hmm_activities = 'NA'
            matched_hmm_id_no_underscore = each_hmm_id
            if '_' in matched_hmm_id_no_underscore:
                matched_hmm_id_no_underscore = matched_hmm_id_no_underscore.split('_')[0]
            if matched_hmm_id_no_underscore in fam_to_activities_dict:
                each_hmm_activities = fam_to_activities_dict[matched_hmm_id_no_underscore]

            pwd_annotation_results_stats_TotalDepth_handle.write('%s\t%s\t%s\n' % (each_hmm_id, each_hmm_TotalDepth, each_hmm_activities))
            total_depth_identified += each_hmm_TotalDepth
        pwd_annotation_results_stats_TotalDepth_handle.close()

        #################### get summary of annotation results TotalDepth pct ####################

        AnnotateNorm(file_in=pwd_annotation_results_stats_TotalDepth, skip_header=True, value_column=2, Divisor_value=total_depth_identified, file_out=pwd_annotation_results_stats_TotalDepth_pct, file_out_header='Family\tTotalDepth_pct\tActivities\n')


def get_dbCAN_annot_df(annotation_dir, annotation_df_absolute_num, annotation_df_percentage, with_depth):

    annotation_folder_re = '%s/*_dbCAN_wd' % annotation_dir
    annotation_folder_list = [os.path.basename(file_name) for file_name in glob.glob(annotation_folder_re)]

    cazy_family_to_num_dict = {}
    cazy_family_to_num_pct_dict = {}
    all_identified_cazy_family = set()
    for annotation_folder in annotation_folder_list:

        annotation_folder_basename = annotation_folder.split('_dbCAN_wd')[0]

        if with_depth is False:
            pwd_annotation_stats_file =     '%s/%s/%s_dbCAN_stats_GeneNumber.txt'       % (annotation_dir, annotation_folder, annotation_folder_basename)
            pwd_annotation_stats_file_pct = '%s/%s/%s_dbCAN_stats_GeneNumber_pct.txt'   % (annotation_dir, annotation_folder, annotation_folder_basename)
        else:
            pwd_annotation_stats_file =     '%s/%s/%s_dbCAN_stats_TotalDepth.txt'       % (annotation_dir, annotation_folder, annotation_folder_basename)
            pwd_annotation_stats_file_pct = '%s/%s/%s_dbCAN_stats_TotalDepth_pct.txt'   % (annotation_dir, annotation_folder, annotation_folder_basename)

        current_cazy_family_to_num_dict = {}
        for cazy_family in open(pwd_annotation_stats_file):
            if not cazy_family.startswith('Family'):
                cazy_family_split = cazy_family.strip().split('\t')
                if with_depth is False:
                    current_cazy_family_to_num_dict[cazy_family_split[0]] = int(cazy_family_split[1])
                else:
                    current_cazy_family_to_num_dict[cazy_family_split[0]] = float(cazy_family_split[1])
                all_identified_cazy_family.add(cazy_family_split[0])

        current_cazy_family_to_num_pct_dict = {}
        for cazy_family in open(pwd_annotation_stats_file_pct):
            if not cazy_family.startswith('Family'):
                cazy_family_split = cazy_family.strip().split('\t')
                current_cazy_family_to_num_pct_dict[cazy_family_split[0]] = float(cazy_family_split[1])
                all_identified_cazy_family.add(cazy_family_split[0])

        cazy_family_to_num_dict[annotation_folder_basename] = current_cazy_family_to_num_dict
        cazy_family_to_num_pct_dict[annotation_folder_basename] = current_cazy_family_to_num_pct_dict

    all_identified_cazy_family_list = sorted([i for i in all_identified_cazy_family])

    annotation_df_absolute_num_handle = open(annotation_df_absolute_num, 'w')
    annotation_df_percentage_handle = open(annotation_df_percentage, 'w')
    annotation_df_absolute_num_handle.write('\t%s\n' % '\t'.join(all_identified_cazy_family_list))
    annotation_df_percentage_handle.write('\t%s\n' % '\t'.join(all_identified_cazy_family_list))
    for annotation_folder in sorted(annotation_folder_list):
        annotation_folder_basename = annotation_folder.split('_dbCAN_wd')[0]
        current_annotation_cazy_num_dict = cazy_family_to_num_dict[annotation_folder_basename]
        current_annotation_cazy_num_pct_dict = cazy_family_to_num_pct_dict[annotation_folder_basename]

        current_annotation_cazy_num_list = []
        current_annotation_cazy_num_list_pct = []
        for identified_cazy_family in all_identified_cazy_family_list:
            identified_cazy_family_num = 0
            identified_cazy_family_num_pct = 0
            if identified_cazy_family in current_annotation_cazy_num_dict:
                identified_cazy_family_num = current_annotation_cazy_num_dict[identified_cazy_family]
                identified_cazy_family_num_pct = current_annotation_cazy_num_pct_dict[identified_cazy_family]
            current_annotation_cazy_num_list.append(identified_cazy_family_num)
            current_annotation_cazy_num_list_pct.append(identified_cazy_family_num_pct)

        # write out
        annotation_df_absolute_num_handle.write('%s\t%s\n' % (annotation_folder_basename, '\t'.join([str(i) for i in current_annotation_cazy_num_list])))
        annotation_df_percentage_handle.write('%s\t%s\n' % (annotation_folder_basename, '\t'.join([str(i) for i in current_annotation_cazy_num_list_pct])))
    annotation_df_absolute_num_handle.close()
    annotation_df_percentage_handle.close()


def dbCAN(args):

    file_in =           args['i']
    file_extension =    args['x']
    sequence_type =     args['m']
    depth_file =        args['d']
    DB_dir =            args['db_dir']
    num_threads =       args['t']

    pwd_hmmscan_parser =    '%s/hmmscan-parser.sh'            % DB_dir
    pwd_dbCAN_fam_HMMs =    '%s/dbCAN-fam-HMMs.txt'           % DB_dir
    CAZyDB_fam_activities = '%s/CAZyDB.fam-activities.txt'    % DB_dir

    CAZyDB_fam_activities_07312019 = '%s/CAZyDB.07312019.fam-activities.txt'    % DB_dir
    if (os.path.isfile(CAZyDB_fam_activities_07312019) is True) and (os.path.isfile(CAZyDB_fam_activities) is False):
        os.system('mv %s %s' % (CAZyDB_fam_activities_07312019, CAZyDB_fam_activities))

    ############################################ check whether db file exist ###########################################

    # check whether db file exist
    unfound_inputs = []
    for each_input in [pwd_hmmscan_parser, pwd_dbCAN_fam_HMMs]:
        if (not os.path.isfile(each_input)) and (not os.path.isdir(each_input)):
            unfound_inputs.append(each_input)
    if len(unfound_inputs) > 0:
        for each_unfound in unfound_inputs:
            print('%s not found' % each_unfound)
        exit()

    # store CAZyDB.fam-activities.txt in dict
    fam_to_activities_dict = {}
    for each_fam in open(CAZyDB_fam_activities):
        each_fam_split = each_fam.strip().split('	  ')
        if len(each_fam_split) == 2:
            fam_id = each_fam_split[0]
            fam_activities = each_fam_split[1]
            fam_to_activities_dict[fam_id] = fam_activities

    ################################################## if input is file ################################################

    # if input is file
    if os.path.isfile(file_in) is True:

        # check whether depth file exist
        if depth_file is not None:
            if os.path.isfile(depth_file) is False:
                print(datetime.now().strftime(time_format) + 'specified depth file not found, program exited!')
                exit()

        print(datetime.now().strftime(time_format) + 'Running dbCAN for 1 file with %s cores' % (num_threads))

        file_in_path, file_in_basename, file_in_ext = sep_path_basename_ext(file_in)
        dbCAN_worker([file_in, pwd_hmmscan_parser, pwd_dbCAN_fam_HMMs, sequence_type, file_in_path, fam_to_activities_dict, depth_file])

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

            # check whether depth file exist
            if depth_file is not None:

                if os.path.isfile(depth_file) is True:
                    print(datetime.now().strftime(time_format) + 'please provide the folder containing individual depth files (with extension .depth) for each of your input sequence file.')
                    print(datetime.now().strftime(time_format) + 'a single file (not folder) detected, program exited!')
                    exit()

                if os.path.isdir(depth_file) is False:
                    print(datetime.now().strftime(time_format) + 'provided depth folder not found, program exited!')
                    exit()

                if os.path.isdir(depth_file) is True:

                    undetected_depth_file = []
                    for input_seq_file in input_file_name_list:
                        input_seq_file_basename = '.'.join(input_seq_file.split('.')[:-1])
                        input_seq_file_depth = '%s/%s.depth' % (depth_file, input_seq_file_basename)
                        if os.path.isfile(input_seq_file_depth) is False:
                            undetected_depth_file.append(input_seq_file_depth)

                    if len(undetected_depth_file) > 0:
                        print(datetime.now().strftime(time_format) + 'the following depth files not found, program exited!')
                        print(','.join(undetected_depth_file))
                        exit()

            ################################################### define file name ###################################################

            if '/' in file_in:
                file_in_folder_name = file_in.split('/')[-1]
            else:
                file_in_folder_name = file_in

            output_folder                   = '%s_dbCAN_wd'                 % file_in_folder_name
            annotation_df_GeneNumber        = '%s/%s_GeneNumber.txt'        % (output_folder, file_in_folder_name)
            annotation_df_GeneNumber_pct    = '%s/%s_GeneNumber_pct.txt'    % (output_folder, file_in_folder_name)
            annotation_df_TotalDepth        = '%s/%s_TotalDepth.txt'        % (output_folder, file_in_folder_name)
            annotation_df_TotalDepth_pct    = '%s/%s_TotalDepth_pct.txt'    % (output_folder, file_in_folder_name)

            # create output folder
            force_create_folder(output_folder)

            ######################################################### main #########################################################

            print(datetime.now().strftime(time_format) + 'Running dbCAN for %s input files with %s cores' % (len(input_file_name_list), num_threads))

            list_for_multiple_arguments_dbCAN = []
            for input_file in input_file_name_list:

                input_file_basename = '.'.join(input_file.split('.')[:-1])
                pwd_input_file = '%s/%s' % (file_in, input_file)

                # get path to current depth file
                if depth_file is None:
                    input_file_depth = None
                else:
                    input_file_depth = '%s/%s.depth' % (depth_file, input_file_basename)

                list_for_multiple_arguments_dbCAN.append([pwd_input_file, pwd_hmmscan_parser, pwd_dbCAN_fam_HMMs, sequence_type, output_folder, fam_to_activities_dict, input_file_depth])

            # run COG annotaion files with multiprocessing
            pool = mp.Pool(processes=num_threads)
            pool.map(dbCAN_worker, list_for_multiple_arguments_dbCAN)
            pool.close()
            pool.join()

            ######################################################### get dataframe #########################################################

            get_dbCAN_annot_df(output_folder, annotation_df_GeneNumber, annotation_df_GeneNumber_pct, with_depth=False)
            if depth_file is not None:
                get_dbCAN_annot_df(output_folder, annotation_df_TotalDepth, annotation_df_TotalDepth_pct, with_depth=True)

            # report
            print(datetime.now().strftime(time_format) + 'Data matrix exported to:')
            print(datetime.now().strftime(time_format) + annotation_df_GeneNumber.split('/')[-1])
            print(datetime.now().strftime(time_format) + annotation_df_GeneNumber_pct.split('/')[-1])
            if depth_file is not None:
                print(datetime.now().strftime(time_format) + annotation_df_TotalDepth.split('/')[-1])
                print(datetime.now().strftime(time_format) + annotation_df_TotalDepth_pct.split('/')[-1])

    print(datetime.now().strftime(time_format) + 'Done!')


if __name__ == '__main__':

    dbCAN_parser = argparse.ArgumentParser(usage=dbCAN_parser_usage)
    dbCAN_parser.add_argument('-i',         required=True,                          help='path to input sequences (in multi-fasta format)')
    dbCAN_parser.add_argument('-x',         required=False,                         help='file extension')
    dbCAN_parser.add_argument('-m',         required=False, default='P',            help='input sequence type, "N/n" for "nucleotide", "P/p" for "protein"')
    dbCAN_parser.add_argument('-d',         required=False, default=None,           help='gene depth file/folder')
    dbCAN_parser.add_argument('-db_dir',    required=True,                          help='db folder')
    dbCAN_parser.add_argument('-t',         required=False, type=int, default=1,    help='number of threads')
    args = vars(dbCAN_parser.parse_args())
    dbCAN(args)
