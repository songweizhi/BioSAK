import os
import glob
import shutil
import argparse
from Bio import SeqIO
from datetime import datetime
from Bio.SeqRecord import SeqRecord
import multiprocessing as mp
from BioSAK.BioSAK_config import config_dict


COG_parser_usage = '''
====================================== COG2003 example commands ======================================

# module needed
module load python/3.7.3
module load perl/5.20.1
module load blast+/2.6.0

# annotate protein sequences
BioSAK COG2003 -m P -t 6 -db_dir /srv/scratch/z5039045/DB/COG_DB -i recipient.faa
BioSAK COG2003 -m P -t 6 -db_dir /srv/scratch/z5039045/DB/COG_DB -i faa_files -x faa

# annotate DNA sequences
BioSAK COG2003 -m N -t 6 -db_dir /srv/scratch/z5039045/DB/COG_DB -i recipient.ffn
BioSAK COG2003 -m N -t 6 -db_dir /srv/scratch/z5039045/DB/COG_DB -i ffn_files -x ffn

# Prepare DB files:
cd $db_dir
wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid.tbl.gz
wget ftp://ftp.ncbi.nih.gov/pub/COG/COG/fun.txt
wget ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog
gunzip Cog_LE.tar.gz
gunzip cddid.tbl.gz

# How it works:
https://github.com/aleimba/bac-genomics-scripts/tree/master/cdd2cog

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


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_ext = os.path.splitext(file_name)

    return file_path, file_basename, file_ext


def dna2aa(dna_file, aa_file):
    query_aa_handle = open(aa_file, 'w')
    for each in SeqIO.parse(dna_file, 'fasta'):
        each_aa = each.seq.translate()
        each_aa_record = SeqRecord(each_aa)
        each_aa_record.id = each.id
        each_aa_record.description = each.description
        SeqIO.write(each_aa_record, query_aa_handle, 'fasta')
    query_aa_handle.close()


def Annotation_COG_worker(argument_list):

    pwd_input_file = argument_list[0]
    pwd_db =         argument_list[1]
    pwd_cdd2cog =    argument_list[2]
    pwd_cddid =      argument_list[3]
    pwd_fun =        argument_list[4]
    pwd_whog =       argument_list[5]
    sequence_type =  argument_list[6]
    output_folder =  argument_list[7]

    input_seq_no_path, input_seq_no_ext, input_seq_ext = sep_path_basename_ext(pwd_input_file)
    current_output_folder = '%s/%s_COG_wd' % (output_folder, input_seq_no_ext)

    rpsblast_output = '%s_COG.tab' % (input_seq_no_ext)
    pwd_rpsblast_output = '%s/%s_COG.tab' % (current_output_folder, input_seq_no_ext)

    force_create_folder(current_output_folder)


    input_seq_aa = ''
    if (sequence_type == 'N') or (sequence_type == 'n'):
        input_seq_aa = '%s_aa.fasta' % input_seq_no_ext
        dna2aa(pwd_input_file, input_seq_aa)
    elif (sequence_type == 'P') or (sequence_type == 'p'):
        input_seq_aa = pwd_input_file
    else:
        print('Specified input sequence type unrecognizable, program exited!')
        exit()

    # run rpsblast
    os.system('rpsblast -query %s -db %s -out %s -evalue 1e-2 -outfmt 6 -num_threads 1' % (input_seq_aa, pwd_db, pwd_rpsblast_output))

    # run cdd2cog.perl
    current_wd = os.getcwd()
    os.chdir(current_output_folder)
    os.system('perl %s -r %s -c %s -f %s -w %s' % (pwd_cdd2cog, rpsblast_output, pwd_cddid, pwd_fun, pwd_whog))
    os.chdir(current_wd)

    # rename output files
    each_op_file_re = '%s/results/*.txt' % current_output_folder
    each_op_file_list = [os.path.basename(file_name) for file_name in glob.glob(each_op_file_re)]
    for each_op_file in each_op_file_list:

        pwd_each_op_file = '%s/results/%s' % (current_output_folder, each_op_file)
        each_op_file_renamed = '%s_%s' % (input_seq_no_ext, each_op_file)
        pwd_each_op_file_renamed = '%s/%s' % (current_output_folder, each_op_file_renamed)

        os.system('mv %s %s' % (pwd_each_op_file, pwd_each_op_file_renamed))

    os.system('rm -r %s/results' % current_output_folder)


def get_COG_annot_df(annotation_dir, stats_level, annotation_df_absolute_num, annotation_df_percentage):

    annotation_dir_re = '%s/*_COG_wd' % annotation_dir
    annotation_folder_list = [os.path.basename(file_name) for file_name in glob.glob(annotation_dir_re)]

    cog_num_dict = {}
    all_identified_cog = set()
    for annotation_folder in annotation_folder_list:

        annotation_folder_basename = annotation_folder.split('_COG_wd')[0]

        pwd_annotation_stats_file = ''
        if stats_level == 'cog_id':
            pwd_annotation_stats_file = '%s/%s/%s_cog_stats.txt' % (annotation_dir, annotation_folder, annotation_folder_basename)
        if stats_level == 'cog_cate':
            pwd_annotation_stats_file = '%s/%s/%s_func_stats.txt' % (annotation_dir, annotation_folder, annotation_folder_basename)

        current_cog_to_num_dict = {}
        for cog in open(pwd_annotation_stats_file):
            ko_split = cog.strip().split('\t')
            current_cog_to_num_dict[ko_split[0]] = int(ko_split[2])
            all_identified_cog.add(ko_split[0])

        cog_num_dict[annotation_folder_basename] = current_cog_to_num_dict

    all_identified_cog_list = sorted([i for i in all_identified_cog])

    annotation_df_absolute_num_handle = open(annotation_df_absolute_num, 'w')
    annotation_df_percentage_handle = open(annotation_df_percentage, 'w')
    annotation_df_absolute_num_handle.write('\t%s\n' % '\t'.join(all_identified_cog_list))
    annotation_df_percentage_handle.write('\t%s\n' % '\t'.join(all_identified_cog_list))
    for annotation_folder in sorted(annotation_folder_list):

        annotation_folder_basename = annotation_folder.split('_COG_wd')[0]
        current_cog_num_dict = cog_num_dict[annotation_folder_basename]

        current_cog_num_list = []
        for identified_cog in all_identified_cog_list:

            # get num list
            identified_cog_num = 0
            if identified_cog in current_cog_num_dict:
                identified_cog_num = current_cog_num_dict[identified_cog]
            current_cog_num_list.append(identified_cog_num)

        # get percentage
        current_cog_num_sum = sum(current_cog_num_list)
        current_cog_num_list_percentage = [float("{0:.2f}".format(i * 100 / current_cog_num_sum)) for i in current_cog_num_list]

        # write out
        annotation_df_absolute_num_handle.write('%s\t%s\n' % (annotation_folder_basename, '\t'.join([str(i) for i in current_cog_num_list])))
        annotation_df_percentage_handle.write('%s\t%s\n' % (annotation_folder_basename, '\t'.join([str(i) for i in current_cog_num_list_percentage])))

    annotation_df_absolute_num_handle.close()
    annotation_df_percentage_handle.close()


def COG2003(args, config_dict):

    file_in =           args['i']
    file_extension =    args['x']
    sequence_type =     args['m']
    DB_dir =            args['db_dir']
    num_threads =       args['t']

    pwd_cdd2cog = config_dict['cdd2cog_perl']

    pwd_whog =      '%s/whog'       % DB_dir
    pwd_cddid =     '%s/cddid.tbl'  % DB_dir
    pwd_db =        '%s/Cog'        % DB_dir
    pwd_fun =       '%s/fun.txt'    % DB_dir

    time_format = '[%Y-%m-%d %H:%M:%S] '


    ########################################### check input is file or folder ##########################################

    # if input is file
    if os.path.isfile(file_in) is True:

        file_in_path, file_in_basename, file_in_ext = sep_path_basename_ext(file_in)

        Annotation_COG_worker([file_in, pwd_db, pwd_cdd2cog, pwd_cddid, pwd_fun, pwd_whog, sequence_type, file_in_path])

    # if input is folder
    else:

        ############################################# check whether file exist #############################################

        # check whether db file exist
        unfound_inputs = []
        for each_input in [pwd_fun, pwd_cddid, pwd_whog]:
            if (not os.path.isfile(each_input)) and (not os.path.isdir(each_input)):
                unfound_inputs.append(each_input)
        if len(unfound_inputs) > 0:
            for each_unfound in unfound_inputs:
                print('%s not found' % each_unfound)
            exit()


        # check whether input folder exist
        if os.path.isdir(file_in) is False:
            print(datetime.now().strftime(time_format) + 'input folder not found, program exited')
            exit()
        else:
            # check whether input genome exist
            input_file_re = '%s/*.%s' % (file_in, file_extension)
            input_file_name_list = [os.path.basename(file_name) for file_name in glob.glob(input_file_re)]

            if len(input_file_name_list) == 0:
                print(datetime.now().strftime(time_format) + 'input file not found, program exited')
                exit()


            ################################################### define file name ###################################################

            if '/' in file_in:
                file_in_folder_name = file_in.split('/')[-1]
            else:
                file_in_folder_name = file_in

            output_folder = '%s_COG_wd' % file_in_folder_name

            # create output folder
            force_create_folder(output_folder)

            annotation_df_cog_id_num =          '%s/%s_cog_id_num.txt'           % (output_folder, file_in_folder_name)
            annotation_df_cog_id_percent =      '%s/%s_cog_id_percentage.txt'    % (output_folder, file_in_folder_name)
            annotation_df_cog_cate_num =        '%s/%s_cog_cate_num.txt'         % (output_folder, file_in_folder_name)
            annotation_df_cog_cate_percent =    '%s/%s_cog_cate_percentage.txt'  % (output_folder, file_in_folder_name)


            ######################################################### main #########################################################

            print(datetime.now().strftime(time_format) + 'Running COG annotation for %s input files with %s cores' % (len(input_file_name_list), num_threads))

            list_for_multiple_arguments_COG = []
            for input_file in input_file_name_list:

                pwd_input_file = '%s/%s' % (file_in, input_file)
                list_for_multiple_arguments_COG.append([pwd_input_file, pwd_db, pwd_cdd2cog, pwd_cddid, pwd_fun, pwd_whog, sequence_type, output_folder])

            # run COG annotaion files with multiprocessing
            pool = mp.Pool(processes=num_threads)
            pool.map(Annotation_COG_worker, list_for_multiple_arguments_COG)
            pool.close()
            pool.join()


            # get dataframe
            get_COG_annot_df(output_folder, 'cog_id', annotation_df_cog_id_num, annotation_df_cog_id_percent)
            get_COG_annot_df(output_folder, 'cog_cate', annotation_df_cog_cate_num, annotation_df_cog_cate_percent)


        ################################################## Final report ####################################################

    print(datetime.now().strftime(time_format) + 'Done!')


if __name__ == '__main__':

    COG_parser = argparse.ArgumentParser()

    # arguments for COG_parser
    COG_parser.add_argument('-i',               required=True,  help='path to input sequences (in multi-fasta format)')
    COG_parser.add_argument('-x',               required=False, help='file extension')
    COG_parser.add_argument('-m',               required=True,  help='The type of input sequences, "N/n" for "nucleotide", "P/p" for "protein"')
    COG_parser.add_argument('-db_dir',          required=True,  help='folder holds Cog, whog, fun.txt and cddid.tbl')
    COG_parser.add_argument('-t',               required=False, type=int, default=1, help='number of threads')

    args = vars(COG_parser.parse_args())

    COG2003(args, config_dict)
