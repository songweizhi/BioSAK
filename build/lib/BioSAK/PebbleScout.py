import os
import glob
import argparse
import pandas as pd
from Bio import SeqIO
import multiprocessing as mp


PebbleScout_usage = '''
=================== PebbleScout example commands ===================

BioSAK PebbleScout -i gnm_dir -x fna -db meta,meta_vol2 -o op_dir

====================================================================
'''


def sep_path_basename_ext(file_in):
    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]
    return f_name, f_path, f_base, f_ext


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):

    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


def pebblescout_results_to_data_matrix(file_base_list, db_list, pebblescout_op_dir, min_cov_pct, op_txt, op_txt_t, metadata_txt):

    subject_to_biosample_dict = dict()
    biosample_title_dict = dict()
    subject_id_set_all = set()
    biosample_subject_id_set_all = set()
    pebblescout_op_dod = dict()
    for f_base in file_base_list:
        current_gnm_dict = dict()
        for each_db in db_list:
            op_file = '%s/%s_%s.txt' % (pebblescout_op_dir, f_base, each_db)
            if os.path.isfile(op_file):
                for each_line in open(op_file):
                    if not each_line.startswith('QueryID\tSubjectID'):
                        each_line_split = each_line.strip().split('\t')
                        subject_id      = each_line_split[1]
                        pct_coverage    = each_line_split[3]
                        bioSample_id    = each_line_split[5]
                        bioSample_title = each_line_split[6]
                        if (bioSample_id != '""') and (bioSample_title != '""'):
                            if float(pct_coverage) >= min_cov_pct:
                                current_gnm_dict[subject_id] = pct_coverage
                                subject_to_biosample_dict[subject_id] = bioSample_id
                                biosample_title_dict[bioSample_id] = bioSample_title
                                subject_id_set_all.add(subject_id)
                                biosample_subject_id_set_all.add('%s_%s' % (bioSample_id, subject_id))
        pebblescout_op_dod[f_base] = current_gnm_dict

    subject_id_set_all_with_biosample_sorted = sorted(list(biosample_subject_id_set_all))
    subject_desc_list_sorted = []
    for each_biosample_subject in subject_id_set_all_with_biosample_sorted:
        biosample_id    = each_biosample_subject.split('_')[0]
        biosample_title = biosample_title_dict.get(biosample_id, 'NA')
        subject_desc_list_sorted.append(biosample_title)

    # write out to file
    op_txt_handle = open(op_txt, 'w')
    op_txt_handle.write('Genome\t%s\n' % ('\t'.join(subject_id_set_all_with_biosample_sorted)))
    for each_gnm in sorted(list(pebblescout_op_dod.keys())):
        current_gnm_dict = pebblescout_op_dod[each_gnm]
        value_list = [each_gnm]
        for each_biosample_subject in subject_id_set_all_with_biosample_sorted:
            subject_id = each_biosample_subject.split('_')[1]
            value_list.append(current_gnm_dict.get(subject_id, '0'))
        op_txt_handle.write('%s\n' % ('\t'.join(value_list)))
    op_txt_handle.write('Description\t%s\n' % ('\t'.join(subject_desc_list_sorted)))
    op_txt_handle.close()

    transpose_csv(op_txt, op_txt_t, '\t', 0, 0)
    os.system('rm -r %s' % op_txt)

    # metadata_txt_handle = open(metadata_txt, 'w')
    # for each_biosample_subject in subject_id_set_all_with_biosample_sorted:
    #     bioSample_id    = each_biosample_subject.split('_')[0]
    #     subject_id      = each_biosample_subject.split('_')[1]
    #     bioSample_title = biosample_title_dict[bioSample_id]
    #     metadata_txt_handle.write('%s\t%s\t%s\n' % (subject_id, bioSample_id, bioSample_title))
    # metadata_txt_handle.close()


def PebbleScout(args):

    file_dir            = args['i']
    file_ext            = args['x']
    db_str              = args['db']
    # min_cov_pct         = args['min_cov']
    num_threads         = args['t']
    force_create_op_dir = args['f']
    op_dir              = args['o']
    gnm_id_txt          = args['id']
    op_prefix           = args['p']

    if op_prefix is None:
        prefix_str = ''
    else:
        prefix_str = '%s_' % op_prefix

    cov_pct_list        = [1, 3, 5, 10, 20, 30, 50, 0]

    # define file name
    gnm_dir                 = '%s/genomes_concatenated'             % op_dir
    pebblescout_op_dir      = '%s/pebblescout_op'                   % op_dir
    cmd_txt                 = '%s/commands.txt'                     % op_dir
    metadata_txt            = '%s/metadata.txt'                     % op_dir

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % gnm_dir)
    os.system('mkdir %s' % pebblescout_op_dir)

    db_list = db_str.split(',')

    gnm_id_set = set()
    if gnm_id_txt is not None:
        if os.path.isfile(gnm_id_txt) is True:
            for each_gnm in open(gnm_id_txt):
                gnm_id_set.add(each_gnm.strip().split()[0])
        else:
            print('%s not found, program exited!' % gnm_id_txt)
            exit()

    # concatenate sequences in input genomes
    file_re = '%s/*.%s' % (file_dir, file_ext)
    file_list = glob.glob(file_re)
    file_base_list_to_process = []
    for each_gnm in file_list:
        f_name, f_path, f_base, f_ext = sep_path_basename_ext(each_gnm)

        to_process = ''
        if len(gnm_id_set) == 0:
            to_process = True
        else:
            if (f_base in gnm_id_set) or (f_name in gnm_id_set):
                to_process = True

        if to_process is True:
            file_base_list_to_process.append(f_base)
            current_gnm_after_cat = '%s/%s.%s' % (gnm_dir, f_base, f_ext)
            concatenates_seq_str = ''
            for each_seq in SeqIO.parse(each_gnm, 'fasta'):
                concatenates_seq_str = concatenates_seq_str + str(each_seq.seq).strip()
            current_gnm_after_cat_handle = open(current_gnm_after_cat, 'w')
            current_gnm_after_cat_handle.write('>%s\n' % f_base)
            current_gnm_after_cat_handle.write(concatenates_seq_str + '\n')
            current_gnm_after_cat_handle.close()

    # get pebblescout commands
    cmd_txt_handle = open(cmd_txt, 'w')
    cmd_set = set()
    for f_base in file_base_list_to_process:
        for each_db in db_list:
            pebblescout_cmd = 'curl -s -F "fasta=@%s/%s.%s" "https://pebblescout.ncbi.nlm.nih.gov/sra-cl-be/sra-cl-be.cgi?db=%s&m=2&rettype=pebblescout&download=yes" > %s/%s_%s.txt' % (gnm_dir, f_base, file_ext, each_db, pebblescout_op_dir, f_base, each_db)
            cmd_set.add(pebblescout_cmd)
            cmd_txt_handle.write(pebblescout_cmd + '\n')
    cmd_txt_handle.close()

    # run pebblescout with multi-processing
    print('Running %s commands with %s cores' % (len(file_base_list_to_process), num_threads))
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, sorted([i for i in cmd_set]))
    pool.close()
    pool.join()

    # get data matrix
    for each_cov_pct in cov_pct_list:
        df_txt       = '%s/%sPebbleScout_all_tmp.txt'   % (op_dir, prefix_str)
        df_txt_t     = '%s/%sPebbleScout_all.txt'       % (op_dir, prefix_str)
        if each_cov_pct > 0:
            df_txt   = '%s/%sPebbleScout_cov%s_tmp.txt' % (op_dir, prefix_str, each_cov_pct)
            df_txt_t = '%s/%sPebbleScout_cov%s.txt'     % (op_dir, prefix_str, each_cov_pct)
        pebblescout_results_to_data_matrix(file_base_list_to_process, db_list, pebblescout_op_dir, each_cov_pct, df_txt, df_txt_t, metadata_txt)

    print('Done!')


if __name__ == '__main__':

    PebbleScout_parser = argparse.ArgumentParser(usage=PebbleScout_usage)
    PebbleScout_parser.add_argument('-p',           required=False, default=None,                   help='output prefix')
    PebbleScout_parser.add_argument('-i',           required=True,                                  help='input fasta file/dir')
    PebbleScout_parser.add_argument('-x',           required=False, default=None,                   help='file extension')
    PebbleScout_parser.add_argument('-id',          required=False, default=None,                   help='id of genomes to process')
    PebbleScout_parser.add_argument('-db',          required=False, default='meta,meta_vol2',       help='query database, default is meta and meta_vol2')
    # PebbleScout_parser.add_argument('-min_cov',     required=False, type=str, default='1,3,5,10',   help='minimum %coverage to include in the datamatrix, default is 1,3,5,10')
    PebbleScout_parser.add_argument('-t',           required=False, type=int, default=1,            help='number of core, default is 1')
    PebbleScout_parser.add_argument('-f',           required=False, action="store_true",            help='force overwrite')
    PebbleScout_parser.add_argument('-o',           required=True,                                  help='output directory')
    args = vars(PebbleScout_parser.parse_args())
    PebbleScout(args)
