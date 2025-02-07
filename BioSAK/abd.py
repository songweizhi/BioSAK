import os
import glob
import math
import random
import argparse
import pandas as pd
from Bio import SeqIO
from time import sleep
from Bio.Seq import Seq
import multiprocessing as mp
from functools import reduce


abd_usage = '''
=========================== abd example commands ===========================

BioSAK abd -i sample.txt -o op_dir -t 36 -r reference.fa -mask
BioSAK abd -i sample.txt -o op_dir -t 36 -r reference_masked.fa

# sample.txt file format (tab separated)
Seawater1   path/to/Seawater1_1.fastq.gz    path/to/Seawater1_2.fastq.gz
Seawater2	path/to/Seawater2_1.fastq.gz    path/to/Seawater2_2.fastq.gz
Aphrocallistes  AphrocallistesBeatrix_subset.fastq

# Note
1. Input reads need to be in pair.
2. No underscore in genome name.
3. Sequence id format: gnm1_1, gnm1_2, gnm1_3 ...

============================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def Reads_simulator(op_dir, op_prefix, pwd_genome_file, read_number, read_depth, read_length, insert_size):

    # read in reference sequences
    ref_seq_total_len = 0
    ref_seq_dict = {}
    ref_seq_len_dict = {}
    for each_seq in SeqIO.parse(pwd_genome_file, 'fasta'):
        seq_id = each_seq.id
        seq_str = str(each_seq.seq)
        ref_seq_dict[seq_id] = seq_str
        ref_seq_len_dict[seq_id] = len(seq_str)
        ref_seq_total_len += len(seq_str)

    # calculate the number of reads to simulate
    if (read_number is None) and (read_depth is not None):
        read_number = round(ref_seq_total_len*read_depth/(read_length*2))
        print('The number of read pairs to simulate: %s' % read_number)

    # get the number of reads to simulate from each ref seq
    ref_to_read_num_dict = {}
    current_ref_index = 1
    total_assigned_read_num = 0
    for each_ref_seq in ref_seq_len_dict:
        if current_ref_index == len(ref_seq_len_dict):
            current_ref_seq_read_num = read_number - total_assigned_read_num
        else:
            current_ref_seq_len_pct = ref_seq_len_dict[each_ref_seq]/ref_seq_total_len
            current_ref_seq_read_num = round(current_ref_seq_len_pct*read_number)
            total_assigned_read_num += current_ref_seq_read_num
        ref_to_read_num_dict[each_ref_seq] = current_ref_seq_read_num
        current_ref_index += 1

    # simulate reads
    output_r1 = '%s/%s_R1.fasta' % (op_dir, op_prefix)
    output_r2 = '%s/%s_R2.fasta' % (op_dir, op_prefix)
    output_r1_handle = open(output_r1, 'w')
    output_r2_handle = open(output_r2, 'w')

    overall_n = 1
    for each_ref_seq in ref_to_read_num_dict:
        ref_seq_str = ref_seq_dict[each_ref_seq]
        ref_seq_read_num = ref_to_read_num_dict[each_ref_seq]
        sequence_length = len(ref_seq_str)
        fragment_length = 2 * read_length + insert_size

        n = 1
        while n <= ref_seq_read_num:
            rdm_num = random.randint(1, sequence_length)
            current_fragment = ''
            if (rdm_num + fragment_length) <= sequence_length:
                current_fragment = ref_seq_str[rdm_num - 1: rdm_num + fragment_length - 1]
            elif (rdm_num + fragment_length) > sequence_length:
                seq_part_1_seq = ref_seq_str[rdm_num - 1:]
                seq_part_2_seq = ref_seq_str[:fragment_length - sequence_length + rdm_num - 1]
                current_fragment = seq_part_1_seq + seq_part_2_seq
            current_fragment_r1 = current_fragment[:read_length]
            current_fragment_r2 = current_fragment[-read_length:]
            current_fragment_r2_reverse_complement = str(Seq(current_fragment_r2).reverse_complement())
            current_read_r1_id = '%s_%s.1' % (op_prefix, overall_n)
            current_read_r2_id = '%s_%s.2' % (op_prefix, overall_n)
            output_r1_handle.write('>%s\n' % current_read_r1_id)
            output_r1_handle.write('%s\n'  % current_fragment_r1)
            output_r2_handle.write('>%s\n' % current_read_r2_id)
            output_r2_handle.write('%s\n'  % current_fragment_r2_reverse_complement)
            n += 1
            overall_n += 1
    output_r1_handle.close()
    output_r2_handle.close()
    print('Simulated reads exported to %s and %s' % (output_r1, output_r2))


def fq2fa(fq_in, fa_out):
    fa_out_handle = open(fa_out, 'w')
    for each_long_read in SeqIO.parse(fq_in, 'fastq'):
        fa_out_handle.write('>%s\n%s\n' % (each_long_read.id, each_long_read.seq))
    fa_out_handle.close()


# paired reads will be simulated from the unpaired long reads
def get_abd2_mapping_worker_long(arg_list):

    op_prefix               = arg_list[0]
    fq_long                 = arg_list[1]
    ref_seq                 = arg_list[2]
    op_dir                  = arg_list[3]
    num_threads             = arg_list[4]
    keep_bam_file           = arg_list[5]
    simulate_num            = arg_list[6]
    simulate_depth          = arg_list[7]
    simulate_len            = arg_list[8]
    simulate_insert_size    = arg_list[9]
    ref_index_dir           = arg_list[10]

    # check input files
    if os.path.isfile(ref_seq) is False:
        print('%s not found, program exited!' % ref_seq)
        exit()

    if os.path.isfile(fq_long) is False:
        print('%s not found, program exited!' % fq_long)
        exit()

    os.system('mkdir %s' % op_dir)

    fq_nonpaired_name, _, _, _ = sep_path_basename_ext(fq_long)
    ref_name, _, _, _ = sep_path_basename_ext(ref_seq)

    perform_decompress              = False
    fq_nonpaired_decompressed       = fq_long
    if fq_long.endswith('.gz'):
        fq_nonpaired_decompressed   = '%s/%s'               % (op_dir, fq_nonpaired_name[:-3])
        gunzip_cmd_long             = 'gunzip -c %s > %s'   % (fq_long, fq_nonpaired_decompressed)
        perform_decompress          = True
        os.system(gunzip_cmd_long)

    # simulate reads here
    long_reads_fasta = '%s/%s_long.fasta' % (op_dir, op_prefix)
    fq2fa(fq_nonpaired_decompressed, long_reads_fasta)
    Reads_simulator(op_dir, op_prefix, long_reads_fasta, simulate_num, simulate_depth, simulate_len, simulate_insert_size)

    r1_fa = '%s/%s_R1.fasta' % (op_dir, op_prefix)
    r2_fa = '%s/%s_R1.fasta' % (op_dir, op_prefix)

    # index reference sequences
    bwa_cmd           = 'bwa mem -5SP -t %s %s/%s %s %s | samblaster > %s/%s.sam'                                                                                           % (num_threads, ref_index_dir, ref_name, r1_fa, r2_fa, op_dir, op_prefix)
    samtools_view_cmd = 'samtools view -@ 32 -bS -h -b %s/%s.sam > %s/%s.bam'                                                                                               % (op_dir, op_prefix, op_dir, op_prefix)
    samtools_sort_cmd = 'samtools sort -@ 32 %s/%s.bam -o %s/%s.sorted.bam'                                                                                                 % (op_dir, op_prefix, op_dir, op_prefix)
    coverm_filter_cmd = 'coverm filter -b %s/%s.sorted.bam --min-read-aligned-percent 0.9 --min-read-percent-identity 0.99 --output-bam-files %s/%s.sorted_filtered.bam'    % (op_dir, op_prefix, op_dir, op_prefix)
    pileup_sh_cmd     = 'pileup.sh in=%s/%s.sorted_filtered.bam out=%s/%s.sorted_filtered.cov rpkm=%s/%s.sorted_filtered.rpkm overwrite=true'                               % (op_dir, op_prefix, op_dir, op_prefix, op_dir, op_prefix)
    seqkit_stat_cmd   = 'seqkit stat %s > %s/%s.stat'                                                                                                                       % (r1_fa, op_dir, op_prefix)

    print(seqkit_stat_cmd)
    os.system(seqkit_stat_cmd)

    print(bwa_cmd)
    os.system(bwa_cmd)
    if perform_decompress is True:
        os.system('rm %s/%s' % (op_dir, fq_nonpaired_name[:-3]))

    os.system('rm %s'       % long_reads_fasta)
    os.system('rm %s'       % r1_fa)
    os.system('rm %s'       % r2_fa)
    os.system('rm %s/%s'    % (op_dir, ref_name))
    os.system('rm %s/%s.*'  % (op_dir, ref_name))

    print(samtools_view_cmd)
    os.system(samtools_view_cmd)
    sleep(1)
    os.system('rm %s/%s.sam' % (op_dir, op_prefix))

    print(samtools_sort_cmd)
    os.system(samtools_sort_cmd)
    sleep(1)
    os.system('rm %s/%s.bam' % (op_dir, op_prefix))

    print(coverm_filter_cmd)
    os.system(coverm_filter_cmd)
    sleep(1)
    os.system('rm %s/%s.sorted.bam' % (op_dir, op_prefix))

    print(pileup_sh_cmd)
    os.system(pileup_sh_cmd)
    if keep_bam_file is False:
        os.system('rm %s/%s.sorted_filtered.bam' % (op_dir, op_prefix))
    else:
        os.system('samtools index %s/%s.sorted_filtered.bam' % (op_dir, op_prefix))


def get_abd2_mapping_worker_paired(arg_list):

    op_prefix       = arg_list[0]
    fq_r1           = arg_list[1]
    fq_r2           = arg_list[2]
    ref_seq         = arg_list[3]
    op_dir          = arg_list[4]
    num_threads     = arg_list[5]
    keep_bam_file   = arg_list[6]
    ref_index_dir   = arg_list[7]

    # check input files
    if os.path.isfile(ref_seq) is False:
        print('%s not found, program exited!' % ref_seq)
        exit()

    if os.path.isfile(fq_r1) is False:
        print('%s not found, program exited!' % fq_r1)
        exit()
    if os.path.isfile(fq_r2) is False:
        print('%s not found, program exited!' % fq_r2)
        exit()

    os.system('mkdir %s' % op_dir)

    fq1_name, _, _, _ = sep_path_basename_ext(fq_r1)
    fq2_name, _, _, _ = sep_path_basename_ext(fq_r2)
    ref_name, _, _, _ = sep_path_basename_ext(ref_seq)

    perform_decompress = False

    fq_r1_decompressed = fq_r1
    fq_r2_decompressed = fq_r2
    if fq_r1.endswith('.gz'):
        fq_r1_decompressed  = '%s/%s'               % (op_dir, fq1_name[:-3])
        fq_r2_decompressed  = '%s/%s'               % (op_dir, fq2_name[:-3])
        gunzip_cmd_r1       = 'gunzip -c %s > %s'   % (fq_r1, fq_r1_decompressed)
        gunzip_cmd_r2       = 'gunzip -c %s > %s'   % (fq_r2, fq_r2_decompressed)
        perform_decompress  = True
        os.system(gunzip_cmd_r1)
        os.system(gunzip_cmd_r2)

    # index reference sequences
    bwa_cmd           = 'bwa mem -5SP -t %s %s/%s %s %s | samblaster > %s/%s.sam'                                                                                           % (num_threads, ref_index_dir, ref_name, fq_r1_decompressed, fq_r2_decompressed, op_dir, op_prefix)
    samtools_view_cmd = 'samtools view -@ 32 -bS -h -b %s/%s.sam > %s/%s.bam'                                                                                               % (op_dir, op_prefix, op_dir, op_prefix)
    samtools_sort_cmd = 'samtools sort -@ 32 %s/%s.bam -o %s/%s.sorted.bam'                                                                                                 % (op_dir, op_prefix, op_dir, op_prefix)
    coverm_filter_cmd = 'coverm filter -b %s/%s.sorted.bam --min-read-aligned-percent 0.9 --min-read-percent-identity 0.99 --output-bam-files %s/%s.sorted_filtered.bam'    % (op_dir, op_prefix, op_dir, op_prefix)
    pileup_sh_cmd     = 'pileup.sh in=%s/%s.sorted_filtered.bam out=%s/%s.sorted_filtered.cov rpkm=%s/%s.sorted_filtered.rpkm overwrite=true'                               % (op_dir, op_prefix, op_dir, op_prefix, op_dir, op_prefix)
    seqkit_stat_cmd   = 'seqkit stat %s > %s/%s.stat'                                                                                                                       % (fq_r1_decompressed, op_dir, op_prefix)

    print(seqkit_stat_cmd)
    os.system(seqkit_stat_cmd)

    print(bwa_cmd)
    os.system(bwa_cmd)
    if perform_decompress is True:
        os.system('rm %s/%s' % (op_dir, fq1_name[:-3]))
        os.system('rm %s/%s' % (op_dir, fq2_name[:-3]))
    os.system('rm %s/%s'   % (op_dir, ref_name))
    os.system('rm %s/%s.*' % (op_dir, ref_name))

    print(samtools_view_cmd)
    os.system(samtools_view_cmd)
    sleep(1)
    os.system('rm %s/%s.sam' % (op_dir, op_prefix))

    print(samtools_sort_cmd)
    os.system(samtools_sort_cmd)
    sleep(1)
    os.system('rm %s/%s.bam' % (op_dir, op_prefix))

    print(coverm_filter_cmd)
    os.system(coverm_filter_cmd)
    sleep(1)
    os.system('rm %s/%s.sorted.bam' % (op_dir, op_prefix))

    print(pileup_sh_cmd)
    os.system(pileup_sh_cmd)
    if keep_bam_file is False:
        os.system('rm %s/%s.sorted_filtered.bam' % (op_dir, op_prefix))
    else:
        os.system('samtools index %s/%s.sorted_filtered.bam' % (op_dir, op_prefix))


def get_abd3_stats(rpkm_stat_cov_dir, genome_rename_file, stats_dir):

    # check input files
    if os.path.isdir(rpkm_stat_cov_dir) is False:
        print('%s not found, program exited!' % rpkm_stat_cov_dir)
        exit()

    rpkm_file_re    = '%s/*.rpkm' % rpkm_stat_cov_dir
    stat_file_re    = '%s/*.stat' % rpkm_stat_cov_dir
    rpkm_file_list  = glob.glob(rpkm_file_re)
    stat_file_list  = glob.glob(stat_file_re)

    if len(rpkm_file_list) == 0:
        print('.rpkm file not found, program exited!')
        exit()

    if len(stat_file_list) == 0:
        print('.stat file not found, program exited!')
        exit()

    # define output file name
    gnm_level_rpkm_dir               = '%s/rpkm_by_genome'              % stats_dir
    combined_rpkm_file               = '%s/rpkm_df.txt'                 % stats_dir
    combined_rpkm_file_renamed       = '%s/rpkm_df.renamed.txt'         % stats_dir
    combined_rpkm_file_log10         = '%s/rpkm_df.log10.txt'           % stats_dir
    combined_rpkm_file_log10_renamed = '%s/rpkm_df.renamed.log10.txt'   % stats_dir

    # create op dir
    os.system('mkdir %s' % stats_dir)
    os.system('mkdir %s' % gnm_level_rpkm_dir)

    # get genome rename dict
    gnm_rename_dict = dict()
    if genome_rename_file is not None:
        for line in open(genome_rename_file):
            line_split = line.strip().split('\t')
            gnm_rename_dict[line_split[0]] = line_split[1]

    # read in stat file
    sample_total_read_num_dict = dict()
    for stat_file in sorted(stat_file_list):
        _, _, stat_base, _ = sep_path_basename_ext(stat_file)
        stat_file_lines = open(stat_file).readlines()
        info_line_split = stat_file_lines[-1].strip().split(' ')
        while '' in info_line_split:
            info_line_split.remove('')
        read_num = int(info_line_split[3].replace(',', ''))
        sample_total_read_num_dict[stat_base] = read_num

    # summarise rpkm on genome level
    sample_id_list = []
    for rpkm_file in sorted(rpkm_file_list):
        rpkm_name, rpkm_path, rpkm_base, rpkm_ext = sep_path_basename_ext(rpkm_file)
        sample_id = rpkm_base.split('.sorted_filtered')[0]
        sample_id_list.append(sample_id)
        op_txt = '%s/%s.cal_bin.rpkm' % (gnm_level_rpkm_dir, sample_id)
        rpkm_df = pd.read_csv(rpkm_file, sep='\t', header=4)
        subset_df_col_list = ['#Name', 'Length', 'Reads']
        sample_total_read_num = sample_total_read_num_dict[sample_id]
        subset_df = rpkm_df.loc[:, subset_df_col_list]
        subset_df['bin_name'] = subset_df['#Name'].str.rsplit('_', n=1, expand=True)[0]
        subset_df['Length'] = subset_df['Length'].astype(int)
        subset_df['Reads'] = subset_df['Reads'].astype(int)

        # calulate RPKM
        op_df = subset_df.groupby(["bin_name"])[['Length', 'Reads']].sum()
        op_df[sample_id] = (op_df['Reads']*1000000000)/(op_df['Length']*sample_total_read_num*2)
        op_df.to_csv(op_txt, sep='\t')

    # combine the results
    rpkm_df_list = []
    for sample_id in sample_id_list:
        gnm_level_rpkm_file = '%s/%s.cal_bin.rpkm' % (gnm_level_rpkm_dir, sample_id)
        gnm_level_rpkm_df = pd.read_csv(gnm_level_rpkm_file, sep='\t')
        gnm_level_rpkm_df = gnm_level_rpkm_df.drop(["Length", "Reads"], axis=1)
        rpkm_df_list.append(gnm_level_rpkm_df)

    # write out combined_rpkm_file
    rpkm_df_combined = reduce(lambda left, right: pd.merge(left, right, on='bin_name',how='outer'), rpkm_df_list)
    rpkm_df_combined.to_csv(combined_rpkm_file, sep='\t',index=False)

    # write out combined_rpkm_file_log10
    combined_rpkm_file_log10_handle = open(combined_rpkm_file_log10, 'w')
    for each_line in open(combined_rpkm_file):
        each_line_split = each_line.strip().split('\t')
        if each_line.startswith('bin_name'):
            combined_rpkm_file_log10_handle.write(each_line)
        else:
            value_list = [each_line_split[0]]
            for each_value in each_line_split[1:]:
                each_value = float(each_value)
                if each_value == 0:
                    value_list.append('na')
                else:
                    each_value_log = math.log(each_value)
                    value_list.append(each_value_log)
            combined_rpkm_file_log10_handle.write('%s\n' % '\t'.join([str(i) for i in value_list]))
    combined_rpkm_file_log10_handle.close()

    # rename gnm id in the final dataframe
    if len(gnm_rename_dict) > 0:

        # rename combined_rpkm_file_renamed
        combined_rpkm_file_renamed_handle = open(combined_rpkm_file_renamed, 'w')
        for each_line in open(combined_rpkm_file):
            if each_line.startswith('bin_name'):
                combined_rpkm_file_renamed_handle.write(each_line)
            else:
                each_line_split = each_line.strip().split('\t')
                gnm_id = each_line_split[0]
                gnm_id_new = gnm_rename_dict.get(gnm_id, gnm_id)
                combined_rpkm_file_renamed_handle.write(each_line.replace(gnm_id, gnm_id_new))
        combined_rpkm_file_renamed_handle.close()

        # rename combined_rpkm_file_log10_renamed
        combined_rpkm_file_log10_renamed_handle = open(combined_rpkm_file_log10_renamed, 'w')
        for each_line in open(combined_rpkm_file_log10):
            if each_line.startswith('bin_name'):
                combined_rpkm_file_log10_renamed_handle.write(each_line)
            else:
                each_line_split = each_line.strip().split('\t')
                gnm_id = each_line_split[0]
                gnm_id_new = gnm_rename_dict.get(gnm_id, gnm_id)
                combined_rpkm_file_log10_renamed_handle.write(each_line.replace(gnm_id, gnm_id_new))
        combined_rpkm_file_log10_renamed_handle.close()


def abd(args):

    input_txt               = args['i']
    ref_seq                 = args['r']
    op_dir                  = args['o']
    keep_bam_file           = args['keep_bam']
    num_threads             = args['t']
    force_overwrite         = args['f']
    rename_txt              = args['rename']

    # parameters for read_simulator
    read_simulator_num          = 5000000
    read_simulator_depth        = None
    read_simulator_len          = 150
    read_simulator_insert_size  = 200

    # check input file
    if rename_txt is not None:
        if os.path.isfile(rename_txt) is False:
            print('%s not found, program exited!' % rename_txt)
            exit()

    # define file name
    ref_seq_index_dir   = '%s/ref_seq_index'        % op_dir
    mapping_wd          = '%s/mapping_wd'           % op_dir
    cov_rpkm_stat_dir   = '%s/cov_rpkm_stat_files'  % op_dir
    stats_dir           = '%s/stats_dir'            % op_dir

    # create op dir
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output directory detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    #################### index reference sequences ####################

    os.system('mkdir %s' % ref_seq_index_dir)

    ref_name, _, _, _ = sep_path_basename_ext(ref_seq)
    index_ref_cmd     = 'cp %s %s/; bwa index %s/%s'    % (ref_seq, ref_seq_index_dir, ref_seq_index_dir, ref_name)
    print(index_ref_cmd)
    os.system(index_ref_cmd)

    #################### mapping ####################

    os.system('mkdir %s' % mapping_wd)

    # run mapping
    sample_lol_long = []
    sample_lol_paired = []
    for each_sample in open(input_txt):
        current_arg_list = each_sample.strip().split('\t')
        if len(current_arg_list) == 3:
            sample_lol_paired.append(current_arg_list)
        elif len(current_arg_list) == 2:
            sample_lol_long.append(current_arg_list)

    #################### mapping paired reads ####################

    if len(sample_lol_paired) > 0:
        threads_per_job_paired = num_threads//len(sample_lol_paired)

        arg_lol= []
        for arg_list in sample_lol_paired:
            sample_id = arg_list[0]
            current_op_dir = '%s/%s' % (mapping_wd, sample_id)
            arg_list.append(ref_seq)
            arg_list.append(current_op_dir)
            arg_list.append(threads_per_job_paired)
            arg_list.append(keep_bam_file)
            arg_list.append(ref_seq_index_dir)
            arg_list.append(mapping_wd)
            arg_lol.append(arg_list)

        # run mapping with multiprocessing
        pool = mp.Pool(processes=num_threads)
        pool.map(get_abd2_mapping_worker_paired, arg_lol)
        pool.close()
        pool.join()

    #################### mapping long reads ####################

    if len(sample_lol_long) > 0:
        threads_per_job_long = num_threads//len(sample_lol_long)

        arg_lol= []
        for arg_list in sample_lol_long:
            sample_id = arg_list[0]
            current_op_dir = '%s/%s' % (mapping_wd, sample_id)
            arg_list.append(ref_seq)
            arg_list.append(current_op_dir)
            arg_list.append(threads_per_job_long)
            arg_list.append(keep_bam_file)
            arg_list.append(read_simulator_num)
            arg_list.append(read_simulator_depth)
            arg_list.append(read_simulator_len)
            arg_list.append(read_simulator_insert_size)
            arg_list.append(ref_seq_index_dir)
            arg_list.append(mapping_wd)
            arg_lol.append(arg_list)

        # run mapping with multiprocessing
        pool = mp.Pool(processes=num_threads)
        pool.map(get_abd2_mapping_worker_long, arg_lol)
        pool.close()
        pool.join()

    #################### get stats ####################

    # mkdir
    os.mkdir(cov_rpkm_stat_dir)
    cp_cmd_rpkm = 'cp %s/*/*.rpkm %s'   % (mapping_wd, cov_rpkm_stat_dir)
    cp_cmd_stat = 'cp %s/*/*.stat %s'   % (mapping_wd, cov_rpkm_stat_dir)
    cp_cmd_cov  = 'cp %s/*/*.cov %s'    % (mapping_wd, cov_rpkm_stat_dir)
    os.system(cp_cmd_rpkm)
    os.system(cp_cmd_stat)
    os.system(cp_cmd_cov)

    # get stats
    get_abd3_stats(cov_rpkm_stat_dir, rename_txt, stats_dir)

    #################### final report ####################

    # final report
    print('Estimated abundance exported to:')
    print('%s/rpkm_df.txt'       % stats_dir)
    print('%s/rpkm_df.log10.txt' % stats_dir)
    print('Done!')


if __name__ == '__main__':

    abd_parser = argparse.ArgumentParser()
    abd_parser.add_argument('-i',           required=True,                          help='input txt')
    abd_parser.add_argument('-r',           required=True,                          help='reference, need to be masked')
    abd_parser.add_argument('-o',           required=True,                          help='output directory')
    abd_parser.add_argument('-rename',      required=False, default=None,           help='rename genome id in the final dataframe')
    abd_parser.add_argument('-keep_bam',    required=False, action="store_true",    help='do not delete bam file')
    abd_parser.add_argument('-t',           required=False, type=int, default=1,    help='number of threads, default is 1')
    abd_parser.add_argument('-f',           required=False, action="store_true",    help='force overwrite')
    args = vars(abd_parser.parse_args())
    abd(args)
