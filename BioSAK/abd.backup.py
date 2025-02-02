import os
import glob
import math
import argparse
import pandas as pd
from time import sleep
import multiprocessing as mp
from functools import reduce


abd_usage = '''
=========================== abd example commands ===========================

BioSAK abd -i sample.txt -r reference.fa -mask -t 36
BioSAK abd -i sample.txt -r reference_masked.fa -t 36

# input txt file format (tab separated)
Seawater1	path/to/Seawater1_1.fastq.gz	path/to/Seawater1_2.fastq.gz
Seawater2	path/to/Seawater2_1.fastq.gz	path/to/Seawater2_2.fastq.gz

============================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def get_abd1_mask(fasta_file, op_dir, num_threads):

    _, _, f_base, _ = sep_path_basename_ext(fasta_file)

    # get rRNA regions
    shell_cmd_1_metaxa2_ssu_cmd     = 'metaxa2 --plus --mode m --cpu %s --multi_thread T --table T -g ssu --not_found T -i %s -o %s/%s.metaxa2_ssu'             % (num_threads, fasta_file, op_dir, f_base)
    shell_cmd_2_metaxa2_lsu_cmd     = 'metaxa2 --plus --mode m --cpu %s --multi_thread T --table T -g lsu --not_found T -i %s -o %s/%s.metaxa2_lsu'             % (num_threads, fasta_file, op_dir, f_base)
    shell_cmd_3                     = 'cut -f 1,9,10 %s/%s.metaxa2_ssu.extraction.results > %s/%s.masked_metaxa_ssu.bed'                                        % (op_dir, f_base, op_dir, f_base)
    shell_cmd_4                     = 'cut -f 1,9,10 %s/%s.metaxa2_lsu.extraction.results > %s/%s.masked_metaxa_lsu.bed'                                        % (op_dir, f_base, op_dir, f_base)
    shell_cmd_5                     = 'cat %s/%s.masked_metaxa_ssu.bed %s/%s.masked_metaxa_lsu.bed > %s/%s.masked_metaxa.bed'                                   % (op_dir, f_base, op_dir, f_base, op_dir, f_base)
    shell_cmd_6_barrnap_bac_cmd     = 'barrnap --kingdom bac --threads %s --reject 0.3 %s > %s/%s.barrnap_bac.gff'                                              % (num_threads, fasta_file, op_dir, f_base)
    shell_cmd_7_barrnap_arc_cmd     = 'barrnap --kingdom arc --threads %s --reject 0.3 %s > %s/%s.barrnap_arc.gff'                                              % (num_threads, fasta_file, op_dir, f_base)
    shell_cmd_8_barrnap_euk_cmd     = 'barrnap --kingdom euk --threads %s --reject 0.3 %s > %s/%s.barrnap_euk.gff'                                              % (num_threads, fasta_file, op_dir, f_base)
    shell_cmd_9                     = 'cut -f 1,4,5 %s/%s.barrnap_bac.gff > %s/%s.masked_barrnap_bac.bed'                                                       % (op_dir, f_base, op_dir, f_base)
    shell_cmd_10                    = 'cut -f 1,4,5 %s/%s.barrnap_arc.gff > %s/%s.masked_barrnap_arc.bed'                                                       % (op_dir, f_base, op_dir, f_base)
    shell_cmd_11                    = 'cut -f 1,4,5 %s/%s.barrnap_euk.gff > %s/%s.masked_barrnap_euk.bed'                                                       % (op_dir, f_base, op_dir, f_base)
    shell_cmd_12                    = 'cat %s/%s.masked_barrnap_bac.bed %s/%s.masked_barrnap_arc.bed %s/%s.masked_barrnap_euk.bed > %s/%s.masked_barrnap.bed'   % (op_dir, f_base, op_dir, f_base, op_dir, f_base, op_dir, f_base)
    shell_cmd_13                    = 'cat %s/%s.masked_barrnap.bed %s/%s.masked_metaxa.bed > %s/%s.masked_rrna.bed'                                            % (op_dir, f_base, op_dir, f_base, op_dir, f_base)

    # get tRNA regions
    shell_cmd_14_trna_scan_arc_cmd  = 'tRNAscan-SE -A -b %s/%s.ar.bedformat --thread 36 %s'                                                                     % (op_dir, f_base, fasta_file)
    shell_cmd_15_trna_scan_bac_cmd  = 'tRNAscan-SE -B -b %s/%s.bac.bedformat --thread 36 %s'                                                                    % (op_dir, f_base, fasta_file)
    shell_cmd_16                    = 'cat %s/%s.ar.bedformat %s/%s.bac.bedformat > %s/%s.bedformat'                                                            % (op_dir, f_base, op_dir, f_base, op_dir, f_base)
    shell_cmd_17                    = 'cut -f 1,2,3 %s/%s.bedformat > %s/%s.masked_trnascan.bed'                                                                % (op_dir, f_base, op_dir, f_base)

    # find duplicated areas
    shell_cmd_18                    = 'dustmasker -in %s -out %s/%s.lowcom.out -outfmt acclist'                                                                 % (fasta_file, op_dir, f_base)
    shell_cmd_19                    = "cut -d '>' -f 2 %s/%s.lowcom.out > %s/%s.masked_dustmasker.bed"                                                          % (op_dir, f_base, op_dir, f_base)

    # cover the above area with NNNN
    shell_cmd_20                    = "cat %s/%s.masked_rrna.bed %s/%s.masked_trnascan.bed %s/%s.masked_dustmasker.bed > %s/%s.mask.bed"                        % (op_dir, f_base, op_dir, f_base, op_dir, f_base, op_dir, f_base)
    shell_cmd_21                    = "awk -F'\t' '$2>=0' %s/%s.mask.bed > %s/%s.mask_0.bed"                                                                    % (op_dir, f_base, op_dir, f_base)
    shell_cmd_22                    = "bedtools maskfasta -fi %s -bed %s/%s.mask_0.bed -fo %s/%s.masked.fasta"                                                  % (fasta_file, op_dir, f_base, op_dir, f_base)

    # run cmds
    print(shell_cmd_1_metaxa2_ssu_cmd)
    os.system(shell_cmd_1_metaxa2_ssu_cmd)

    print(shell_cmd_2_metaxa2_lsu_cmd)
    os.system(shell_cmd_2_metaxa2_lsu_cmd)

    print(shell_cmd_3)
    os.system(shell_cmd_3)

    print(shell_cmd_4)
    os.system(shell_cmd_4)

    print(shell_cmd_5)
    os.system(shell_cmd_5)

    print(shell_cmd_6_barrnap_bac_cmd)
    os.system(shell_cmd_6_barrnap_bac_cmd)

    print(shell_cmd_7_barrnap_arc_cmd)
    os.system(shell_cmd_7_barrnap_arc_cmd)

    print(shell_cmd_8_barrnap_euk_cmd)
    os.system(shell_cmd_8_barrnap_euk_cmd)

    print(shell_cmd_9)
    os.system(shell_cmd_9)

    print(shell_cmd_10)
    os.system(shell_cmd_10)

    print(shell_cmd_11)
    os.system(shell_cmd_11)

    print(shell_cmd_12)
    os.system(shell_cmd_12)

    print(shell_cmd_13)
    os.system(shell_cmd_13)

    print(shell_cmd_14_trna_scan_arc_cmd)
    os.system(shell_cmd_14_trna_scan_arc_cmd)

    print(shell_cmd_15_trna_scan_bac_cmd)
    os.system(shell_cmd_15_trna_scan_bac_cmd)

    print(shell_cmd_16)
    os.system(shell_cmd_16)

    print(shell_cmd_17)
    os.system(shell_cmd_17)

    print(shell_cmd_18)
    os.system(shell_cmd_18)

    print(shell_cmd_19)
    os.system(shell_cmd_19)

    print(shell_cmd_20)
    os.system(shell_cmd_20)

    print(shell_cmd_21)
    os.system(shell_cmd_21)

    print(shell_cmd_22)
    os.system(shell_cmd_22)

    print('Masking reference sequences finished!')


def get_abd2_mapping_worker(arg_list):

    if len(arg_list) == 6:
        op_prefix       = arg_list[0]
        fq_r1           = arg_list[1]
        fq_r2           = arg_list[2]
        ref_seq         = arg_list[3]
        op_dir          = arg_list[4]
        num_threads     = arg_list[5]
    elif len(arg_list) == 5:
        op_prefix       = arg_list[0]
        fq_nonpaired    = arg_list[1]
        ref_seq         = arg_list[2]
        op_dir          = arg_list[3]
        num_threads     = arg_list[4]

    # check input files
    if os.path.isfile(ref_seq) is False:
        print('%s not found, program exited!' % ref_seq)
        exit()

    if len(arg_list) == 6:
        if os.path.isfile(fq_r1) is False:
            print('%s not found, program exited!' % fq_r1)
            exit()
        if os.path.isfile(fq_r2) is False:
            print('%s not found, program exited!' % fq_r2)
            exit()

    if len(arg_list) == 5:
        if os.path.isfile(fq_nonpaired) is False:
            print('%s not found, program exited!' % fq_nonpaired)
            exit()

    os.system('mkdir %s' % op_dir)

    if len(arg_list) == 5:
        fq_nonpaired_name, _, _, _ = sep_path_basename_ext(fq_nonpaired)
    if len(arg_list) == 6:
        fq1_name, _, _, _ = sep_path_basename_ext(fq_r1)
        fq2_name, _, _, _ = sep_path_basename_ext(fq_r2)
    ref_name, _, _, _ = sep_path_basename_ext(ref_seq)

    perform_decompress          = False
    if len(arg_list) == 5:
        fq_nonpaired_decompressed = fq_nonpaired
        if fq_nonpaired.endswith('.gz'):
            fq_nonpaired_decompressed   = '%s/%s'               % (op_dir, fq_nonpaired_name[:-3])
            gunzip_cmd_nonpaired        = 'gunzip -c %s > %s'   % (fq_nonpaired, fq_nonpaired_decompressed)
            perform_decompress          = True
            os.system(gunzip_cmd_nonpaired)

    if len(arg_list) == 6:
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
    index_ref_cmd     = 'cp %s %s/; bwa index %s/%s'                                                                                                                        % (ref_seq, op_dir, op_dir, ref_name)

    if len(arg_list) == 6:
        bwa_cmd           = 'bwa mem -5SP -t %s %s/%s %s %s | samblaster > %s/%s.sam'                                                                                       % (num_threads, op_dir, ref_name, fq_r1_decompressed, fq_r2_decompressed, op_dir, op_prefix)
    if len(arg_list) == 5:
        bwa_cmd           = 'bwa mem -5SP -t %s %s/%s %s | samblaster > %s/%s.sam'                                                                                          % (num_threads, op_dir, ref_name, fq_nonpaired_decompressed, op_dir, op_prefix)

    samtools_view_cmd = 'samtools view -@ 32 -bS -h -b %s/%s.sam > %s/%s.bam'                                                                                               % (op_dir, op_prefix, op_dir, op_prefix)
    samtools_sort_cmd = 'samtools sort -@ 32 %s/%s.bam -o %s/%s.sorted.bam'                                                                                                 % (op_dir, op_prefix, op_dir, op_prefix)
    coverm_filter_cmd = 'coverm filter -b %s/%s.sorted.bam --min-read-aligned-percent 0.9 --min-read-percent-identity 0.99 --output-bam-files %s/%s.sorted_filtered.bam'    % (op_dir, op_prefix, op_dir, op_prefix)
    pileup_sh_cmd     = 'pileup.sh in=%s/%s.sorted_filtered.bam out=%s/%s.sorted_filtered.cov rpkm=%s/%s.sorted_filtered.rpkm overwrite=true'                               % (op_dir, op_prefix, op_dir, op_prefix, op_dir, op_prefix)
    if len(arg_list) == 6:
        seqkit_stat_cmd   = 'seqkit stat %s > %s/%s.stat'                                                                                                                       % (fq_r1_decompressed, op_dir, op_prefix)
    if len(arg_list) == 5:
        seqkit_stat_cmd   = 'seqkit stat %s > %s/%s.stat'                                                                                                                       % (fq_nonpaired_decompressed, op_dir, op_prefix)

    print(seqkit_stat_cmd)
    os.system(seqkit_stat_cmd)

    print(index_ref_cmd)
    os.system(index_ref_cmd)

    print(bwa_cmd)
    os.system(bwa_cmd)
    if perform_decompress is True:
        if len(arg_list) == 6:
            os.system('rm %s/%s' % (op_dir, fq1_name[:-3]))
            os.system('rm %s/%s' % (op_dir, fq2_name[:-3]))
        if len(arg_list) == 5:
            os.system('rm %s/%s' % (op_dir, fq_nonpaired_name[:-3]))

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
    sleep(1)
    os.system('rm %s/%s.sorted_filtered.bam' % (op_dir, op_prefix))


def get_abd3_stats(rpkm_stat_cov_dir, stats_dir):

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
    gnm_level_rpkm_dir       = '%s/rpkm_by_genome'    % stats_dir
    combined_rpkm_file       = '%s/df_rpkm.txt'       % stats_dir
    combined_rpkm_file_log10 = '%s/df_rpkm.log10.txt' % stats_dir

    # create op dir
    os.system('mkdir %s' % stats_dir)
    os.system('mkdir %s' % gnm_level_rpkm_dir)

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


def abd(args):

    input_txt           = args['i']
    ref_seq             = args['r']
    op_dir              = args['o']
    run_mask            = args['mask']
    non_paired_reads    = args['nonpaired']
    num_threads         = args['t']
    force_overwrite     = args['f']

    # define file name
    mask_wd             = '%s/mask_wd'              % op_dir
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

    #################### mask reference sequences ####################

    ref_seq_masked = ref_seq
    if run_mask is True:
        os.system('mkdir %s' % mask_wd)
        get_abd1_mask(ref_seq, mask_wd, num_threads)
        _, _, ref_seq_base, _ = sep_path_basename_ext(ref_seq)
        ref_seq_masked = '%s/%s.masked.fasta' % (mask_wd, ref_seq_base)

    #################### mapping ####################

    # run mapping
    sample_lol = []
    for each_sample in open(input_txt):
        current_arg_list = each_sample.strip().split('\t')
        sample_lol.append(current_arg_list)

    threads_per_job = num_threads//len(sample_lol)

    arg_lol= []
    for arg_list in sample_lol:
        sample_id = arg_list[0]
        current_op_dir = '%s/%s_mapping_wd' % (op_dir, sample_id)
        arg_list.append(ref_seq_masked)
        arg_list.append(current_op_dir)
        arg_list.append(threads_per_job)
        arg_lol.append(arg_list)

    # run mapping with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(get_abd2_mapping_worker, arg_lol)
    pool.close()
    pool.join()

    #################### get stats ####################

    # mkdir
    os.mkdir(cov_rpkm_stat_dir)
    cp_cmd = 'cp %s/*mapping_wd/* %s' % (op_dir, cov_rpkm_stat_dir)
    os.system(cp_cmd)

    # get stats
    get_abd3_stats(cov_rpkm_stat_dir, stats_dir)

    #################### final report ####################

    # final report
    print('Estimated abundance exported to:')
    print('%s/df_rpkm.txt'       % stats_dir)
    print('%s/df_rpkm.log10.txt' % stats_dir)

    print('Done!')


if __name__ == '__main__':

    get_abd_parser = argparse.ArgumentParser()
    get_abd_parser.add_argument('-i',           required=True,                          help='input txt')
    get_abd_parser.add_argument('-r',           required=True,                          help='reference, need to be masked if you do not specify -mask')
    get_abd_parser.add_argument('-o',           required=True,                          help='output directory')
    get_abd_parser.add_argument('-mask',        required=False, action="store_true",    help='specify to mask reference sequence')
    get_abd_parser.add_argument('-nonpaired',   required=False, action="store_true",    help='specify for nonpaired reads, including long reads')
    get_abd_parser.add_argument('-t',           required=False, type=int, default=1,    help='number of threads, default is 1')
    get_abd_parser.add_argument('-f',           required=False, action="store_true",    help='force overwrite')
    args = vars(get_abd_parser.parse_args())
    abd(args)
