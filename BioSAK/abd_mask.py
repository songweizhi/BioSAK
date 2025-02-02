import os
import glob
import argparse
from Bio import SeqIO
import multiprocessing as mp


abd_mask_usage = '''
====================== abd_mask example commands ======================

BioSAK abd_mask -i gnm_dir -x fna -o mask_wd -p dRep99_406 -t 36 -f

=======================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def rename_and_cat_gnm(gnm_dir, gnm_ext, combined_renamed_fna, rename_txt):

    gnm_file_re = '%s/*.%s' % (gnm_dir, gnm_ext)
    gnm_file_list = glob.glob(gnm_file_re)

    gnm_rename_txt_handle = open(rename_txt, 'w')
    combined_renamed_fa_handle = open(combined_renamed_fna, 'w')
    for each_gnm in sorted(gnm_file_list):
        _, _, gnm_base, _ = sep_path_basename_ext(each_gnm)
        gnm_id_no_underscore = gnm_base.replace('_', '')
        gnm_rename_txt_handle.write('%s\t%s\n' % (gnm_base, gnm_id_no_underscore))

        seq_index = 1
        for each_seq in SeqIO.parse(each_gnm, 'fasta'):
            combined_renamed_fa_handle.write('>%s_%s\n%s\n' % (gnm_id_no_underscore, seq_index, str(each_seq.seq)))
            seq_index += 1

    gnm_rename_txt_handle.close()
    combined_renamed_fa_handle.close()


def abd_mask(args):

    seq_dir         = args['i']
    seq_ext         = args['x']
    op_dir          = args['o']
    op_prefix       = args['p']
    num_threads     = args['t']
    force_overwrite = args['f']

    ########################################## define file and directory name ##########################################

    renamed_combined_fna    = '%s/%s.fna'       % (op_dir, op_prefix)
    renamed_combined_txt    = '%s/%s.txt'       % (op_dir, op_prefix)
    tRNAscan_wd             = '%s/tRNAscan_wd'  % op_dir

    ####################################################################################################################

    seq_file_re     = '%s/*.%s' % (seq_dir, seq_ext)
    seq_file_list   = glob.glob(seq_file_re)
    if len(seq_file_list) == 0:
        print('No sequence files found in %s, program exited!' % seq_dir)
        exit()

    # create op dir
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output directory detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    # rename and combine sequences
    print('rename and combine input genomes')
    rename_and_cat_gnm(seq_dir, seq_ext, renamed_combined_fna, renamed_combined_txt)

    ####################################################################################################################

    # get rRNA regions
    shell_cmd_1_metaxa2_ssu_cmd     = 'metaxa2 --plus --mode m --cpu %s --multi_thread T --table T -g ssu --not_found T -i %s -o %s/%s.metaxa2_ssu'             % (num_threads, renamed_combined_fna, op_dir, op_prefix)
    shell_cmd_2_metaxa2_lsu_cmd     = 'metaxa2 --plus --mode m --cpu %s --multi_thread T --table T -g lsu --not_found T -i %s -o %s/%s.metaxa2_lsu'             % (num_threads, renamed_combined_fna, op_dir, op_prefix)
    shell_cmd_3                     = 'cut -f 1,9,10 %s/%s.metaxa2_ssu.extraction.results > %s/%s.masked_metaxa_ssu.bed'                                        % (op_dir, op_prefix, op_dir, op_prefix)
    shell_cmd_4                     = 'cut -f 1,9,10 %s/%s.metaxa2_lsu.extraction.results > %s/%s.masked_metaxa_lsu.bed'                                        % (op_dir, op_prefix, op_dir, op_prefix)
    shell_cmd_5                     = 'cat %s/%s.masked_metaxa_ssu.bed %s/%s.masked_metaxa_lsu.bed > %s/%s.masked_metaxa.bed'                                   % (op_dir, op_prefix, op_dir, op_prefix, op_dir, op_prefix)
    shell_cmd_6_barrnap_bac_cmd     = 'barrnap --kingdom bac --threads %s --reject 0.3 %s > %s/%s.barrnap_bac.gff'                                              % (num_threads, renamed_combined_fna, op_dir, op_prefix)
    shell_cmd_7_barrnap_arc_cmd     = 'barrnap --kingdom arc --threads %s --reject 0.3 %s > %s/%s.barrnap_arc.gff'                                              % (num_threads, renamed_combined_fna, op_dir, op_prefix)
    shell_cmd_8_barrnap_euk_cmd     = 'barrnap --kingdom euk --threads %s --reject 0.3 %s > %s/%s.barrnap_euk.gff'                                              % (num_threads, renamed_combined_fna, op_dir, op_prefix)
    shell_cmd_9                     = 'cut -f 1,4,5 %s/%s.barrnap_bac.gff > %s/%s.masked_barrnap_bac.bed'                                                       % (op_dir, op_prefix, op_dir, op_prefix)
    shell_cmd_10                    = 'cut -f 1,4,5 %s/%s.barrnap_arc.gff > %s/%s.masked_barrnap_arc.bed'                                                       % (op_dir, op_prefix, op_dir, op_prefix)
    shell_cmd_11                    = 'cut -f 1,4,5 %s/%s.barrnap_euk.gff > %s/%s.masked_barrnap_euk.bed'                                                       % (op_dir, op_prefix, op_dir, op_prefix)
    shell_cmd_12                    = 'cat %s/%s.masked_barrnap_bac.bed %s/%s.masked_barrnap_arc.bed %s/%s.masked_barrnap_euk.bed > %s/%s.masked_barrnap.bed'   % (op_dir, op_prefix, op_dir, op_prefix, op_dir, op_prefix, op_dir, op_prefix)
    shell_cmd_13                    = 'cat %s/%s.masked_barrnap.bed %s/%s.masked_metaxa.bed > %s/%s.masked_rrna.bed'                                            % (op_dir, op_prefix, op_dir, op_prefix, op_dir, op_prefix)

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

    ################################################# get tRNA regions #################################################

    os.system('mkdir %s' % tRNAscan_wd)

    threads_per_job = num_threads/(len(seq_file_list)*2)
    if threads_per_job == 0:
        threads_per_job = 1

    tRNAscan_cmd_list = []
    for each_gnm in seq_file_list:
        _, _, gnm_base, _ = sep_path_basename_ext(each_gnm)
        trnascan_arc_cmd = 'tRNAscan-SE -A -b %s/%s.ar.bedformat --thread %s %s'  % (tRNAscan_wd, gnm_base, threads_per_job, each_gnm)
        trnascan_bac_cmd = 'tRNAscan-SE -B -b %s/%s.bac.bedformat --thread %s %s' % (tRNAscan_wd, gnm_base, threads_per_job, each_gnm)
        tRNAscan_cmd_list.append(trnascan_arc_cmd)
        tRNAscan_cmd_list.append(trnascan_bac_cmd)

    # run tRNAscan-SE with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, tRNAscan_cmd_list)
    pool.close()
    pool.join()

    # combine the results
    shell_cmd_16 = 'cat %s/%s.ar.bedformat %s/%s.bac.bedformat > %s/%s.bedformat'   % (tRNAscan_wd, op_prefix, tRNAscan_wd, op_prefix, op_dir, op_prefix)
    shell_cmd_17 = 'cut -f 1,2,3 %s/%s.bedformat > %s/%s.masked_trnascan.bed'       % (op_dir, op_prefix, op_dir, op_prefix)
    print(shell_cmd_16)
    os.system(shell_cmd_16)
    print(shell_cmd_17)
    os.system(shell_cmd_17)

    ####################################################################################################################

    # find duplicated areas
    shell_cmd_18    = 'dustmasker -in %s -out %s/%s.lowcom.out -outfmt acclist'         % (renamed_combined_fna, op_dir, op_prefix)
    shell_cmd_19    = "cut -d '>' -f 2 %s/%s.lowcom.out > %s/%s.masked_dustmasker.bed"  % (op_dir, op_prefix, op_dir, op_prefix)
    print(shell_cmd_18)
    os.system(shell_cmd_18)
    print(shell_cmd_19)
    os.system(shell_cmd_19)

    # cover the above area with NNNN
    shell_cmd_20    = "cat %s/%s.masked_rrna.bed %s/%s.masked_trnascan.bed %s/%s.masked_dustmasker.bed > %s/%s.mask.bed"    % (op_dir, op_prefix, op_dir, op_prefix, op_dir, op_prefix, op_dir, op_prefix)
    shell_cmd_21    = "awk -F'\t' '$2>=0' %s/%s.mask.bed > %s/%s.mask_0.bed"                                                % (op_dir, op_prefix, op_dir, op_prefix)
    shell_cmd_22    = "bedtools maskfasta -fi %s -bed %s/%s.mask_0.bed -fo %s/%s.masked.fna"                                % (renamed_combined_fna, op_dir, op_prefix, op_dir, op_prefix)
    print(shell_cmd_20)
    os.system(shell_cmd_20)
    print(shell_cmd_21)
    os.system(shell_cmd_21)
    print(shell_cmd_22)
    os.system(shell_cmd_22)

    ####################################################################################################################

    print('Masked sequences exported to %s/%s.masked.fna' % (op_dir, op_prefix))
    print('Done!')


if __name__ == '__main__':

    abd_mask_parser = argparse.ArgumentParser()
    abd_mask_parser.add_argument('-i', required=True,                               help='sequence folder')
    abd_mask_parser.add_argument('-x', required=True,                               help='sequence file extension')
    abd_mask_parser.add_argument('-o', required=True,                               help='output directory')
    abd_mask_parser.add_argument('-p', required=False, default='renamed_combined',  help='output prefix, default is renamed_combined')
    abd_mask_parser.add_argument('-t', required=False, type=int, default=1,         help='number of threads, default is 1')
    abd_mask_parser.add_argument('-f', required=False, action="store_true",         help='force overwrite')
    args = vars(abd_mask_parser.parse_args())
    abd_mask(args)
