import os
import argparse


get_abd1_mask_usage = '''
================================= get_abd1_mask example commands =================================

BioSAK get_abd1_mask -r symbionts_dRep99_sponge_coral_66.fna -o get_abd1_mask_wd

==================================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def get_abd1_mask(args):

    fasta_file      = args['r']
    op_dir          = args['o']
    force_overwrite = args['f']
    num_threads     = args['t']

    # create op dir
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output directory detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    # fasta_file = 'symbionts_dRep99_sponge_coral_66.fna'
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


if __name__ == '__main__':

    abund_parser = argparse.ArgumentParser()
    abund_parser.add_argument('-r', required=True,                          help='reference sequence')
    abund_parser.add_argument('-o', required=True,                          help='output directory')
    abund_parser.add_argument('-t', required=False, type=int, default=1,    help='number of threads, default is 1')
    abund_parser.add_argument('-f', required=False, action="store_true",    help='force overwrite')
    args = vars(abund_parser.parse_args())
    get_abd1_mask(args)
