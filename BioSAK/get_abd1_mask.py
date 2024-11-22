import os


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


########################################################################################################################

fasta_file = 'symbionts_dRep99_sponge_coral_66.fna'


########################################################################################################################

_, _, f_base, _ = sep_path_basename_ext(fasta_file)

# get rRNA regions
shell_cmd_1_metaxa2_ssu_cmd     = 'metaxa2 --plus --mode m --cpu 12 --multi_thread T --table T -g ssu --not_found T -i %s -o %s.metaxa2_ssu'    % (fasta_file, f_base)
shell_cmd_2_metaxa2_lsu_cmd     = 'metaxa2 --plus --mode m --cpu 12 --multi_thread T --table T -g lsu --not_found T -i %s -o %s.metaxa2_lsu'    % (fasta_file, f_base)
shell_cmd_3                     = 'cut -f 1,9,10 %s.metaxa2_ssu.extraction.results > %s.masked_metaxa_ssu.bed'                                  % (f_base, f_base)
shell_cmd_4                     = 'cut -f 1,9,10 %s.metaxa2_lsu.extraction.results > %s.masked_metaxa_lsu.bed'                                  % (f_base, f_base)
shell_cmd_5                     = 'cat %s.masked_metaxa_ssu.bed %s.masked_metaxa_lsu.bed > %s.masked_metaxa.bed'                                % (f_base, f_base, f_base)
shell_cmd_6_barrnap_bac_cmd     = 'barrnap --kingdom bac --threads 12 --reject 0.3 %s > %s.barrnap_bac.gff'                                     % (fasta_file, f_base)
shell_cmd_7_barrnap_arc_cmd     = 'barrnap --kingdom arc --threads 12 --reject 0.3 %s > %s.barrnap_arc.gff'                                     % (fasta_file, f_base)
shell_cmd_8_barrnap_euk_cmd     = 'barrnap --kingdom euk --threads 12 --reject 0.3 %s > %s.barrnap_euk.gff'                                     % (fasta_file, f_base)
shell_cmd_9                     = 'cut -f 1,4,5 %s.barrnap_bac.gff > %s.masked_barrnap_bac.bed'                                                 % (f_base, f_base)
shell_cmd_10                    = 'cut -f 1,4,5 %s.barrnap_arc.gff > %s.masked_barrnap_arc.bed'                                                 % (f_base, f_base)
shell_cmd_11                    = 'cut -f 1,4,5 %s.barrnap_euk.gff > %s.masked_barrnap_euk.bed'                                                 % (f_base, f_base)
shell_cmd_12                    = 'cat %s.masked_barrnap_bac.bed %s.masked_barrnap_arc.bed %s.masked_barrnap_euk.bed > %s.masked_barrnap.bed'   % (f_base, f_base, f_base, f_base)
shell_cmd_13                    = 'cat %s.masked_barrnap.bed %s.masked_metaxa.bed > %s.masked_rrna.bed'                                         % (f_base, f_base, f_base)

# get tRNA regions
shell_cmd_14_trna_scan_arc_cmd  = 'tRNAscan-SE -A -b %s.ar.bedformat --thread 36 %s'                                                            % (f_base, fasta_file)
shell_cmd_15_trna_scan_bac_cmd  = 'tRNAscan-SE -B -b %s.bac.bedformat --thread 36 %s'                                                           % (f_base, fasta_file)
shell_cmd_16                    = 'cat %s.ar.bedformat %s.bac.bedformat > %s.bedformat'                                                         % (f_base, f_base, f_base)
shell_cmd_17                    = 'cut -f 1,2,3 %s.bedformat > %s.masked_trnascan.bed'                                                          % (f_base, f_base)

# find duplicated areas
shell_cmd_18                    = 'dustmasker -in %s -out %s.lowcom.out -outfmt acclist'                                                        % (fasta_file, f_base)
shell_cmd_19                    = "cut -d '>' -f 2 %s.lowcom.out > %s.masked_dustmasker.bed"                                                    % (f_base, f_base)

# cover the above area with NNNN
shell_cmd_20                    = "cat %s.masked_rrna.bed %s.masked_trnascan.bed %s.masked_dustmasker.bed > %s.mask.bed"                        % (f_base, f_base, f_base, f_base)
shell_cmd_21                    = "awk -F'\t' '$2>=0' %s.mask.bed > %s.mask_0.bed"                                                              % (f_base, f_base)
shell_cmd_22                    = "bedtools maskfasta -fi %s -bed %s.mask_0.bed -fo %s.masked.fasta"                                            % (fasta_file, f_base, f_base)

# run cmds
os.system(shell_cmd_1_metaxa2_ssu_cmd)
os.system(shell_cmd_2_metaxa2_lsu_cmd)
os.system(shell_cmd_3)
os.system(shell_cmd_4)
os.system(shell_cmd_5)
os.system(shell_cmd_6_barrnap_bac_cmd)
os.system(shell_cmd_7_barrnap_arc_cmd)
os.system(shell_cmd_8_barrnap_euk_cmd)
os.system(shell_cmd_9)
os.system(shell_cmd_10)
os.system(shell_cmd_11)
os.system(shell_cmd_12)
os.system(shell_cmd_13)
os.system(shell_cmd_14_trna_scan_arc_cmd)
os.system(shell_cmd_15_trna_scan_bac_cmd)
os.system(shell_cmd_16)
os.system(shell_cmd_17)
os.system(shell_cmd_18)
os.system(shell_cmd_19)
os.system(shell_cmd_20)
os.system(shell_cmd_21)
os.system(shell_cmd_22)
