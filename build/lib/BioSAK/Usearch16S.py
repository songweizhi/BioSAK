import os
import glob
import math
import argparse
import subprocess
from Bio import SeqIO
import multiprocessing as mp


Usearch16S_usage = '''
====================================== Usearch16S example commands ======================================

# This is a wrapper for the following Usearch modules, please refer to the output cmds.txt for details:
-fastx_uniques, -unoise3, -uchime2_ref, -otutab

# classify OTU with best hit approach
BioSAK Usearch16S -i CleanData -x fna -o op_dir -t 12 -f -r ssu_all_r220.blca.fa
BioSAK Usearch16S -i CleanData -x fna -o op_dir -t 12 -f -r SILVA_138.2_SSURef_NR99_tax_silva.blca.fa

# classify OTU with BLCA
BioSAK Usearch16S -i CleanData -x fna -o op_dir -t 12 -f -blca -r ssu_all_r220.blca.fa -c ssu_all_r220.blca.tax
BioSAK Usearch16S -i CleanData -x fna -o op_dir -t 12 -f -blca -r SILVA_138.2_SSURef_NR99_tax_silva.blca.fa -c SILVA_138.2_SSURef_NR99_tax_silva.blca.tax

# GTDB SSU (BLCA compatible)
-r /Users/songweizhi/DB/BLCA/GTDB_SSU/ssu_all_r220.blca.fa
-c /Users/songweizhi/DB/BLCA/GTDB_SSU/ssu_all_r220.blca.tax

# SILVA SSU (BLCA compatible)
-r /Users/songweizhi/DB/BLCA/SILVA_SSU/SILVA_138.2_SSURef_NR99_tax_silva.blca.fa
-c /Users/songweizhi/DB/BLCA/SILVA_SSU/SILVA_138.2_SSURef_NR99_tax_silva.blca.tax

=========================================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def check_executables(program_list):

    not_detected_programs = []
    for needed_program in program_list:

        if subprocess.call(['which', needed_program], stdout=open(os.devnull, 'wb')) != 0:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not detected, program exited!' % ','.join(not_detected_programs))
        exit()


def best_hit(file_in, file_out):

    file_out_handle = open(file_out, 'w')
    best_hit_line = ''
    best_hit_query_id = ''
    best_hit_score = 0
    for blast_hit in open(file_in):
        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]
        bit_score = float(blast_hit_split[11])

        if best_hit_query_id == '':
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

        elif (query_id == best_hit_query_id) and (bit_score > best_hit_score):
            best_hit_score = bit_score
            best_hit_line = blast_hit

        elif query_id != best_hit_query_id:
            file_out_handle.write(best_hit_line)
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

    file_out_handle.write(best_hit_line)
    file_out_handle.close()


def parse_blca_op(blca_output):

    blca_op_name, blca_op_path, blca_op_base, blca_op_ext = sep_path_basename_ext(blca_output)
    output_file_1 = '%s/%s.1.txt' % (blca_op_path, blca_op_base)
    output_file_2 = '%s/%s.2.txt' % (blca_op_path, blca_op_base)

    # read in input file
    s16_taxon_blca_dict = {}
    for each_16s_taxon in open(blca_output):
        each_16s_taxon_split = each_16s_taxon.strip().split('\t')
        s16_taxon_blca_dict[each_16s_taxon_split[0]] = each_16s_taxon_split[1]

    taxon_dict_formatted_with_num = {}
    taxon_dict_formatted_no_num = {}
    for each_16s in s16_taxon_blca_dict:
        taxon_blca_raw = s16_taxon_blca_dict[each_16s]
        formatted_taxon_str_with_num = 'Unclassified'
        formatted_taxon_str_no_num = 'Unclassified'
        if taxon_blca_raw != 'Unclassified':
            split2 = taxon_blca_raw.strip().split(';')[:-1]
            formatted_taxon_list_with_num2 = []
            formatted_taxon_list_no_num2 = []
            element_index = 0
            while element_index < len(split2):
                taxon_blca = split2[element_index]
                confidence_value = split2[element_index + 1]
                confidence_value = float("{0:.2f}".format(float(confidence_value)))
                taxon_rank_blca = taxon_blca.split(':')[0]

                taxon_rank_gtdb = ''
                if taxon_rank_blca == 'superkingdom':
                    taxon_rank_gtdb = 'd'
                elif taxon_rank_blca == 'phylum':
                    taxon_rank_gtdb = 'p'
                elif taxon_rank_blca == 'class':
                    taxon_rank_gtdb = 'c'
                elif taxon_rank_blca == 'order':
                    taxon_rank_gtdb = 'o'
                elif taxon_rank_blca == 'family':
                    taxon_rank_gtdb = 'f'
                elif taxon_rank_blca == 'genus':
                    taxon_rank_gtdb = 'g'
                elif taxon_rank_blca == 'species':
                    taxon_rank_gtdb = 's'

                taxon_name_blca  = ':'.join(taxon_blca.split(':')[1:]).replace(':', '')
                current_rank_with_num = '%s__%s(%s)' % (taxon_rank_gtdb, taxon_name_blca, confidence_value)
                current_rank_no_num   = '%s__%s' % (taxon_rank_gtdb, taxon_name_blca)
                formatted_taxon_list_with_num2.append(current_rank_with_num)
                formatted_taxon_list_no_num2.append(current_rank_no_num)
                element_index += 2

            formatted_taxon_str_with_num = ';'.join(formatted_taxon_list_with_num2)
            formatted_taxon_str_no_num = ';'.join(formatted_taxon_list_no_num2)

        formatted_taxon_str_with_numno_space = '_'.join(formatted_taxon_str_with_num.split(' '))
        formatted_taxon_str_no_num_no_space = '_'.join(formatted_taxon_str_no_num.split(' '))

        taxon_dict_formatted_with_num[each_16s] = formatted_taxon_str_with_numno_space
        taxon_dict_formatted_no_num[each_16s] = formatted_taxon_str_no_num_no_space

    output_file_1_handle = open(output_file_1, 'w')
    output_file_2_handle = open(output_file_2, 'w')
    for each_seq in taxon_dict_formatted_with_num:
        output_file_1_handle.write('%s\t%s\n' % (each_seq, taxon_dict_formatted_with_num[each_seq]))
        output_file_2_handle.write('%s\t%s\n' % (each_seq, taxon_dict_formatted_no_num[each_seq]))
    output_file_1_handle.close()
    output_file_2_handle.close()


def split_fasta(fasta_in, per_file_seq_num, num_of_file, output_dir):

    fasta_in_name, fasta_in_path, fasta_in_basename, fasta_in_ext = sep_path_basename_ext(fasta_in)

    if (per_file_seq_num is None) and (num_of_file is None):
        print('Please either use the "-ns" flag to define the number of sequence per file or the "-nf" flag to specify the number of files to be generated.')
        exit()
    elif (per_file_seq_num is not None) and (num_of_file is not None):
        print('Options "-ns" and "-nf" are incompatible, only specify one of them.')
        exit()
    elif (per_file_seq_num is not None) and (num_of_file is None):

        if os.path.isdir(output_dir) is True:
            print('Output folder detected, please provide a different name!')
            exit()
        else:
            os.mkdir(output_dir)

        n = 1
        for seq_record in SeqIO.parse(fasta_in, 'fasta'):
            file_index = (n - 1) // per_file_seq_num + 1
            pwd_sub_file = '%s/%s_%s.fa' % (output_dir, fasta_in_basename, file_index)
            if per_file_seq_num == 1:
                pwd_sub_file = '%s/%s.fa' % (output_dir, seq_record.id)
            # write out sequence
            with open(pwd_sub_file, 'a') as pwd_sub_file_handle:
                pwd_sub_file_handle.write('>%s\n' % seq_record.id)
                pwd_sub_file_handle.write('%s\n'  % str(seq_record.seq))
            n += 1
    elif (per_file_seq_num is None) and (num_of_file is not None):

        if os.path.isdir(output_dir) is True:
            print('Output folder detected, please provide a different name!')
            exit()
        else:
            os.mkdir(output_dir)

        # get total number of sequences
        total_seq_num = 0
        for each_seq in SeqIO.parse(fasta_in, 'fasta'):
            total_seq_num += 1
        seq_num_per_file = math.ceil(total_seq_num/num_of_file)

        # write out
        seq_index = 1
        for seq_record in SeqIO.parse(fasta_in, 'fasta'):
            file_index = math.ceil(seq_index/seq_num_per_file)
            pwd_sub_file = '%s/%s_%s.fa' % (output_dir, fasta_in_basename, file_index)

            # write out sequence
            with open(pwd_sub_file, 'a') as pwd_sub_file_handle:
                pwd_sub_file_handle.write('>%s\n' % seq_record.id)
                pwd_sub_file_handle.write('%s\n'  % str(seq_record.seq))
            seq_index += 1

def blca(fa_dir, fa_ext, output_folder, ref_seq, ref_tax, num_threads):

    pwd_current_file  = os.path.realpath(__file__)
    current_file_path = '/'.join(pwd_current_file.split('/')[:-1])
    blca_main_py      = '%s/blca_main.py' % (current_file_path)

    seq_file_re = '%s/*.%s' % (fa_dir, fa_ext)
    seq_file_list = glob.glob(seq_file_re)

    # prepare blca commands
    blca_cmd_list = []
    for seq_file in seq_file_list:
        blca_cmd = 'python3 %s -q %s -r %s -p %s -i %s -o %s' % (blca_main_py, ref_seq, ref_tax, 1, seq_file, output_folder)
        blca_cmd_list.append(blca_cmd)

    # run blca with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, blca_cmd_list)
    pool.close()
    pool.join()

    # parse blca output
    blca_op_re = '%s/*.blca.out' % output_folder
    blca_op_list = glob.glob(blca_op_re)
    for blca_op in blca_op_list:
        parse_blca_op(blca_op)


def Usearch16S(args):

    clean_data_dir      = args['i']
    clean_data_ext      = args['x']
    op_dir              = args['o']
    run_blca            = args['blca']
    ref_seq             = args['r']
    ref_tax             = args['c']
    num_threads         = args['t']
    force_create_op_dir = args['f']
    len_range_str       = args['l']

    if len_range_str.startswith('-'):
        min_len = 0
        max_len = int(len_range_str[1:])
    elif len_range_str.endswith('-'):
        min_len = int(len_range_str[:-1])
        max_len = 9999999
    else:
        min_len = int(len_range_str.split('-')[0])
        max_len = int(len_range_str.split('-')[1])

    if run_blca is True:
        check_executables(['blastn', 'blastdbcmd', 'clustalo', 'muscle', 'usearch'])

    # define file name
    cmd_txt                     = '%s/cmds.txt'                                 % op_dir
    dir_clean_data_renamed      = '%s/s01_CleanData_renamed_%s_%sbp'            % (op_dir, min_len, max_len)
    dir_DereplicatedData        = '%s/s02_DereplicatedData'                     % op_dir
    unoise_nc_fa                = '%s/s06_AllSamples_unoise_nc.fasta'           % op_dir
    unoise_nc_fa_split_dir      = '%s/s06_AllSamples_unoise_nc_split'           % op_dir
    s07_OtuTable                = '%s/s07_AllSamples_unoise_otu_table.txt'      % op_dir
    s09_FinalOtuTable           = '%s/s09_AllSamples_unoise_otu_tax.txt'        % op_dir
    s10_blca_op_dir             = '%s/s10_AllSamples_BLCA_classifications'      % op_dir
    combined_blca_blastn_tmp    = '%s/s08_AllSamples_unoise_nc.blca.blastn.tmp' % op_dir
    combined_blca_blastn_sorted = '%s/s08_AllSamples_unoise_nc.blca.blastn'     % op_dir
    combined_blca_out_tmp       = '%s/s08_AllSamples_unoise_nc.blca.out.tmp'    % op_dir
    combined_blca_out_sorted    = '%s/s08_AllSamples_unoise_nc.blca.out'        % op_dir
    combined_blca_1_txt_tmp     = '%s/s08_AllSamples_unoise_nc.blca.1.txt.tmp'  % op_dir
    combined_blca_1_txt_sorted  = '%s/s08_AllSamples_unoise_nc.blca.1.txt'      % op_dir
    combined_blca_2_txt_tmp     = '%s/s08_AllSamples_unoise_nc.blca.2.txt.tmp'  % op_dir
    combined_blca_2_txt_sorted  = '%s/s08_AllSamples_unoise_nc.blca.2.txt'      % op_dir

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % dir_clean_data_renamed)
    os.system('mkdir %s' % dir_DereplicatedData)

    clean_data_re = '%s/*.%s' % (clean_data_dir, clean_data_ext)
    clean_data_list = glob.glob(clean_data_re)
    if len(clean_data_list) == 0:
        print('Input sequences not found, program exited!')
        exit()

    file_index = 1
    for each_file in sorted(clean_data_list):

        fa_name, _, fa_base, _  = sep_path_basename_ext(each_file)
        pwd_fa_renamed          = '%s/%s'                           % (dir_clean_data_renamed, fa_name)
        pwd_fa_uniques          = '%s/%s_uniques.fasta'             % (dir_DereplicatedData, fa_base)

        print('Running usearch -fastx_uniques on %s/%s: %s' % (file_index, len(clean_data_list), fa_base))
        file_index += 1

        # add sample id to contig id
        pwd_fa_renamed_handle = open(pwd_fa_renamed, 'w')
        for each_seq in SeqIO.parse(each_file, 'fasta'):
            seq_id = each_seq.id
            seq_len = len(each_seq.seq)
            if max_len >= seq_len >= min_len:
                pwd_fa_renamed_handle.write('>%s;sample=%s;\n' % (seq_id, fa_base))
                pwd_fa_renamed_handle.write(str(each_seq.seq) + '\n')
        pwd_fa_renamed_handle.close()

        # run dereplicate
        cmd_fastx_uniques = 'usearch -fastx_uniques %s -fastaout %s -sizeout &>/dev/null' % (pwd_fa_renamed, pwd_fa_uniques)
        with open(cmd_txt, 'a') as cmd_txt_handle:
            cmd_txt_handle.write(cmd_fastx_uniques + '\n')
        os.system(cmd_fastx_uniques)

    # combine sequences from all samples
    cmd_2 = 'cat %s/*.fasta > %s/s03_AllSamples.fasta' % (dir_DereplicatedData, op_dir)
    with open(cmd_txt, 'a') as cmd_txt_handle:
        cmd_txt_handle.write(cmd_2 + '\n')
    print(cmd_2)
    os.system(cmd_2)

    # dereplicate AllSamples.fasta
    cmd_3  = 'usearch -fastx_uniques %s/s03_AllSamples.fasta -fastaout %s/s04_AllSamples_uniques.fasta -sizein -sizeout -strand both &>/dev/null' % (op_dir, op_dir)
    with open(cmd_txt, 'a') as cmd_txt_handle:
        cmd_txt_handle.write(cmd_3 + '\n')
    print(cmd_3)
    os.system(cmd_3)

    # Generating unique sequences using UNOISE
    cmd_4  = 'usearch -unoise3 %s/s04_AllSamples_uniques.fasta -zotus %s/s05_AllSamples_denoised.fasta &>/dev/null' % (op_dir, op_dir)
    with open(cmd_txt, 'a') as cmd_txt_handle:
        cmd_txt_handle.write(cmd_4 + '\n')
    print(cmd_4)
    os.system(cmd_4)

    # Chimera Removal
    cmd_5  = 'usearch -uchime2_ref %s/s05_AllSamples_denoised.fasta -db %s -strand plus -mode high_confidence -notmatched %s/s06_AllSamples_unoise_nc.fasta &>/dev/null' % (op_dir, ref_seq, op_dir)
    if run_blca is True:
        cmd_5  = 'usearch -uchime2_ref %s/s05_AllSamples_denoised.fasta -db %s -strand plus -mode high_confidence -notmatched %s/s06_AllSamples_unoise_nc.fasta &>/dev/null' % (op_dir, ref_seq, op_dir)
    with open(cmd_txt, 'a') as cmd_txt_handle:
        cmd_txt_handle.write(cmd_5 + '\n')
    print(cmd_5)
    os.system(cmd_5)

    # generate OTU table, an OTU table is made by the otutab command
    cmd_6 = 'usearch -otutab %s/s03_AllSamples.fasta -db %s/s06_AllSamples_unoise_nc.fasta -id 0.97 -otutabout %s/s07_AllSamples_unoise_otu_table.txt &>/dev/null' % (op_dir, op_dir, op_dir)
    with open(cmd_txt, 'a') as cmd_txt_handle:
        cmd_txt_handle.write(cmd_6 + '\n')
    print(cmd_6)
    os.system(cmd_6)
    print('usearch -otutab finished')

    # Mapping of OTUs on Reference Database
    if run_blca is False:
        cmd_7  = 'blastn -query %s/s06_AllSamples_unoise_nc.fasta -outfmt 6 -out %s/s08_AllSamples_unoise_nc.txt -db %s -evalue 1e-20 -num_threads %s' % (op_dir, op_dir, ref_seq, num_threads)
        with open(cmd_txt, 'a') as cmd_txt_handle:
            cmd_txt_handle.write(cmd_7 + '\n')
        print(cmd_7)
        os.system(cmd_7)

        # keep best hit
        blast_op          = '%s/s08_AllSamples_unoise_nc_blastn.txt'          % op_dir
        blast_op_best_hit = '%s/s08_AllSamples_unoise_nc_blastn_best_hit.txt' % op_dir
        best_hit(blast_op, blast_op_best_hit)

        # Merge Table and Taxonomy
        silva_ref_to_tax_dict = dict()
        for each_seq in SeqIO.parse(ref_seq, 'fasta'):
            seq_id = each_seq.id
            seq_tax = ' '.join(each_seq.description.split(' ')[1:])
            silva_ref_to_tax_dict[seq_id] = seq_tax

        otu_best_hit_dict = dict()
        for each in open(blast_op_best_hit):
            each_split = each.strip().split('\t')
            info_list = [each_split[1], each_split[2], each_split[3], each_split[10]]
            otu_best_hit_dict[each_split[0]] = info_list

        s09_FinalOtuTable_handle = open(s09_FinalOtuTable, 'w')
        for each_line in open(s07_OtuTable):
            if each_line.startswith('#'):
                s09_FinalOtuTable_handle.write('ID\tIdentity\tAlignment_Length\tEvalue\tTaxonomy\n')
            else:
                each_line_split = each_line.strip().split('\t')
                otu_id = each_line_split[0]
                best_hit_info_list = otu_best_hit_dict.get(otu_id, ['na', 'na', 'na', 'na'])
                best_hit_info_list_without_ref_id = []
                best_hit_info_list_without_ref_id.append(best_hit_info_list[1])
                best_hit_info_list_without_ref_id.append(best_hit_info_list[2])
                best_hit_info_list_without_ref_id.append(best_hit_info_list[3])
                best_hit_info_str = '\t'.join(best_hit_info_list_without_ref_id)
                ref_id = best_hit_info_list[0]
                ref_tax = silva_ref_to_tax_dict.get(ref_id, 'na')
                s09_FinalOtuTable_handle.write('%s\t%s\t%s\n' % (otu_id, best_hit_info_str, ref_tax))
        s09_FinalOtuTable_handle.close()

    else:
        # split fasta file
        print('split fasta')
        split_fasta(unoise_nc_fa, None, num_threads, unoise_nc_fa_split_dir)

        # create blca output directory
        os.system('mkdir %s' % s10_blca_op_dir)

        # run blca
        print('running blca')
        blca(unoise_nc_fa_split_dir, 'fa', s10_blca_op_dir, ref_seq, ref_tax, num_threads)

        # combine results
        os.system('cat %s/*.blastn > %s'       % (s10_blca_op_dir, combined_blca_blastn_tmp))
        os.system('cat %s/*.blca.out > %s'     % (s10_blca_op_dir, combined_blca_out_tmp))
        os.system('cat %s/*.blca.1.txt > %s'   % (s10_blca_op_dir, combined_blca_1_txt_tmp))
        os.system('cat %s/*.blca.2.txt > %s'   % (s10_blca_op_dir, combined_blca_2_txt_tmp))

        # sort results
        os.system('cat %s | sort > %s' % (combined_blca_blastn_tmp, combined_blca_blastn_sorted))
        os.system('cat %s | sort > %s' % (combined_blca_out_tmp, combined_blca_out_sorted))
        os.system('cat %s | sort > %s' % (combined_blca_1_txt_tmp, combined_blca_1_txt_sorted))
        os.system('cat %s | sort > %s' % (combined_blca_2_txt_tmp, combined_blca_2_txt_sorted))

        # remove tmp files
        os.system('rm -r %s' % s10_blca_op_dir)
        os.system('rm %s'    % combined_blca_blastn_tmp)
        os.system('rm %s'    % combined_blca_out_tmp)
        os.system('rm %s'    % combined_blca_1_txt_tmp)
        os.system('rm %s'    % combined_blca_2_txt_tmp)

    # final report
    print('OTU table exported to: %s' % s07_OtuTable)
    if run_blca is False:
        print('OTU taxonomy exported to: %s' % s09_FinalOtuTable)
    else:
        print('OTU taxonomy exported to: %s and %s' % (combined_blca_1_txt_sorted, combined_blca_2_txt_sorted))
    print('Done')


if __name__ == '__main__':

    Usearch16S_parser = argparse.ArgumentParser(usage=Usearch16S_usage)
    Usearch16S_parser.add_argument('-i',       required=True,                        help='path to input sequences')
    Usearch16S_parser.add_argument('-x',       required=True,                        help='file extension')
    Usearch16S_parser.add_argument('-o',       required=True,                        help='output directory')
    Usearch16S_parser.add_argument('-l',       required=False,                       help='allowed length range, default is 0-9999999')
    Usearch16S_parser.add_argument('-r',       required=False,                       help='reference sequences')
    Usearch16S_parser.add_argument('-blca',    required=False, action="store_true",  help='perform classification with BLCA')
    Usearch16S_parser.add_argument('-c',       required=False,                       help='reference sequence taxonomy')
    Usearch16S_parser.add_argument('-t',       required=False, type=int, default=1,  help='number of threads, default is 1')
    Usearch16S_parser.add_argument('-f',       required=False, action="store_true",  help='force overwrite')
    args = vars(Usearch16S_parser.parse_args())
    Usearch16S(args)
