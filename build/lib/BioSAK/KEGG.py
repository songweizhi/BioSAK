import os
import glob
import argparse
from Bio import SeqIO
from time import sleep
import multiprocessing as mp
from datetime import datetime
from BioSAK.global_functions import time_format
from BioSAK.global_functions import force_create_folder
from BioSAK.global_functions import sep_path_basename_ext
from BioSAK.global_functions import get_gene_list_TotalDepth
from BioSAK.global_functions import AnnotateNorm


KEGG_parser_usage = '''
======================================== KEGG example commands =======================================

# Dependencies: blast+ or diamond

# annotation with NCBI blastp (default, for small dataset)
BioSAK KEGG -db_dir path/to/your/KEGG_db_dir -t 6 -seq_in input.faa -d input.depth

# annotation with Diamond blastp (for big dataset)
BioSAK KEGG -db_dir path/to/your/KEGG_db_dir -t 12 -seq_in faa_folder -x faa -d depth_files -diamond

# get summary for BlastKOALA/GhostKOALA produced results
BioSAK KEGG -db_dir path/to/your/KEGG_db_dir -t 9 -ko_in user_ko.txt
BioSAK KEGG -db_dir path/to/your/KEGG_db_dir -t 9 -ko_in user_ko_folder -x txt

# Prepare DB files, you need to have the following three files in your KEGG_db_dir:
1. Sequence file, only needed for "-seq_in" mode, DECOMPRESS and RENAME it to kegg_db_seq.fasta
   e.g. prokaryotes.pep.gz (https://www.kegg.jp/kegg/download/Readme/README.fasta)
2. seq2ko file, only needed for "-seq_in" mode, DECOMPRESS and RENAME it to kegg_db_seq2ko.txt
   e.g. prokaryotes.dat.gz (https://www.kegg.jp/kegg/download/Readme/README.fasta)
3. ko00001.keg
   https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir=

# seq2ko file format(tab separated)
aaa:Acav_4596	K05375
zpr:ZPR_0691	K21572

# How it works:
1. KEGG module uses Blast+/Diamond to get the best hits of query genes in the database with user defined e-value cutoff (default 0.001).
2. The TotalDepth of a KO is calculated by summing up the depth of all genes assigned to it.
3. The percentage of GeneNumber/TotalDepth of genes assigned to a KO is calculated by dividing them 
   by the total number/depth of genes with KO assignment (default) or by all genes in a genome ("-pct_by_all"). 

# Note!!!
1. If you run KEGG annotation for multiple files in a batch manner and want to have their depth info incorporated into the results, 
   you need to provide a folder containing individual depth files for each of your input sequence file.
   Name of the depth file needs to be exactly the same as its corresponding sequence file, except the extension which is ".depth".
2. Diamond requires quite a lot of memory for sequence comparison, especially for huge db file (e.g. KEGG db).
   Remember to request sufficient memory (e.g. 90 or 120gb) in your job script and specify a small number (e.g. -t 6) 
   of jobs executing in parallel. Otherwise, you may see some of your query genomes with no gene been annotated.

# Depth file format (one gene per line, tab separated)
gene_1	30
gene_2	10.58

# To do:
1. level C stats: separate stats for Pathway, Brite and the rests

======================================================================================================
'''


def keep_blast_hit_with_highest_bit_score(file_in, file_out):
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


def run_blast_worker(argument_list):

    pwd_input_file =        argument_list[0]
    run_blast =             argument_list[1]
    run_diamond =           argument_list[2]
    KEGG_DB_seq =           argument_list[3]
    KEGG_DB_seq_diamond =   argument_list[4]
    op_dir =                argument_list[5]
    evalue_cutoff =         argument_list[6]
    threads_num =           argument_list[7]

    ################################################### define file name ###################################################

    input_file_path, in_file_basename, input_file_ext = sep_path_basename_ext(pwd_input_file)

    blast_results =           '%s/%s_KEGG_wd/%s_blast.tab'              % (op_dir, in_file_basename, in_file_basename)
    blast_results_best_hit =  '%s/%s_KEGG_wd/%s_blast_best_hits.tab'    % (op_dir, in_file_basename, in_file_basename)

    # create output folder
    force_create_folder('%s/%s_KEGG_wd' % (op_dir, in_file_basename))


    ########################################## blast against KEGG database (Shan) ##########################################

    if run_blast is True:

        if run_diamond is False:
            blastp_cmd = 'blastp -query %s -db %s -out %s -outfmt 6 -evalue %s -num_alignments 10 -num_threads %s' % (pwd_input_file, KEGG_DB_seq, blast_results, evalue_cutoff, threads_num)
            os.system(blastp_cmd)

        else:
            diamond_cmd = 'diamond blastp -q %s --db %s --out %s --outfmt 6 --evalue %s --block-size 1 --threads %s --quiet' % (pwd_input_file, KEGG_DB_seq_diamond, blast_results, evalue_cutoff, threads_num)
            os.system(diamond_cmd)

        # only keep the best hit
        keep_blast_hit_with_highest_bit_score(blast_results, blast_results_best_hit)


def write_out_stats_GeneNumber(identified_ko_list, ko_to_gene_member_dict, ko_description_dict, stats_file_GeneNumber):

    stats_file_GeneNumber_handle = open(stats_file_GeneNumber, 'w')
    stats_file_GeneNumber_handle.write('KO\tGeneNumber\tDescription\n')
    for ko in identified_ko_list:
        ko_GeneNumber = len(ko_to_gene_member_dict[ko])
        stats_file_GeneNumber_handle.write('%s\t%s\t%s\n' % (ko[2:], ko_GeneNumber, ko_description_dict[ko[2:]]))
    stats_file_GeneNumber_handle.close()


def write_out_stats_TotalDepth(identified_ko_list, ko_to_gene_member_dict, gene_depth_dict, ko_description_dict, stats_file_TotalDepth):

    stats_file_TotalDepth_handle = open(stats_file_TotalDepth, 'w')
    stats_file_TotalDepth_handle.write('KO\tTotalDepth\tDescription\n')
    for ko in identified_ko_list:
        ko_gene_total_depth = 0
        for each_gene in ko_to_gene_member_dict[ko]:
            each_gene_depth = gene_depth_dict[each_gene]
            ko_gene_total_depth += each_gene_depth
        ko_TotalDepth = float("{0:.2f}".format(ko_gene_total_depth))
        stats_file_TotalDepth_handle.write('%s\t%s\t%s\n' % (ko[2:], ko_TotalDepth, ko_description_dict[ko[2:]]))
    stats_file_TotalDepth_handle.close()


def parse_blast_op_worker(argument_list):

    pwd_input_file =        argument_list[0]
    run_blast =             argument_list[1]
    As_description_dict =   argument_list[2]
    Bs_description_dict =   argument_list[3]
    Cs_description_dict =   argument_list[4]
    Ds_description_dict =   argument_list[5]
    D2ABCD_dict =           argument_list[6]
    db_seq_to_KO_dict =     argument_list[7]
    op_dir =                argument_list[8]
    depth_file =            argument_list[9]
    pct_by_all =            argument_list[10]

    ################################################### define file name ###################################################

    input_file_path, in_file_basename, input_file_ext = sep_path_basename_ext(pwd_input_file)

    blast_results_best_hit =             '%s/%s_KEGG_wd/%s_blast_best_hits.tab'             % (op_dir, in_file_basename, in_file_basename)
    KO_assignment_file_D =               '%s/%s_KEGG_wd/%s_ko_assignment_D.txt'             % (op_dir, in_file_basename, in_file_basename)
    KO_assignment_file_DCBA =            '%s/%s_KEGG_wd/%s_ko_assignment_ABCD.txt'          % (op_dir, in_file_basename, in_file_basename)
    stats_file_A_GeneNumber =            '%s/%s_KEGG_wd/%s_ko_stats_A.txt'                  % (op_dir, in_file_basename, in_file_basename)
    stats_file_B_GeneNumber =            '%s/%s_KEGG_wd/%s_ko_stats_B.txt'                  % (op_dir, in_file_basename, in_file_basename)
    stats_file_C_GeneNumber =            '%s/%s_KEGG_wd/%s_ko_stats_C.txt'                  % (op_dir, in_file_basename, in_file_basename)
    stats_file_D_GeneNumber =            '%s/%s_KEGG_wd/%s_ko_stats_D.txt'                  % (op_dir, in_file_basename, in_file_basename)
    stats_file_A_TotalDepth =            '%s/%s_KEGG_wd/%s_ko_stats_A_depth.txt'            % (op_dir, in_file_basename, in_file_basename)
    stats_file_B_TotalDepth =            '%s/%s_KEGG_wd/%s_ko_stats_B_depth.txt'            % (op_dir, in_file_basename, in_file_basename)
    stats_file_C_TotalDepth =            '%s/%s_KEGG_wd/%s_ko_stats_C_depth.txt'            % (op_dir, in_file_basename, in_file_basename)
    stats_file_D_TotalDepth =            '%s/%s_KEGG_wd/%s_ko_stats_D_depth.txt'            % (op_dir, in_file_basename, in_file_basename)
    stats_file_A_GeneNumber_pct =        '%s/%s_KEGG_wd/%s_ko_stats_A_pct.txt'              % (op_dir, in_file_basename, in_file_basename)
    stats_file_B_GeneNumber_pct =        '%s/%s_KEGG_wd/%s_ko_stats_B_pct.txt'              % (op_dir, in_file_basename, in_file_basename)
    stats_file_C_GeneNumber_pct =        '%s/%s_KEGG_wd/%s_ko_stats_C_pct.txt'              % (op_dir, in_file_basename, in_file_basename)
    stats_file_D_GeneNumber_pct =        '%s/%s_KEGG_wd/%s_ko_stats_D_pct.txt'              % (op_dir, in_file_basename, in_file_basename)
    stats_file_A_TotalDepth_pct =        '%s/%s_KEGG_wd/%s_ko_stats_A_depth_pct.txt'        % (op_dir, in_file_basename, in_file_basename)
    stats_file_B_TotalDepth_pct =        '%s/%s_KEGG_wd/%s_ko_stats_B_depth_pct.txt'        % (op_dir, in_file_basename, in_file_basename)
    stats_file_C_TotalDepth_pct =        '%s/%s_KEGG_wd/%s_ko_stats_C_depth_pct.txt'        % (op_dir, in_file_basename, in_file_basename)
    stats_file_D_TotalDepth_pct =        '%s/%s_KEGG_wd/%s_ko_stats_D_depth_pct.txt'        % (op_dir, in_file_basename, in_file_basename)
    stats_file_A_GeneNumber_pct_by_all = '%s/%s_KEGG_wd/%s_ko_stats_A_pct_by_all.txt'       % (op_dir, in_file_basename, in_file_basename)
    stats_file_B_GeneNumber_pct_by_all = '%s/%s_KEGG_wd/%s_ko_stats_B_pct_by_all.txt'       % (op_dir, in_file_basename, in_file_basename)
    stats_file_C_GeneNumber_pct_by_all = '%s/%s_KEGG_wd/%s_ko_stats_C_pct_by_all.txt'       % (op_dir, in_file_basename, in_file_basename)
    stats_file_D_GeneNumber_pct_by_all = '%s/%s_KEGG_wd/%s_ko_stats_D_pct_by_all.txt'       % (op_dir, in_file_basename, in_file_basename)
    stats_file_A_TotalDepth_pct_by_all = '%s/%s_KEGG_wd/%s_ko_stats_A_depth_pct_by_all.txt' % (op_dir, in_file_basename, in_file_basename)
    stats_file_B_TotalDepth_pct_by_all = '%s/%s_KEGG_wd/%s_ko_stats_B_depth_pct_by_all.txt' % (op_dir, in_file_basename, in_file_basename)
    stats_file_C_TotalDepth_pct_by_all = '%s/%s_KEGG_wd/%s_ko_stats_C_depth_pct_by_all.txt' % (op_dir, in_file_basename, in_file_basename)
    stats_file_D_TotalDepth_pct_by_all = '%s/%s_KEGG_wd/%s_ko_stats_D_depth_pct_by_all.txt' % (op_dir, in_file_basename, in_file_basename)

    ################################################# parse blast results ##################################################

    if run_blast is True:

        # store blast results in dict
        query_to_db_seq_dict = {}
        for each_query in open(blast_results_best_hit):
            each_query_split = each_query.strip().split('\t')
            query_id = each_query_split[0]
            db_seq = each_query_split[1]
            query_to_db_seq_dict[query_id] = db_seq

        # get all query sequence id
        query_seq_id_list = []
        for each_seq in SeqIO.parse(pwd_input_file, 'fasta'):
            query_seq_id_list.append(str(each_seq.id))

        # get ko id at level D for all query genes
        KO_assignment_file_handle = open(KO_assignment_file_D, 'w')
        for each_query_seq in sorted(query_seq_id_list):

            if each_query_seq in query_to_db_seq_dict:
                db_hit_id = query_to_db_seq_dict[each_query_seq]

                if db_hit_id in db_seq_to_KO_dict:
                    db_hit_id_ko = db_seq_to_KO_dict[db_hit_id]
                    if ',' in db_hit_id_ko:
                        db_hit_id_ko_split = db_hit_id_ko.split(',')
                        for each_db_hit_id_ko in db_hit_id_ko_split:
                            KO_assignment_file_handle.write('%s\t%s\n' % (each_query_seq, each_db_hit_id_ko))
                    else:
                        KO_assignment_file_handle.write('%s\t%s\n' % (each_query_seq, db_hit_id_ko))
                else:
                    KO_assignment_file_handle.write('%s\n' % (each_query_seq))
            else:
                KO_assignment_file_handle.write('%s\n' % (each_query_seq))
        KO_assignment_file_handle.close()

    else:
        KO_assignment_file_D = pwd_input_file

    # get ko id at all levels for all query genes
    ko_assign_ABCD_handle = open(KO_assignment_file_DCBA, 'w')
    ko_assign_ABCD_handle.write('Gene_id\tko_A\tko_B\tko_C\tko_D\tDesc_A\tDesc_B\tDesc_C\tDesc_D\n')
    query_seq_id_all = set()
    genes_with_ko = set()
    for query_gene in open(KO_assignment_file_D):
        query_gene_split = query_gene.strip().split('\t')
        gene_ID = query_gene_split[0]

        if len(query_gene_split) == 1:
            query_seq_id_all.add(query_gene_split[0])
            ko_assign_ABCD_handle.write('%s\n' % gene_ID)

        if len(query_gene_split) == 2:
            query_seq_id_all.add(query_gene_split[0])
            genes_with_ko.add(query_gene_split[0])
            KO_ID = query_gene_split[1]
            if KO_ID in D2ABCD_dict:
                KO_ID_ABCD = D2ABCD_dict[KO_ID]

                if len(KO_ID_ABCD) == 1:
                    KO_DCBA_list = KO_ID_ABCD[0].split('|')[::-1]
                    KO_DCBA_list_only_id = [i.split('_')[1] for i in KO_DCBA_list]
                    desc_A = As_description_dict[KO_DCBA_list_only_id[3]]
                    desc_B = Bs_description_dict[KO_DCBA_list_only_id[2]]
                    desc_C = Cs_description_dict[KO_DCBA_list_only_id[1]]
                    desc_D = Ds_description_dict[KO_DCBA_list_only_id[0]]
                    ko_assign_ABCD_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (gene_ID, '\t'.join(KO_DCBA_list[::-1]), desc_A, desc_B, desc_C, desc_D))

                if len(KO_ID_ABCD) > 1:
                    for each_ABCD in KO_ID_ABCD:
                        each_KO_DCBA_list = each_ABCD.split('|')[::-1]
                        each_KO_DCBA_list_only_id = [i.split('_')[1] for i in each_KO_DCBA_list]
                        each_desc_A = As_description_dict[each_KO_DCBA_list_only_id[3]]
                        each_desc_B = Bs_description_dict[each_KO_DCBA_list_only_id[2]]
                        each_desc_C = Cs_description_dict[each_KO_DCBA_list_only_id[1]]
                        each_desc_D = Ds_description_dict[each_KO_DCBA_list_only_id[0]]
                        ko_assign_ABCD_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (gene_ID, '\t'.join(each_KO_DCBA_list[::-1]), each_desc_A, each_desc_B, each_desc_C, each_desc_D))

    ko_assign_ABCD_handle.close()

    ##################################################### Get summary ######################################################

    # read in depth info
    gene_depth_dict = {}
    if depth_file is not None:
        for each_depth in open(depth_file):
            each_depth_split = each_depth.strip().split('\t')
            gene_depth_dict[each_depth_split[0]] = float(each_depth_split[1])

    # get total number and depth of all genes in one file
    total_depth_for_all_query_genes = 0
    genes_with_ko_TotalDepth = 0
    if depth_file is not None:
        for gene in query_seq_id_all:
            gene_depth = gene_depth_dict[gene]
            total_depth_for_all_query_genes += gene_depth

        genes_with_ko_TotalDepth = get_gene_list_TotalDepth(genes_with_ko, gene_depth_dict)

    identified_ko_A_list = []
    identified_ko_B_list = []
    identified_ko_C_list = []
    identified_ko_D_list = []
    ko_A_to_gene_member_dict = {}
    ko_B_to_gene_member_dict = {}
    ko_C_to_gene_member_dict = {}
    ko_D_to_gene_member_dict = {}
    ko_NA_to_gene_member_list = []
    for each_query in open(KO_assignment_file_DCBA):
        if not each_query.startswith('Gene_id'):
            each_query_split = each_query.strip().split('\t')
            query_id = each_query_split[0]

            if len(each_query_split) == 1:
                ko_NA_to_gene_member_list.append(query_id)

            if len(each_query_split) > 1:
                query_ko_A = each_query_split[1]
                query_ko_B = each_query_split[2]
                query_ko_C = each_query_split[3]
                query_ko_D = each_query_split[4]

                if query_ko_A not in identified_ko_A_list:
                    identified_ko_A_list.append(query_ko_A)
                if query_ko_B not in identified_ko_B_list:
                    identified_ko_B_list.append(query_ko_B)
                if query_ko_C not in identified_ko_C_list:
                    identified_ko_C_list.append(query_ko_C)
                if query_ko_D not in identified_ko_D_list:
                    identified_ko_D_list.append(query_ko_D)

                if query_ko_A not in ko_A_to_gene_member_dict:
                    ko_A_to_gene_member_dict[query_ko_A] = [query_id]
                else:
                    if query_id not in ko_A_to_gene_member_dict[query_ko_A]:
                        ko_A_to_gene_member_dict[query_ko_A].append(query_id)

                if query_ko_B not in ko_B_to_gene_member_dict:
                    ko_B_to_gene_member_dict[query_ko_B] = [query_id]
                else:
                    if query_id not in ko_B_to_gene_member_dict[query_ko_B]:
                        ko_B_to_gene_member_dict[query_ko_B].append(query_id)

                if query_ko_C not in ko_C_to_gene_member_dict:
                    ko_C_to_gene_member_dict[query_ko_C] = [query_id]
                else:
                    if query_id not in ko_C_to_gene_member_dict[query_ko_C]:
                        ko_C_to_gene_member_dict[query_ko_C].append(query_id)

                if query_ko_D not in ko_D_to_gene_member_dict:
                    ko_D_to_gene_member_dict[query_ko_D] = [query_id]
                else:
                    if query_id not in ko_D_to_gene_member_dict[query_ko_D]:
                        ko_D_to_gene_member_dict[query_ko_D].append(query_id)

    #################### write out GeneNumber and TotalDepth stats ####################

    write_out_stats_GeneNumber(identified_ko_A_list, ko_A_to_gene_member_dict, As_description_dict, stats_file_A_GeneNumber)
    write_out_stats_GeneNumber(identified_ko_B_list, ko_B_to_gene_member_dict, Bs_description_dict, stats_file_B_GeneNumber)
    write_out_stats_GeneNumber(identified_ko_C_list, ko_C_to_gene_member_dict, Cs_description_dict, stats_file_C_GeneNumber)
    write_out_stats_GeneNumber(identified_ko_D_list, ko_D_to_gene_member_dict, Ds_description_dict, stats_file_D_GeneNumber)
    if depth_file is not None:
        write_out_stats_TotalDepth(identified_ko_A_list, ko_A_to_gene_member_dict, gene_depth_dict, As_description_dict, stats_file_A_TotalDepth)
        write_out_stats_TotalDepth(identified_ko_B_list, ko_B_to_gene_member_dict, gene_depth_dict, Bs_description_dict, stats_file_B_TotalDepth)
        write_out_stats_TotalDepth(identified_ko_C_list, ko_C_to_gene_member_dict, gene_depth_dict, Cs_description_dict, stats_file_C_TotalDepth)
        write_out_stats_TotalDepth(identified_ko_D_list, ko_D_to_gene_member_dict, gene_depth_dict, Ds_description_dict, stats_file_D_TotalDepth)

    #################### write out GeneNumber and TotalDepth stats (pct) ####################

    AnnotateNorm(stats_file_A_GeneNumber, True, 2, len(genes_with_ko), stats_file_A_GeneNumber_pct, 'KO\tGeneNumber_pct\tDescription\n')
    AnnotateNorm(stats_file_B_GeneNumber, True, 2, len(genes_with_ko), stats_file_B_GeneNumber_pct, 'KO\tGeneNumber_pct\tDescription\n')
    AnnotateNorm(stats_file_C_GeneNumber, True, 2, len(genes_with_ko), stats_file_C_GeneNumber_pct, 'KO\tGeneNumber_pct\tDescription\n')
    AnnotateNorm(stats_file_D_GeneNumber, True, 2, len(genes_with_ko), stats_file_D_GeneNumber_pct, 'KO\tGeneNumber_pct\tDescription\n')
    if depth_file is not None:
        AnnotateNorm(stats_file_A_TotalDepth, True, 2, genes_with_ko_TotalDepth, stats_file_A_TotalDepth_pct, 'KO\tTotalDepth_pct\tDescription\n')
        AnnotateNorm(stats_file_B_TotalDepth, True, 2, genes_with_ko_TotalDepth, stats_file_B_TotalDepth_pct, 'KO\tTotalDepth_pct\tDescription\n')
        AnnotateNorm(stats_file_C_TotalDepth, True, 2, genes_with_ko_TotalDepth, stats_file_C_TotalDepth_pct, 'KO\tTotalDepth_pct\tDescription\n')
        AnnotateNorm(stats_file_D_TotalDepth, True, 2, genes_with_ko_TotalDepth, stats_file_D_TotalDepth_pct, 'KO\tTotalDepth_pct\tDescription\n')

    #################### write out GeneNumber and TotalDepth stats (pct_by_all) ####################

    if pct_by_all is True:
        AnnotateNorm(stats_file_A_GeneNumber, True, 2, len(query_seq_id_all), stats_file_A_GeneNumber_pct_by_all, 'KO\tGeneNumber_pct_by_all\tDescription\n')
        AnnotateNorm(stats_file_B_GeneNumber, True, 2, len(query_seq_id_all), stats_file_B_GeneNumber_pct_by_all, 'KO\tGeneNumber_pct_by_all\tDescription\n')
        AnnotateNorm(stats_file_C_GeneNumber, True, 2, len(query_seq_id_all), stats_file_C_GeneNumber_pct_by_all, 'KO\tGeneNumber_pct_by_all\tDescription\n')
        AnnotateNorm(stats_file_D_GeneNumber, True, 2, len(query_seq_id_all), stats_file_D_GeneNumber_pct_by_all, 'KO\tGeneNumber_pct_by_all\tDescription\n')
        if depth_file is not None:
            AnnotateNorm(stats_file_A_TotalDepth, True, 2, total_depth_for_all_query_genes, stats_file_A_TotalDepth_pct_by_all, 'KO\tTotalDepth_pct_by_all\tDescription\n')
            AnnotateNorm(stats_file_B_TotalDepth, True, 2, total_depth_for_all_query_genes, stats_file_B_TotalDepth_pct_by_all, 'KO\tTotalDepth_pct_by_all\tDescription\n')
            AnnotateNorm(stats_file_C_TotalDepth, True, 2, total_depth_for_all_query_genes, stats_file_C_TotalDepth_pct_by_all, 'KO\tTotalDepth_pct_by_all\tDescription\n')
            AnnotateNorm(stats_file_D_TotalDepth, True, 2, total_depth_for_all_query_genes, stats_file_D_TotalDepth_pct_by_all, 'KO\tTotalDepth_pct_by_all\tDescription\n')


def get_KEGG_annot_df(annotation_dir, stats_level, annotation_df_absolute_num, annotation_df_pct, annotation_df_pct_by_all, ABCD_description_dict, with_depth, pct_by_all, include_ko_fun):

    annotation_dir_re = '%s/*_KEGG_wd' % annotation_dir
    annotation_folder_list = [os.path.basename(file_name) for file_name in glob.glob(annotation_dir_re)]

    ko_num_dict = {}
    ko_num_pct_dict = {}
    ko_num_pct_by_all_dict = {}
    all_identified_ko = set()
    for annotation_folder in annotation_folder_list:

        annotation_folder_basename = annotation_folder.split('_KEGG_wd')[0]

        if with_depth is False:
            pwd_annotation_stats_file =             '%s/%s/%s_ko_stats_%s.txt'                  % (annotation_dir, annotation_folder, annotation_folder_basename, stats_level)
            pwd_annotation_stats_file_pct =         '%s/%s/%s_ko_stats_%s_pct.txt'              % (annotation_dir, annotation_folder, annotation_folder_basename, stats_level)
            pwd_annotation_stats_file_pct_by_all =  '%s/%s/%s_ko_stats_%s_pct_by_all.txt'       % (annotation_dir, annotation_folder, annotation_folder_basename, stats_level)
        else:
            pwd_annotation_stats_file =             '%s/%s/%s_ko_stats_%s_depth.txt'            % (annotation_dir, annotation_folder, annotation_folder_basename, stats_level)
            pwd_annotation_stats_file_pct =         '%s/%s/%s_ko_stats_%s_depth_pct.txt'        % (annotation_dir, annotation_folder, annotation_folder_basename, stats_level)
            pwd_annotation_stats_file_pct_by_all =  '%s/%s/%s_ko_stats_%s_depth_pct_by_all.txt' % (annotation_dir, annotation_folder, annotation_folder_basename, stats_level)

        current_ko_to_num_dict = {}
        for ko in open(pwd_annotation_stats_file):
            if not ko.startswith('KO\t'):
                ko_split = ko.strip().split('\t')
                if with_depth is False:
                    current_ko_to_num_dict[ko_split[0]] = int(ko_split[1])
                else:
                    current_ko_to_num_dict[ko_split[0]] = float(ko_split[1])
                all_identified_ko.add(ko_split[0])

        current_ko_to_num_pct_dict = {}
        for ko in open(pwd_annotation_stats_file_pct):
            if not ko.startswith('KO\t'):
                ko_split = ko.strip().split('\t')
                current_ko_to_num_pct_dict[ko_split[0]] = float(ko_split[1])
                all_identified_ko.add(ko_split[0])

        if pct_by_all is True:
            current_ko_to_num_pct_by_all_dict = {}
            for ko in open(pwd_annotation_stats_file_pct_by_all):
                if not ko.startswith('KO\t'):
                    ko_split = ko.strip().split('\t')
                    current_ko_to_num_pct_by_all_dict[ko_split[0]] = float(ko_split[1])
                    all_identified_ko.add(ko_split[0])

        ko_num_dict[annotation_folder_basename] = current_ko_to_num_dict
        ko_num_pct_dict[annotation_folder_basename] = current_ko_to_num_pct_dict
        if pct_by_all is True:
            ko_num_pct_by_all_dict[annotation_folder_basename] = current_ko_to_num_pct_by_all_dict

    all_identified_ko_list = sorted([i for i in all_identified_ko])

    all_identified_ko_list_with_description = []
    for each_ko in all_identified_ko_list:
        ko_with_desc = '%s_%s' % (each_ko, ABCD_description_dict[each_ko])
        ko_with_desc = ko_with_desc.replace(' ', '_')
        all_identified_ko_list_with_description.append(ko_with_desc)

    list_to_use = all_identified_ko_list
    if include_ko_fun is True:
        list_to_use = all_identified_ko_list_with_description

    annotation_df_absolute_num_handle = open(annotation_df_absolute_num, 'w')
    annotation_df_absolute_num_handle.write('\t%s\n' % '\t'.join(list_to_use))
    annotation_df_percentage_handle = open(annotation_df_pct, 'w')
    annotation_df_percentage_handle.write('\t%s\n' % '\t'.join(list_to_use))
    if pct_by_all is True:
        annotation_df_percentage_by_all_handle = open(annotation_df_pct_by_all, 'w')
        annotation_df_percentage_by_all_handle.write('\t%s\n' % '\t'.join(list_to_use))
    for annotation_folder in sorted(annotation_folder_list):

        annotation_folder_basename = annotation_folder.split('_KEGG_wd')[0]
        current_ko_num_dict = ko_num_dict[annotation_folder_basename]
        current_ko_num_dict_pct = ko_num_pct_dict[annotation_folder_basename]
        if pct_by_all is True:
            current_ko_num_dict_pct_by_all = ko_num_pct_by_all_dict[annotation_folder_basename]

        current_ko_num_list = []
        current_ko_num_list_pct = []
        current_ko_num_list_pct_by_all = []
        for identified_ko in all_identified_ko_list:

            # get num list
            identified_ko_num = 0
            identified_ko_num_pct = 0
            identified_ko_num_pct_by_all = 0
            if identified_ko in current_ko_num_dict:
                identified_ko_num = current_ko_num_dict[identified_ko]
                identified_ko_num_pct = current_ko_num_dict_pct[identified_ko]
                if pct_by_all is True:
                    identified_ko_num_pct_by_all = current_ko_num_dict_pct_by_all[identified_ko]

            current_ko_num_list.append(identified_ko_num)
            current_ko_num_list_pct.append(identified_ko_num_pct)
            if pct_by_all is True:
                current_ko_num_list_pct_by_all.append(identified_ko_num_pct_by_all)

        # write out
        annotation_df_absolute_num_handle.write('%s\t%s\n' % (annotation_folder_basename, '\t'.join([str(i) for i in current_ko_num_list])))
        annotation_df_percentage_handle.write('%s\t%s\n' % (annotation_folder_basename, '\t'.join([str(i) for i in current_ko_num_list_pct])))
        if pct_by_all is True:
            annotation_df_percentage_by_all_handle.write('%s\t%s\n' % (annotation_folder_basename, '\t'.join([str(i) for i in current_ko_num_list_pct_by_all])))

    annotation_df_absolute_num_handle.close()
    annotation_df_percentage_handle.close()
    if pct_by_all is True:
        annotation_df_percentage_by_all_handle.close()


def Annotation_KEGG(args):

    input_file_faa      = args['seq_in']
    input_file_user_ko  = args['ko_in']
    file_extension      = args['x']
    depth_file          = args['d']
    pct_by_all          = args['pct_by_all']
    KEGG_DB_folder      = args['db_dir']
    run_diamond         = args['diamond']
    num_threads         = args['t']
    evalue_cutoff       = args['evalue']
    include_desc        = args['desc']

    run_blast = None
    if (input_file_faa is not None) and (input_file_user_ko is None):
        run_blast = True
    elif (input_file_faa is None) and (input_file_user_ko is not None):
        run_blast = False
    else:
        print(datetime.now().strftime(time_format) + 'Please provide input file with either "-seq_in" or "-ko_in", do not provide both')
        exit()

    if run_blast is True:
        input_file_folder = input_file_faa
    else:
        input_file_folder = input_file_user_ko

    # check whether input file/folder exist
    if (os.path.isfile(input_file_folder) is False) and (os.path.isdir(input_file_folder) is False):
        print(datetime.now().strftime(time_format) + 'input file/folder not found, program exited')
        exit()

    if run_blast is True:
        print(datetime.now().strftime(time_format) + 'Input sequence file detected, will run blastp/diamond first')
        sleep(0.5)
    else:
        print(datetime.now().strftime(time_format) + 'Annotation results provided, blastp/diamond skipped')
        sleep(0.5)

    ################################################# define file name #################################################

    KEGG_DB_seq         = '%s/kegg_db_seq.fasta'        % KEGG_DB_folder
    KEGG_DB_seq_diamond = '%s/kegg_db_seq.fasta.dmnd'   % KEGG_DB_folder
    KEGG_DB_seq2ko      = '%s/kegg_db_seq2ko.txt'       % KEGG_DB_folder
    KEGG_DB_ko          = '%s/ko00001.keg'              % KEGG_DB_folder

    ########################################## check whether diamond db exist ##########################################

    if (run_blast is True) and (run_diamond is True):
        if os.path.isfile(KEGG_DB_seq_diamond) is False:
            print(datetime.now().strftime(time_format) + 'DB file not found, making diamond db with %s' % KEGG_DB_seq)

            if os.path.isfile(KEGG_DB_seq) is True:
                diamond_makedb_cmd = 'diamond makedb --in %s --db %s --quiet' % (KEGG_DB_seq, KEGG_DB_seq_diamond)
                os.system(diamond_makedb_cmd)
            else:
                print(datetime.now().strftime(time_format) + '%s not found, program exited' % KEGG_DB_seq)
                exit()

    ########################################### check whether blast+ db exist ##########################################

    if (run_blast is True) and (run_diamond is False):

        unfound_db_index_file = []
        for db_index in ['phr', 'pin', 'pnd', 'pni', 'pog', 'psd', 'psi', 'psq']:
            pwd_db_index = '%s/kegg_db_seq.fasta.%s' % (KEGG_DB_folder, db_index)
            if not os.path.isfile(pwd_db_index):
                unfound_db_index_file.append(db_index)
        if len(unfound_db_index_file) > 0:
            print(datetime.now().strftime(time_format) + 'blast db index not found, runing makeblastdb first')
            makeblastdb_cmd = 'makeblastdb -in %s -dbtype prot -parse_seqids -logfile %s.log' % (KEGG_DB_seq, KEGG_DB_seq)
            os.system(makeblastdb_cmd)
            print(datetime.now().strftime(time_format) + 'makeblastdb finished')

    ######################################### Run blastp with multiprocessing ##########################################

    # check whether the input file is a file or folder
    if os.path.isfile(input_file_folder) is True:
        input_file_path, input_file_basename, input_file_ext = sep_path_basename_ext(input_file_folder)
        run_blast_worker([input_file_folder, run_blast, run_diamond, KEGG_DB_seq, KEGG_DB_seq_diamond, input_file_path, evalue_cutoff, num_threads])

    if os.path.isdir(input_file_folder) is True:

        # create output folder
        output_folder = '%s_KEGG_wd' % input_file_folder
        force_create_folder(output_folder)

        # check whether input genome exist
        input_file_re = '%s/*.%s' % (input_file_folder, file_extension)
        input_file_name_list = [os.path.basename(file_name) for file_name in glob.glob(input_file_re)]

        if len(input_file_name_list) == 0:
            print(datetime.now().strftime(time_format) + 'input file not found, program exited')
            exit()

        # run blastp with multiprocessing
        if run_blast is True:
            print(datetime.now().strftime(time_format) + 'Running Blast/Diamond for %s input files with %s cores' % (len(input_file_name_list), num_threads))

        list_for_multiple_arguments_blast = []
        for input_file in input_file_name_list:
            pwd_input_file = '%s/%s' % (input_file_folder, input_file)
            list_for_multiple_arguments_blast.append([pwd_input_file, run_blast, run_diamond, KEGG_DB_seq, KEGG_DB_seq_diamond, output_folder, evalue_cutoff, 1])

        # run blastp with multiprocessing
        pool = mp.Pool(processes=num_threads)
        pool.map(run_blast_worker, list_for_multiple_arguments_blast)
        pool.close()
        pool.join()

    ############################################## Read in KEGG DB files ###############################################

    print(datetime.now().strftime(time_format) + 'Read in KEGG DB files')

    As_description_dict = {}
    Bs_description_dict = {}
    Cs_description_dict = {}
    Ds_description_dict = {}
    D2ABCD_dict = {}
    current_A = ''
    current_B = ''
    current_C = ''
    for each_line in open(KEGG_DB_ko):
        if each_line[0] in ['A', 'B', 'C', 'D']:
            each_line_split = each_line.strip().split(' ')

            if each_line[0] == 'A':
                current_A_id = each_line_split[0]
                current_A_description = ' '.join(each_line_split[1:])
                current_A = current_A_id
                As_description_dict[current_A_id] = current_A_description

            elif each_line[0] == 'B':
                if len(each_line_split) > 1:
                    current_B_id = each_line_split[2]
                    current_B_description = ' '.join(each_line_split[3:])
                    current_B = current_B_id
                    Bs_description_dict[current_B_id] = current_B_description

            elif each_line[0] == 'C':
                current_C_id = each_line_split[4]
                current_C_description = ' '.join(each_line_split[5:])
                current_C = current_C_id
                Cs_description_dict[current_C_id] = current_C_description

            elif each_line[0] == 'D':
                current_D_id = each_line_split[6]
                current_D_description = ' '.join(each_line_split[7:])
                Ds_description_dict[current_D_id] = current_D_description
                ABCD_value = 'A_%s|B_%s|C_%s|D_%s' % (current_A, current_B, current_C, current_D_id)
                if current_D_id not in D2ABCD_dict:
                    D2ABCD_dict[current_D_id] = [ABCD_value]
                elif (current_D_id in D2ABCD_dict) and (ABCD_value not in D2ABCD_dict[current_D_id]):
                    D2ABCD_dict[current_D_id].append(ABCD_value)

    ABCD_description_dict = {}
    for each_a in As_description_dict:
        ABCD_description_dict[each_a] = As_description_dict[each_a]
    for each_b in Bs_description_dict:
        ABCD_description_dict[each_b] = Bs_description_dict[each_b]
    for each_c in Cs_description_dict:
        ABCD_description_dict[each_c] = Cs_description_dict[each_c]
    for each_d in Ds_description_dict:
        ABCD_description_dict[each_d] = Ds_description_dict[each_d]

    # get db_seq_to_KO_dict
    db_seq_to_KO_dict = {}
    if run_blast is True:
        for each_hit in open(KEGG_DB_seq2ko):
            each_hit_split = each_hit.strip().split('\t')
            db_seq = each_hit_split[0]
            hit_id_KO = each_hit_split[1]
            if hit_id_KO != '':
                db_seq_to_KO_dict[db_seq] = hit_id_KO

    ########################################################################################################################

    # check whether the input file is a file or folder
    if os.path.isfile(input_file_folder) is True:

        # check whether depth file exist
        if depth_file is not None:
            if os.path.isfile(depth_file) is False:
                print(datetime.now().strftime(time_format) + 'specified depth file not found, program exited!')
                exit()

        print(datetime.now().strftime(time_format) + 'Running KEGG annotation for 1 file with %s cores' % (num_threads))
        input_file_path, input_file_basename, input_file_ext = sep_path_basename_ext(input_file_folder)
        parse_blast_op_worker([input_file_folder, run_blast, As_description_dict, Bs_description_dict, Cs_description_dict, Ds_description_dict, D2ABCD_dict, db_seq_to_KO_dict, input_file_path, depth_file, pct_by_all])


    if os.path.isdir(input_file_folder) is True:

        input_file_re = '%s/*.%s' % (input_file_folder, file_extension)
        input_file_name_list = [os.path.basename(file_name) for file_name in glob.glob(input_file_re)]

        # check whether depth file exist
        if depth_file is not None:

            if os.path.isfile(depth_file) is True:
                print(datetime.now().strftime(
                    time_format) + 'please provide the folder containing individual depth files (with extension .depth) for each of your input sequence file.')
                print(datetime.now().strftime(time_format) + 'single depth file (not folder) detected, program exited!')
                exit()

            if os.path.isdir(depth_file) is False:
                print(datetime.now().strftime(time_format) + 'specified depth folder not found, program exited!')
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

        # create output folder
        output_folder = '%s_KEGG_wd' % input_file_folder
        input_folder_name = input_file_folder
        if '/' in input_file_folder:
            input_folder_name = input_file_folder.split('/')[-1]

        # parse blast results with multiprocessing
        if run_blast is True:
            print(datetime.now().strftime(time_format) + 'Parsing Blast/Diamond results for %s input files with %s cores' % (len(input_file_name_list), num_threads))

        list_for_multiple_arguments_parse_blast_op = []
        for input_file in input_file_name_list:

            input_file_basename = '.'.join(input_file.split('.')[:-1])
            pwd_input_file = '%s/%s' % (input_file_folder, input_file)

            # get path to current depth file
            if depth_file is None:
                input_file_depth = None
            else:
                input_file_depth = '%s/%s.depth' % (depth_file, input_file_basename)

            list_for_multiple_arguments_parse_blast_op.append([pwd_input_file, run_blast, As_description_dict, Bs_description_dict, Cs_description_dict, Ds_description_dict, D2ABCD_dict, db_seq_to_KO_dict, output_folder, input_file_depth, pct_by_all])

        # parse blast results with multiprocessing
        pool = mp.Pool(processes=num_threads)
        pool.map(parse_blast_op_worker, list_for_multiple_arguments_parse_blast_op)
        pool.close()
        pool.join()

        ######################################################### get dataframe #########################################################

        print(datetime.now().strftime(time_format) + 'Data matrix exported to:')

        for ko_level in ['A', 'B', 'C', 'D']:
            annotation_df_GeneNumber =            '%s/%s_%s.txt'                    % (output_folder, input_folder_name, ko_level)
            annotation_df_GeneNumber_pct =        '%s/%s_%s_pct.txt'                % (output_folder, input_folder_name, ko_level)
            annotation_df_GeneNumber_pct_by_all = '%s/%s_%s_pct_by_all.txt'         % (output_folder, input_folder_name, ko_level)
            annotation_df_TotalDepth =            '%s/%s_%s_depth.txt'              % (output_folder, input_folder_name, ko_level)
            annotation_df_TotalDepth_pct =        '%s/%s_%s_depth_pct.txt'          % (output_folder, input_folder_name, ko_level)
            annotation_df_TotalDepth_pct_by_all = '%s/%s_%s_depth_pct_by_all.txt'   % (output_folder, input_folder_name, ko_level)

            #################### get GeneNumber df and report ####################

            get_KEGG_annot_df(output_folder, ko_level, annotation_df_GeneNumber, annotation_df_GeneNumber_pct, annotation_df_GeneNumber_pct_by_all, ABCD_description_dict, with_depth=False, pct_by_all=pct_by_all, include_ko_fun=include_desc)

            print(annotation_df_GeneNumber.split('/')[-1])
            print(annotation_df_GeneNumber_pct.split('/')[-1])
            if pct_by_all is True:
                print(annotation_df_GeneNumber_pct_by_all.split('/')[-1])

            #################### get TotalDepth df and report ####################

            if depth_file is not None:
                get_KEGG_annot_df(output_folder, ko_level, annotation_df_TotalDepth, annotation_df_TotalDepth_pct, annotation_df_TotalDepth_pct_by_all, ABCD_description_dict, with_depth=True, pct_by_all=pct_by_all, include_ko_fun=include_desc)

                print(annotation_df_TotalDepth.split('/')[-1])
                print(annotation_df_TotalDepth_pct.split('/')[-1])
                if pct_by_all is True:
                    print(annotation_df_TotalDepth_pct_by_all.split('/')[-1])

    ################################################## Final report ####################################################

    print(datetime.now().strftime(time_format) + 'Done!')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(usage=KEGG_parser_usage)
    parser.add_argument('-seq_in',      required=False,                             help='faa file')
    parser.add_argument('-ko_in',       required=False,                             help='annotation results from BlastKOALA/GhostKOALA, normally with name user_ko.txt')
    parser.add_argument('-x',           required=False,                             help='file extension')
    parser.add_argument('-d',           required=False, default=None,               help='gene depth file/folder')
    parser.add_argument('-pct_by_all',  required=False, action='store_true',        help='normalize by all query genes, rather than those with ko assignment')
    parser.add_argument('-db_dir',      required=True,                              help='folder holds sequence, seq2ko and ko00001.keg files')
    parser.add_argument('-diamond',     required=False, action='store_true',        help='run diamond (for big dataset), default is NCBI blastp')
    parser.add_argument('-t',           required=False, default=1,     type=int,    help='number of threads, default: 1')
    parser.add_argument('-evalue',      required=False, default=0.0001, type=float, help='evalue cutoff, default: 0.0001')
    parser.add_argument('-desc',        required=False, action='store_true',        help='include KO functions in the final dataframe')
    args = vars(parser.parse_args())
    Annotation_KEGG(args)
