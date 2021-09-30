#!/usr/bin/env python3

import os
import glob
import argparse
from Bio import SeqIO
import multiprocessing as mp
from datetime import datetime
from Bio.SeqRecord import SeqRecord
from BioSAK.global_functions import time_format
from BioSAK.global_functions import force_create_folder
from BioSAK.global_functions import sep_path_basename_ext
from BioSAK.global_functions import AnnotateNorm
from BioSAK.global_functions import get_gene_list_TotalDepth


COG2020_parser_usage = '''
===================================== COG2020 example commands =====================================

# module needed
module load python/3.7.3
module load blast+/2.11.0
module load diamond/0.9.31

# annotate protein sequences
BioSAK COG2020 -m P -t 6 -db_dir /srv/scratch/z5039045/DB/COG2020 -i genes.faa 
BioSAK COG2020 -m P -t 6 -db_dir /srv/scratch/z5039045/DB/COG2020 -i faa_files -x faa -depth depth_files

# annotate DNA sequences (ORFs)
BioSAK COG2020 -m N -t 6 -db_dir /srv/scratch/z5039045/DB/COG2020 -i genes.ffn -depth gene.depth
BioSAK COG2020 -m N -t 6 -db_dir /srv/scratch/z5039045/DB/COG2020 -i ffn_files -x ffn

# Depth file format (one gene per line, tab separated)
gene_1	30
gene_2	10.58

# Prepare DB files (version 2020):
cd db_dir
wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.fa.gz
wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.cog.csv
wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab
wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab
wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/Readme.2020-11-25.txt
gunzip cog-20.fa.gz
module load blast+/2.11.0
makeblastdb -in cog-20.fa -dbtype prot -parse_seqids -logfile cog-20.fa.log
module load diamond/0.9.31
diamond makedb --in cog-20.fa --db cog-20.fa.dmnd --quiet

# How it works:
1. COG2020 module uses Blast+/Diamond to get the best hits of query genes in the database 
   with users defined e-value cutoff (default 0.001).
2. The TotalDepth of a COG id/category is obtained by summing up the depth of all genes assigned to it.
3. The percentage of GeneNumber/TotalDepth of genes assigned to a COG is calculated by dividing them 
   by the total number/depth of genes with COG assignment (default) or all query genes in a file (if "-pct_by_all" specified). 

# Note!!!
If you run COG2020 for multiple files in a batch manner and want to have their depth info incorporated into the results, 
you need to provide a folder containing individual depth file for each of your input sequence file.
Name of the depth file needs to be exactly the same as its corresponding sequence file, except the extension is ".depth".

====================================================================================================
'''


def dna2aa(dna_file, aa_file):
    query_aa_handle = open(aa_file, 'w')
    for each in SeqIO.parse(dna_file, 'fasta'):
        each_aa = each.seq.translate()
        each_aa_record = SeqRecord(each_aa)
        each_aa_record.id = each.id
        each_aa_record.description = each.description
        SeqIO.write(each_aa_record, query_aa_handle, 'fasta')
    query_aa_handle.close()


def best_hit(args):

    file_in = args['i']
    file_out = args['o']

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


def COG2020_worker(argument_list):

    pwd_input_file =                    argument_list[0]
    pwd_prot2003_2014 =                 argument_list[1]
    protein_id_to_cog_id_dict =         argument_list[2]
    cog_id_to_category_dict =           argument_list[3]
    cog_id_to_description_dict =        argument_list[4]
    cog_category_list =                 argument_list[5]
    cog_category_to_description_dict =  argument_list[6]
    sequence_type =                     argument_list[7]
    output_folder =                     argument_list[8]
    thread_num =                        argument_list[9]
    run_diamond =                       argument_list[10]
    evalue_cutoff =                     argument_list[11]
    depth_file =                        argument_list[12]
    pct_by_all =                        argument_list[13]

    input_seq_no_path, input_seq_no_ext, input_seq_ext = sep_path_basename_ext(pwd_input_file)
    current_output_folder = '%s/%s_COG2020_wd' % (output_folder, input_seq_no_ext)

    pwd_blastp_output =                        '%s/%s_blastp.tab'                             % (current_output_folder, input_seq_no_ext)
    pwd_blastp_output_besthits =               '%s/%s_blastp_besthits.tab'                    % (current_output_folder, input_seq_no_ext)
    pwd_query_to_cog_txt =                     '%s/%s_query_to_cog.txt'                       % (current_output_folder, input_seq_no_ext)

    pwd_cog_stats_GeneNumber =                 '%s/%s_cog_stats_GeneNumber.txt'               % (current_output_folder, input_seq_no_ext)
    pwd_cog_stats_TotalDepth =                 '%s/%s_cog_stats_TotalDepth.txt'               % (current_output_folder, input_seq_no_ext)
    pwd_func_stats_GeneNumber =                '%s/%s_func_stats_GeneNumber.txt'              % (current_output_folder, input_seq_no_ext)
    pwd_func_stats_TotalDepth =                '%s/%s_func_stats_TotalDepth.txt'              % (current_output_folder, input_seq_no_ext)

    pwd_cog_stats_GeneNumber_pct =             '%s/%s_cog_stats_GeneNumber_pct.txt'           % (current_output_folder, input_seq_no_ext)
    pwd_cog_stats_TotalDepth_pct =             '%s/%s_cog_stats_TotalDepth_pct.txt'           % (current_output_folder, input_seq_no_ext)
    pwd_func_stats_GeneNumber_pct =            '%s/%s_func_stats_GeneNumber_pct.txt'          % (current_output_folder, input_seq_no_ext)
    pwd_func_stats_TotalDepth_pct =            '%s/%s_func_stats_TotalDepth_pct.txt'          % (current_output_folder, input_seq_no_ext)

    pwd_cog_stats_GeneNumber_pct_by_all =      '%s/%s_cog_stats_GeneNumber_pct_by_all.txt'    % (current_output_folder, input_seq_no_ext)
    pwd_cog_stats_TotalDepth_pct_by_all =      '%s/%s_cog_stats_TotalDepth_pct_by_all.txt'    % (current_output_folder, input_seq_no_ext)
    pwd_func_stats_GeneNumber_pct_by_all =     '%s/%s_func_stats_GeneNumber_pct_by_all.txt'   % (current_output_folder, input_seq_no_ext)
    pwd_func_stats_TotalDepth_pct_by_all =     '%s/%s_func_stats_TotalDepth_pct_by_all.txt'   % (current_output_folder, input_seq_no_ext)

    force_create_folder(current_output_folder)

    input_seq_aa = ''
    if sequence_type in ['N', 'n']:
        input_seq_aa = '%s_aa.fasta' % input_seq_no_ext
        dna2aa(pwd_input_file, input_seq_aa)
    elif sequence_type in ['P', 'p']:
        input_seq_aa = pwd_input_file
    else:
        print('Specified input sequence type unrecognizable, program exited!')
        exit()

    # run blastp
    if run_diamond is False:
        os.system('blastp -query %s -db %s -out %s -evalue %s -outfmt 6 -show_gis -num_threads %s' % (input_seq_aa, pwd_prot2003_2014, pwd_blastp_output, evalue_cutoff, thread_num))
    else:
        os.system('diamond blastp -q %s --db %s.dmnd --out %s --evalue %s --outfmt 6 --threads %s --quiet' % (input_seq_aa, pwd_prot2003_2014, pwd_blastp_output, evalue_cutoff, thread_num))

    # keep only best hits
    best_hit({'i': pwd_blastp_output, 'o': pwd_blastp_output_besthits})

    # get query_to_ref_protein_dict
    query_to_ref_protein_dict = {}
    for each_hit in open(pwd_blastp_output_besthits):
        each_hit_split = each_hit.strip().split('\t')
        each_hit_query = each_hit_split[0]
        each_hit_subject = each_hit_split[1]
        each_hit_subject_no_dot = '_'.join(each_hit_subject.split('.'))
        query_to_ref_protein_dict[each_hit_query] = each_hit_subject_no_dot

    # get query sequences list
    query_seq_list = []
    for query_seq in SeqIO.parse(pwd_input_file, 'fasta'):
        query_seq_list.append(query_seq.id)

    # export annotation
    cog_id_num_dict = {}
    cog_id_to_gene_member_dict = {}
    cog_cate_num_dict = {}
    cog_cate_to_gene_member_dict = {}
    genes_with_cog = set()
    pwd_query_to_cog_txt_handle = open(pwd_query_to_cog_txt, 'w')
    pwd_query_to_cog_txt_handle.write('Query\tCOG\tCategory\tDescription\n')
    for query_gene in sorted(query_seq_list):

        if query_gene not in query_to_ref_protein_dict:
            pwd_query_to_cog_txt_handle.write('%s\n' % (query_gene))

        else:
            db_protein_id = query_to_ref_protein_dict[query_gene]
            if db_protein_id not in protein_id_to_cog_id_dict:
                pwd_query_to_cog_txt_handle.write('%s\n' % (query_gene))

            else:
                cog_id_list = protein_id_to_cog_id_dict[db_protein_id]
                for cog_id in cog_id_list:
                    cog_cate = cog_id_to_category_dict[cog_id]
                    cog_des = cog_id_to_description_dict[cog_id]
                    pwd_query_to_cog_txt_handle.write('%s\t%s\t%s\t%s\n' % (query_gene, cog_id, cog_cate, cog_des))
                    genes_with_cog.add(query_gene)

                    # update cog_id_num_dict
                    if cog_id not in cog_id_num_dict:
                        cog_id_num_dict[cog_id] = 1
                        cog_id_to_gene_member_dict[cog_id] = [query_gene]
                    else:
                        cog_id_num_dict[cog_id] += 1
                        cog_id_to_gene_member_dict[cog_id].append(query_gene)

                    # update cog_cate_num_dict
                    for each_cog_cate in cog_cate:
                        if each_cog_cate not in cog_cate_num_dict:
                            cog_cate_num_dict[each_cog_cate] = 1
                            cog_cate_to_gene_member_dict[each_cog_cate] = [query_gene]
                        else:
                            cog_cate_num_dict[each_cog_cate] += 1
                            cog_cate_to_gene_member_dict[each_cog_cate].append(query_gene)

    pwd_query_to_cog_txt_handle.close()


    # read in depth info
    gene_depth_dict = {}
    if depth_file is not None:
        for each_depth in open(depth_file):
            each_depth_split = each_depth.strip().split('\t')
            gene_depth_dict[each_depth_split[0]] = float(each_depth_split[1])

    # get TotalDepth of all query genes or genes with cog assignment
    if depth_file is not None:
        genes_with_cog_TotalDepth = get_gene_list_TotalDepth(genes_with_cog, gene_depth_dict)
        total_depth_for_all_query_genes = get_gene_list_TotalDepth(query_seq_list, gene_depth_dict)

    #################### export cog_stats_GeneNumber ####################

    pwd_cog_stats_GeneNumber_handle = open(pwd_cog_stats_GeneNumber, 'w')
    pwd_cog_stats_GeneNumber_handle.write('COG\tGeneNumber\tDescription\n')
    for each_cog_id in cog_id_num_dict:
        each_cog_id_GeneNumber = cog_id_num_dict[each_cog_id]
        pwd_cog_stats_GeneNumber_handle.write('%s\t%s\t%s\n' % (each_cog_id, each_cog_id_GeneNumber, cog_id_to_description_dict[each_cog_id]))
    pwd_cog_stats_GeneNumber_handle.close()

    #################### export cog_stats_TotalDepth ####################

    if depth_file is not None:
        pwd_cog_stats_TotalDepth_handle = open(pwd_cog_stats_TotalDepth, 'w')
        pwd_cog_stats_TotalDepth_handle.write('COG\tTotalDepth\tDescription\n')
        for each_cog_id in cog_id_to_gene_member_dict:
            each_cog_id_gene_member = cog_id_to_gene_member_dict[each_cog_id]
            each_cog_id_TotalDepth = 0
            for each_gene in each_cog_id_gene_member:
                each_gene_depth = gene_depth_dict[each_gene]
                each_cog_id_TotalDepth += each_gene_depth
            each_cog_id_TotalDepth = float("{0:.2f}".format(each_cog_id_TotalDepth))
            pwd_cog_stats_TotalDepth_handle.write('%s\t%s\t%s\n' % (each_cog_id, each_cog_id_TotalDepth, cog_id_to_description_dict[each_cog_id]))
        pwd_cog_stats_TotalDepth_handle.close()

    #################### export func_stats_GeneNumber ####################

    pwd_func_stats_GeneNumber_handle = open(pwd_func_stats_GeneNumber, 'w')
    pwd_func_stats_GeneNumber_handle.write('Category\tGeneNumber\tDescription\n')
    for each_cog_cate in cog_category_list:
        each_cog_cate_GeneNumber = 0
        if each_cog_cate in cog_cate_num_dict:
            each_cog_cate_GeneNumber = cog_cate_num_dict[each_cog_cate]
        pwd_func_stats_GeneNumber_handle.write('%s\t%s\t%s\n' % (each_cog_cate, each_cog_cate_GeneNumber, cog_category_to_description_dict[each_cog_cate]))
    pwd_func_stats_GeneNumber_handle.close()

    #################### export func_stats_TotalDepth ####################

    if depth_file is not None:
        pwd_func_stats_TotalDepth_handle = open(pwd_func_stats_TotalDepth, 'w')
        pwd_func_stats_TotalDepth_handle.write('Category\tTotalDepth\tDescription\n')
        for each_cog_cate in cog_category_list:
            each_cog_cate_TotalDepth = 0
            if each_cog_cate in cog_cate_to_gene_member_dict:
                each_cog_cate_gene_member = cog_cate_to_gene_member_dict[each_cog_cate]
                for each_gene in each_cog_cate_gene_member:
                    each_gene_depth = gene_depth_dict[each_gene]
                    each_cog_cate_TotalDepth += each_gene_depth
            each_cog_cate_TotalDepth = float("{0:.2f}".format(each_cog_cate_TotalDepth))
            pwd_func_stats_TotalDepth_handle.write('%s\t%s\t%s\n' % (each_cog_cate, each_cog_cate_TotalDepth, cog_category_to_description_dict[each_cog_cate]))
        pwd_func_stats_TotalDepth_handle.close()

    #################### get pct files ####################

    AnnotateNorm(file_in=pwd_cog_stats_GeneNumber,  skip_header=True, value_column=2, Divisor_value=len(genes_with_cog), file_out=pwd_cog_stats_GeneNumber_pct,  file_out_header='Category\tGeneNumber_pct\tDescription\n')
    AnnotateNorm(file_in=pwd_func_stats_GeneNumber, skip_header=True, value_column=2, Divisor_value=len(genes_with_cog), file_out=pwd_func_stats_GeneNumber_pct, file_out_header='Category\tGeneNumber_pct\tDescription\n')
    if depth_file is not None:
        AnnotateNorm(file_in=pwd_cog_stats_TotalDepth,  skip_header=True, value_column=2, Divisor_value=genes_with_cog_TotalDepth, file_out=pwd_cog_stats_TotalDepth_pct,  file_out_header='Category\tTotalDepth_pct\tDescription\n')
        AnnotateNorm(file_in=pwd_func_stats_TotalDepth, skip_header=True, value_column=2, Divisor_value=genes_with_cog_TotalDepth, file_out=pwd_func_stats_TotalDepth_pct, file_out_header='Category\tTotalDepth_pct\tDescription\n')
    if pct_by_all is True:
        AnnotateNorm(file_in=pwd_cog_stats_GeneNumber,  skip_header=True, value_column=2, Divisor_value=len(query_seq_list), file_out=pwd_cog_stats_GeneNumber_pct_by_all,  file_out_header='Category\tGeneNumber_pct_by_all\tDescription\n')
        AnnotateNorm(file_in=pwd_func_stats_GeneNumber, skip_header=True, value_column=2, Divisor_value=len(query_seq_list), file_out=pwd_func_stats_GeneNumber_pct_by_all, file_out_header='Category\tGeneNumber_pct_by_all\tDescription\n')
        if depth_file is not None:
            AnnotateNorm(file_in=pwd_cog_stats_TotalDepth,  skip_header=True, value_column=2, Divisor_value=total_depth_for_all_query_genes, file_out=pwd_cog_stats_TotalDepth_pct_by_all,  file_out_header='Category\tTotalDepth_pct_by_all\tDescription\n')
            AnnotateNorm(file_in=pwd_func_stats_TotalDepth, skip_header=True, value_column=2, Divisor_value=total_depth_for_all_query_genes, file_out=pwd_func_stats_TotalDepth_pct_by_all, file_out_header='Category\tTotalDepth_pct_by_all\tDescription\n')


def get_COG_annot_df(annotation_dir, stats_level, annotation_df_absolute_num, annotation_df_percentage, annotation_df_percentage_by_all, with_depth, pct_by_all):

    annotation_dir_re = '%s/*_COG2020_wd' % annotation_dir
    annotation_folder_list = [os.path.basename(file_name) for file_name in glob.glob(annotation_dir_re)]

    cog_num_dict = {}
    cog_num_pct_dict = {}
    cog_num_pct_by_all_dict = {}
    all_identified_cog = set()
    for annotation_folder in annotation_folder_list:

        annotation_folder_basename = annotation_folder.split('_COG2020_wd')[0]

        pwd_annotation_stats_file = ''
        pwd_annotation_stats_file_pct = ''
        if stats_level == 'cog_id':
            if with_depth is False:
                pwd_annotation_stats_file =             '%s/%s/%s_cog_stats_GeneNumber.txt'                 % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct =         '%s/%s/%s_cog_stats_GeneNumber_pct.txt'             % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct_by_all =  '%s/%s/%s_cog_stats_GeneNumber_pct_by_all.txt'      % (annotation_dir, annotation_folder, annotation_folder_basename)
            else:
                pwd_annotation_stats_file =             '%s/%s/%s_cog_stats_TotalDepth.txt'                 % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct =         '%s/%s/%s_cog_stats_TotalDepth_pct.txt'             % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct_by_all =  '%s/%s/%s_cog_stats_TotalDepth_pct_by_all.txt'      % (annotation_dir, annotation_folder, annotation_folder_basename)

        if stats_level == 'cog_cate':
            if with_depth is False:
                pwd_annotation_stats_file =             '%s/%s/%s_func_stats_GeneNumber.txt'                % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct =         '%s/%s/%s_func_stats_GeneNumber_pct.txt'            % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct_by_all =  '%s/%s/%s_func_stats_GeneNumber_pct_by_all.txt'     % (annotation_dir, annotation_folder, annotation_folder_basename)
            else:
                pwd_annotation_stats_file =             '%s/%s/%s_func_stats_TotalDepth.txt'                % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct =         '%s/%s/%s_func_stats_TotalDepth_pct.txt'            % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct_by_all =  '%s/%s/%s_func_stats_TotalDepth_pct_by_all.txt'     % (annotation_dir, annotation_folder, annotation_folder_basename)

        current_cog_to_num_dict = {}
        for cog in open(pwd_annotation_stats_file):
            if (not cog.startswith('Category	')) and (not cog.startswith('COG	')):
                cog_split = cog.strip().split('\t')
                if with_depth is False:
                    current_cog_to_num_dict[cog_split[0]] = int(cog_split[1])
                else:
                    current_cog_to_num_dict[cog_split[0]] = float(cog_split[1])
                all_identified_cog.add(cog_split[0])
        cog_num_dict[annotation_folder_basename] = current_cog_to_num_dict


        current_cog_to_num_pct_dict = {}
        for cog in open(pwd_annotation_stats_file_pct):
            if (not cog.startswith('Category	')) and (not cog.startswith('COG	')):
                cog_split = cog.strip().split('\t')
                current_cog_to_num_pct_dict[cog_split[0]] = float(cog_split[1])
                all_identified_cog.add(cog_split[0])
        cog_num_pct_dict[annotation_folder_basename] = current_cog_to_num_pct_dict


        if pct_by_all is True:
            current_cog_to_num_pct_by_all_dict = {}
            for cog in open(pwd_annotation_stats_file_pct_by_all):
                if (not cog.startswith('Category	')) and (not cog.startswith('COG	')):
                    cog_split = cog.strip().split('\t')
                    current_cog_to_num_pct_by_all_dict[cog_split[0]] = float(cog_split[1])
                    all_identified_cog.add(cog_split[0])
            cog_num_pct_by_all_dict[annotation_folder_basename] = current_cog_to_num_pct_by_all_dict


    all_identified_cog_list = sorted([i for i in all_identified_cog])


    annotation_df_absolute_num_handle = open(annotation_df_absolute_num, 'w')
    annotation_df_percentage_handle = open(annotation_df_percentage, 'w')
    annotation_df_absolute_num_handle.write('\t%s\n' % '\t'.join(all_identified_cog_list))
    annotation_df_percentage_handle.write('\t%s\n' % '\t'.join(all_identified_cog_list))
    for annotation_folder in sorted(annotation_folder_list):

        annotation_folder_basename = annotation_folder.split('_COG2020_wd')[0]
        current_cog_num_dict = cog_num_dict[annotation_folder_basename]
        current_cog_num_pct_dict = cog_num_pct_dict[annotation_folder_basename]

        current_cog_num_list = []
        current_cog_num_list_pct = []
        for identified_cog in all_identified_cog_list:

            # get num list
            identified_cog_num = 0
            identified_cog_num_pct = 0
            if identified_cog in current_cog_num_dict:
                identified_cog_num = current_cog_num_dict[identified_cog]
                identified_cog_num_pct = current_cog_num_pct_dict[identified_cog]

            current_cog_num_list.append(identified_cog_num)
            current_cog_num_list_pct.append(identified_cog_num_pct)

        # write out
        annotation_df_absolute_num_handle.write('%s\t%s\n' % (annotation_folder_basename, '\t'.join([str(i) for i in current_cog_num_list])))
        annotation_df_percentage_handle.write('%s\t%s\n' % (annotation_folder_basename, '\t'.join([str(i) for i in current_cog_num_list_pct])))

    annotation_df_absolute_num_handle.close()
    annotation_df_percentage_handle.close()

    ########## write out pct_by_all_file ##########

    if pct_by_all is True:
        annotation_df_percentage_by_all_handle = open(annotation_df_percentage_by_all, 'w')
        annotation_df_percentage_by_all_handle.write('\t%s\n' % '\t'.join(all_identified_cog_list))
        for annotation_folder in sorted(annotation_folder_list):

            annotation_folder_basename = annotation_folder.split('_COG2020_wd')[0]
            current_cog_num_pct_by_all_dict = cog_num_pct_by_all_dict[annotation_folder_basename]

            current_cog_num_list_pct_by_all = []
            for identified_cog in all_identified_cog_list:

                # get num list
                identified_cog_num_pct_by_all = 0
                if identified_cog in current_cog_num_pct_by_all_dict:
                    identified_cog_num_pct_by_all = current_cog_num_pct_by_all_dict[identified_cog]

                current_cog_num_list_pct_by_all.append(identified_cog_num_pct_by_all)

            # write out
            annotation_df_percentage_by_all_handle.write('%s\t%s\n' % (annotation_folder_basename, '\t'.join([str(i) for i in current_cog_num_list_pct_by_all])))

        annotation_df_percentage_by_all_handle.close()


def COG2020(args):

    file_in =           args['i']
    file_extension =    args['x']
    sequence_type =     args['m']
    depth_file =        args['depth']
    pct_by_all =        args['pct_by_all']
    DB_dir =            args['db_dir']
    num_threads =       args['t']
    run_diamond =       args['diamond']
    evalue_cutoff =     args['evalue']

    pwd_cog_20_fa           = '%s/cog-20.fa'         % DB_dir
    pwd_cog_20_fa_diamond   = '%s/cog-20.fa.dmnd'    % DB_dir
    pwd_cog_20_cog_csv      = '%s/cog-20.cog.csv'    % DB_dir
    pwd_cog_20_def_tab      = '%s/cog-20.def.tab'    % DB_dir
    pwd_fun_20_tab          = '%s/fun-20.tab'        % DB_dir


    ############################################ check whether db file exist ###########################################

    # check whether db file exist
    unfound_inputs = []
    for each_input in [pwd_cog_20_fa, pwd_cog_20_def_tab, pwd_fun_20_tab]:
        if (not os.path.isfile(each_input)) and (not os.path.isdir(each_input)):
            unfound_inputs.append(each_input)
    if len(unfound_inputs) > 0:
        for each_unfound in unfound_inputs:
            print('%s not found' % each_unfound)
        exit()

    if run_diamond is True:
        if os.path.isfile(pwd_cog_20_fa_diamond) is False:
            print(datetime.now().strftime(time_format) + 'DB file for diamond not found, please refers to the help info for diamond db preparation')
            print(datetime.now().strftime(time_format) + 'Program exited!')
            exit()


    ################################################# read db into dict ################################################

    # get protein_to_cog_dict (cog-20.cog.csv)
    protein_to_cog_dict = {}
    for each_line in open(pwd_cog_20_cog_csv):
        each_line_split = each_line.strip().split(',')
        protein_id = each_line_split[2]
        protein_id_no_dot = '_'.join(protein_id.split('.'))
        cog_id = each_line_split[6]
        if protein_id_no_dot not in protein_to_cog_dict:
            protein_to_cog_dict[protein_id_no_dot] = {cog_id}
        else:
            protein_to_cog_dict[protein_id_no_dot].add(cog_id)

    # get cog_id_to_category_dict and cog_id_to_description_dict (cognames2003-2014.tab)
    cog_id_to_category_dict = {}
    cog_id_to_description_dict = {}
    for cog_id_to_cate_des in open(pwd_cog_20_def_tab, encoding='windows-1252'):
        if not cog_id_to_cate_des.startswith('#'):
            cog_id_to_cate_des_split = cog_id_to_cate_des.strip().split('\t')
            cog_id = cog_id_to_cate_des_split[0]
            cog_cate = cog_id_to_cate_des_split[1]
            cog_des = cog_id_to_cate_des_split[2]
            cog_id_to_category_dict[cog_id] = cog_cate
            cog_id_to_description_dict[cog_id] = cog_des

    # get cog_category_to_description_dict (fun2003-2014.tab)
    cog_category_list = []
    cog_category_to_description_dict = {}
    for cog_category in open(pwd_fun_20_tab):
        if not cog_category.startswith('#'):
            cog_category_split = cog_category.strip().split('\t')
            cog_category_list.append(cog_category_split[0])
            cog_category_to_description_dict[cog_category_split[0]] = cog_category_split[1]


    ################################################## if input is file ################################################

    # if input is file
    if os.path.isfile(file_in) is True:

        # check whether depth file exist
        if depth_file is not None:
            if os.path.isfile(depth_file) is False:
                print(datetime.now().strftime(time_format) + 'specified depth file not found, program exited!')
                exit()

        print(datetime.now().strftime(time_format) + 'Running COG annotation for 1 file with %s cores' % (num_threads))

        file_in_path, file_in_basename, file_in_ext = sep_path_basename_ext(file_in)

        COG2020_worker([file_in,
                        pwd_cog_20_fa,
                        protein_to_cog_dict,
                        cog_id_to_category_dict,
                        cog_id_to_description_dict,
                        cog_category_list,
                        cog_category_to_description_dict,
                        sequence_type,
                        file_in_path,
                        num_threads,
                        run_diamond,
                        evalue_cutoff,
                        depth_file,
                        pct_by_all])


    ################################################ if input is folder ################################################

    # if input is folder
    else:

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

            output_folder =                             '%s_COG2020_wd'                             % file_in_folder_name

            # create output folder
            force_create_folder(output_folder)


            ######################################################### main #########################################################

            print(datetime.now().strftime(time_format) + 'Running COG annotation for %s files with %s cores' % (len(input_file_name_list), num_threads))

            list_for_multiple_arguments_COG = []
            for input_file in input_file_name_list:

                input_file_basename = '.'.join(input_file.split('.')[:-1])
                pwd_input_file = '%s/%s' % (file_in, input_file)


                # get path to current depth file
                if depth_file is None:
                    input_file_depth = None
                else:
                    input_file_depth = '%s/%s.depth' % (depth_file, input_file_basename)

                list_for_multiple_arguments_COG.append([pwd_input_file,
                                                        pwd_cog_20_fa,
                                                        protein_to_cog_dict,
                                                        cog_id_to_category_dict,
                                                        cog_id_to_description_dict,
                                                        cog_category_list,
                                                        cog_category_to_description_dict,
                                                        sequence_type,
                                                        output_folder,
                                                        1,
                                                        run_diamond,
                                                        evalue_cutoff,
                                                        input_file_depth,
                                                        pct_by_all])

            # run COG annotaion files with multiprocessing
            pool = mp.Pool(processes=num_threads)
            pool.map(COG2020_worker, list_for_multiple_arguments_COG)
            pool.close()
            pool.join()

            ######################################################### get dataframe #########################################################

            annotation_df_cog_cate_GeneNumber =             '%s/%s_COG2020_cate_GeneNumber.txt'             % (output_folder, file_in_folder_name)
            annotation_df_cog_cate_GeneNumber_pct =         '%s/%s_COG2020_cate_GeneNumber_pct.txt'         % (output_folder, file_in_folder_name)
            annotation_df_cog_cate_GeneNumber_pct_by_all =  '%s/%s_COG2020_cate_GeneNumber_pct_by_all.txt'  % (output_folder, file_in_folder_name)

            annotation_df_cog_cate_TotalDepth =             '%s/%s_COG2020_cate_TotalDepth.txt'             % (output_folder, file_in_folder_name)
            annotation_df_cog_cate_TotalDepth_pct =         '%s/%s_COG2020_cate_TotalDepth_pct.txt'         % (output_folder, file_in_folder_name)
            annotation_df_cog_cate_TotalDepth_pct_by_all =  '%s/%s_COG2020_cate_TotalDepth_pct_by_all.txt'  % (output_folder, file_in_folder_name)

            annotation_df_cog_id_GeneNumber =               '%s/%s_COG2020_id_GeneNumber.txt'               % (output_folder, file_in_folder_name)
            annotation_df_cog_id_GeneNumber_pct =           '%s/%s_COG2020_id_GeneNumber_pct.txt'           % (output_folder, file_in_folder_name)
            annotation_df_cog_id_GeneNumber_pct_by_all =    '%s/%s_COG2020_id_GeneNumber_pct_by_all.txt'    % (output_folder, file_in_folder_name)

            annotation_df_cog_id_TotalDepth =               '%s/%s_COG2020_id_TotalDepth.txt'               % (output_folder, file_in_folder_name)
            annotation_df_cog_id_TotalDepth_pct =           '%s/%s_COG2020_id_TotalDepth_pct.txt'           % (output_folder, file_in_folder_name)
            annotation_df_cog_id_TotalDepth_pct_by_all =    '%s/%s_COG2020_id_TotalDepth_pct_by_all.txt'    % (output_folder, file_in_folder_name)

            print(datetime.now().strftime(time_format) + 'Data matrix exported to:')

            # get df
            get_COG_annot_df(output_folder, 'cog_cate', annotation_df_cog_cate_GeneNumber, annotation_df_cog_cate_GeneNumber_pct, annotation_df_cog_cate_GeneNumber_pct_by_all, with_depth=False, pct_by_all=False)
            get_COG_annot_df(output_folder, 'cog_id', annotation_df_cog_id_GeneNumber, annotation_df_cog_id_GeneNumber_pct, annotation_df_cog_id_GeneNumber_pct_by_all, with_depth=False, pct_by_all=False)
            if pct_by_all is True:
                get_COG_annot_df(output_folder, 'cog_cate', annotation_df_cog_cate_GeneNumber, annotation_df_cog_cate_GeneNumber_pct, annotation_df_cog_cate_GeneNumber_pct_by_all, with_depth=False, pct_by_all=True)
                get_COG_annot_df(output_folder, 'cog_id', annotation_df_cog_id_GeneNumber, annotation_df_cog_id_GeneNumber_pct, annotation_df_cog_id_GeneNumber_pct_by_all, with_depth=False, pct_by_all=True)

            # report
            if pct_by_all is False:
                print(datetime.now().strftime(time_format) + '%s and %s' % (annotation_df_cog_id_GeneNumber.split('/')[-1], annotation_df_cog_id_GeneNumber_pct.split('/')[-1]))
                print(datetime.now().strftime(time_format) + '%s and %s' % (annotation_df_cog_cate_GeneNumber.split('/')[-1], annotation_df_cog_cate_GeneNumber_pct.split('/')[-1]))
            else:
                print(datetime.now().strftime(time_format) + '%s, %s and %s' % (annotation_df_cog_id_GeneNumber.split('/')[-1], annotation_df_cog_id_GeneNumber_pct.split('/')[-1], annotation_df_cog_id_GeneNumber_pct_by_all.split('/')[-1]))
                print(datetime.now().strftime(time_format) + '%s, %s and %s' % (annotation_df_cog_cate_GeneNumber.split('/')[-1], annotation_df_cog_cate_GeneNumber_pct.split('/')[-1], annotation_df_cog_cate_GeneNumber_pct_by_all.split('/')[-1]))

            if depth_file is not None:
                get_COG_annot_df(output_folder, 'cog_cate', annotation_df_cog_cate_TotalDepth, annotation_df_cog_cate_TotalDepth_pct, annotation_df_cog_cate_TotalDepth_pct_by_all, with_depth=True, pct_by_all=False)
                get_COG_annot_df(output_folder, 'cog_id', annotation_df_cog_id_TotalDepth, annotation_df_cog_id_TotalDepth_pct, annotation_df_cog_id_TotalDepth_pct_by_all, with_depth=True, pct_by_all=False)
                if pct_by_all is True:
                    get_COG_annot_df(output_folder, 'cog_cate', annotation_df_cog_cate_TotalDepth, annotation_df_cog_cate_TotalDepth_pct, annotation_df_cog_cate_TotalDepth_pct_by_all, with_depth=True, pct_by_all=True)
                    get_COG_annot_df(output_folder, 'cog_id', annotation_df_cog_id_TotalDepth, annotation_df_cog_id_TotalDepth_pct, annotation_df_cog_id_TotalDepth_pct_by_all, with_depth=True, pct_by_all=True)

                # report
                if pct_by_all is False:
                    print(datetime.now().strftime(time_format) + '%s and %s' % (annotation_df_cog_id_TotalDepth.split('/')[-1], annotation_df_cog_id_TotalDepth_pct.split('/')[-1]))
                    print(datetime.now().strftime(time_format) + '%s and %s' % (annotation_df_cog_cate_TotalDepth.split('/')[-1], annotation_df_cog_cate_TotalDepth_pct.split('/')[-1]))
                else:
                    print(datetime.now().strftime(time_format) + '%s, %s and %s' % (annotation_df_cog_id_TotalDepth.split('/')[-1], annotation_df_cog_id_TotalDepth_pct.split('/')[-1], annotation_df_cog_id_TotalDepth_pct_by_all.split('/')[-1]))
                    print(datetime.now().strftime(time_format) + '%s, %s and %s' % (annotation_df_cog_cate_TotalDepth.split('/')[-1], annotation_df_cog_cate_TotalDepth_pct.split('/')[-1], annotation_df_cog_cate_TotalDepth_pct_by_all.split('/')[-1]))

    ################################################## Final report ####################################################

    print(datetime.now().strftime(time_format) + 'Done!')


if __name__ == '__main__':

    COG_parser = argparse.ArgumentParser()

    # arguments for COG_parser
    COG_parser.add_argument('-i',               required=True,                              help='path to input sequences (in multi-fasta format)')
    COG_parser.add_argument('-x',               required=False,                             help='file extension')
    COG_parser.add_argument('-m',               required=True,                              help='sequence type, "N/n" for "nucleotide", "P/p" for "protein"')
    COG_parser.add_argument('-depth',           required=False, default=None,               help='gene depth file/folder')
    COG_parser.add_argument('-pct_by_all',      required=False, action='store_true',        help='normalize by all query genes, including those without COG assignment')
    COG_parser.add_argument('-db_dir',          required=True,                              help='DB folder')
    COG_parser.add_argument('-diamond',         required=False, action='store_true',        help='run diamond (for big dataset), default is NCBI blastp')
    COG_parser.add_argument('-t',               required=False, type=int, default=1,        help='number of threads')
    COG_parser.add_argument('-evalue',          required=False, default=0.001, type=float,  help='evalue cutoff, default: 0.001')

    args = vars(COG_parser.parse_args())

    COG2020(args)
