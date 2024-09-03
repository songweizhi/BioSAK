import os
import glob
import argparse
from Bio import SeqIO
import multiprocessing as mp
from datetime import datetime
from Bio.SeqRecord import SeqRecord


arCOG_parser_usage = '''
====================================== arCOG example commands ======================================

# Dependencies: blast+, diamond

# annotate protein sequences
BioSAK arCOG -m P -t 6 -db_dir arCOG18_db_dir -i genes.faa -o op_dir
BioSAK arCOG -m P -t 6 -db_dir arCOG18_db_dir -i faa_files -x faa -d depth_files -o op_dir

# annotate DNA sequences (ORFs)
BioSAK arCOG -m N -t 6 -db_dir arCOG18_db_dir -i genes.ffn -d gene.depth -o op_dir
BioSAK arCOG -m N -t 6 -db_dir arCOG18_db_dir -i ffn_files -x ffn -o op_dir

# Download the following files to your database directory:
https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab
https://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/tmp.ar18/ar18.fa.gz
https://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/tmp.ar18/arCOGdef.tab
https://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/tmp.ar18/ar18.ar14.02.csv
gunzip ar18.fa.gz
makeblastdb -in ar18.fa -dbtype prot -parse_seqids -logfile ar18.fa.log
diamond makedb --in ar18.fa --db ar18.fa.dmnd --quiet

# For how it works, please refers to the help information of the COG2020 module (BioSAK COG2020 -h)

====================================================================================================
'''

time_format = '[%Y-%m-%d %H:%M:%S] '


def sep_path_basename_ext(file_in):

    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'
    file_basename, file_ext = os.path.splitext(file_name)
    return file_path, file_basename, file_ext


def AnnotateNorm(file_in, skip_header, value_column, Divisor_value, file_out, file_out_header):

    file_out_handle = open(file_out, 'w')
    file_out_handle.write(file_out_header)
    line_num = 0
    for each_line in open(file_in):
        each_line_split = each_line.strip().split('\t')
        value_str = each_line_split[value_column - 1]
        if (skip_header is True and line_num > 0) or (skip_header is False):
            value_pct = 0
            if Divisor_value != 0:
                value_pct = float(value_str) * 100 / Divisor_value
            each_line_split[value_column - 1] = str(float("{0:.2f}".format(value_pct)))
            file_out_handle.write('%s\n' % '\t'.join(each_line_split))
        line_num += 1
    file_out_handle.close()


def get_gene_list_depth(gene_list, gene_to_depth_dict):

    total_depth = 0
    for gene in gene_list:
        gene_depth = gene_to_depth_dict[gene]
        total_depth += gene_depth
    return total_depth


def dna2aa(dna_file, aa_file):
    query_aa_handle = open(aa_file, 'w')
    for each in SeqIO.parse(dna_file, 'fasta'):
        each_aa = each.seq.translate()
        each_aa_record = SeqRecord(each_aa)
        each_aa_record.id = each.id
        each_aa_record.description = each.description
        SeqIO.write(each_aa_record, query_aa_handle, 'fasta')
    query_aa_handle.close()


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


def arCOG_worker(argument_list):

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
    current_output_folder = '%s/%s_arCOG_wd' % (output_folder, input_seq_no_ext)

    pwd_blastp_output =                 '%s/%s_blastp.tab'                          % (current_output_folder, input_seq_no_ext)
    pwd_blastp_output_besthits =        '%s/%s_blastp_besthits.tab'                 % (current_output_folder, input_seq_no_ext)
    pwd_query_to_cog_txt =              '%s/%s_query_to_cog.txt'                    % (current_output_folder, input_seq_no_ext)
    pwd_cog_stats_copy =                '%s/%s_arcog_stats.txt'                     % (current_output_folder, input_seq_no_ext)
    pwd_cog_stats_depth =               '%s/%s_arcog_stats_depth.txt'               % (current_output_folder, input_seq_no_ext)
    pwd_func_stats_copy =               '%s/%s_func_stats.txt'                      % (current_output_folder, input_seq_no_ext)
    pwd_func_stats_depth =              '%s/%s_func_stats_depth.txt'                % (current_output_folder, input_seq_no_ext)
    pwd_cog_stats_copy_pct =            '%s/%s_arcog_stats_pct.txt'                 % (current_output_folder, input_seq_no_ext)
    pwd_cog_stats_depth_pct =           '%s/%s_arcog_stats_depth_pct.txt'           % (current_output_folder, input_seq_no_ext)
    pwd_func_stats_copy_pct =           '%s/%s_func_stats_pct.txt'                  % (current_output_folder, input_seq_no_ext)
    pwd_func_stats_depth_pct =          '%s/%s_func_stats_depth_pct.txt'            % (current_output_folder, input_seq_no_ext)
    pwd_cog_stats_copy_pct_by_all =     '%s/%s_arcog_stats_pct_by_all.txt'          % (current_output_folder, input_seq_no_ext)
    pwd_cog_stats_depth_pct_by_all =    '%s/%s_arcog_stats_depth_pct_by_all.txt'    % (current_output_folder, input_seq_no_ext)
    pwd_func_stats_copy_pct_by_all =    '%s/%s_func_stats_pct_by_all.txt'           % (current_output_folder, input_seq_no_ext)
    pwd_func_stats_depth_pct_by_all =   '%s/%s_func_stats_depth_pct_by_all.txt'     % (current_output_folder, input_seq_no_ext)

    os.system('mkdir %s' % current_output_folder)

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
    best_hit(pwd_blastp_output, pwd_blastp_output_besthits)

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
    cog_cate_num_dict = {}
    cog_id_to_gene_member_dict = {}
    cog_cate_to_gene_member_dict = {}
    genes_with_cog = set()
    pwd_query_to_cog_txt_handle = open(pwd_query_to_cog_txt, 'w')
    pwd_query_to_cog_txt_handle.write('Query\tarCOG\tCategory\tDescription\n')
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
                    cog_cate = cog_id_to_category_dict.get(cog_id, 'na')
                    cog_des = cog_id_to_description_dict[cog_id]
                    pwd_query_to_cog_txt_handle.write('%s\t%s\t%s\t%s\n' % (query_gene, cog_id, ''.join(cog_cate), cog_des))
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

    # get depth of all query genes or genes with cog assignment
    if depth_file is not None:
        genes_with_cog_depth = get_gene_list_depth(genes_with_cog, gene_depth_dict)
        total_depth_for_all_query_genes = get_gene_list_depth(query_seq_list, gene_depth_dict)

    #################### export cog_stats_copy ####################

    pwd_cog_stats_copy_handle = open(pwd_cog_stats_copy, 'w')
    pwd_cog_stats_copy_handle.write('arCOG\tcopy\tDescription\n')
    for each_cog_id in cog_id_num_dict:
        each_cog_id_copy = cog_id_num_dict[each_cog_id]
        pwd_cog_stats_copy_handle.write('%s\t%s\t%s\n' % (each_cog_id, each_cog_id_copy, cog_id_to_description_dict[each_cog_id]))
    pwd_cog_stats_copy_handle.close()

    #################### export cog_stats_depth ####################

    if depth_file is not None:
        pwd_cog_stats_depth_handle = open(pwd_cog_stats_depth, 'w')
        pwd_cog_stats_depth_handle.write('COG\tdepth\tDescription\n')
        for each_cog_id in cog_id_to_gene_member_dict:
            each_cog_id_gene_member = cog_id_to_gene_member_dict[each_cog_id]
            each_cog_id_depth = 0
            for each_gene in each_cog_id_gene_member:
                each_gene_depth = gene_depth_dict[each_gene]
                each_cog_id_depth += each_gene_depth
            each_cog_id_depth = float("{0:.2f}".format(each_cog_id_depth))
            pwd_cog_stats_depth_handle.write('%s\t%s\t%s\n' % (each_cog_id, each_cog_id_depth, cog_id_to_description_dict[each_cog_id]))
        pwd_cog_stats_depth_handle.close()

    #################### export func_stats_copy ####################

    pwd_func_stats_copy_handle = open(pwd_func_stats_copy, 'w')
    pwd_func_stats_copy_handle.write('Category\tcopy\tDescription\n')
    for each_cog_cate in cog_category_list:
        each_cog_cate_copy = 0
        if each_cog_cate in cog_cate_num_dict:
            each_cog_cate_copy = cog_cate_num_dict[each_cog_cate]
        pwd_func_stats_copy_handle.write('%s\t%s\t%s\n' % (each_cog_cate, each_cog_cate_copy, cog_category_to_description_dict.get(each_cog_cate, 'na')))
    pwd_func_stats_copy_handle.close()

    #################### export func_stats_depth ####################

    if depth_file is not None:
        pwd_func_stats_depth_handle = open(pwd_func_stats_depth, 'w')
        pwd_func_stats_depth_handle.write('Category\tdepth\tDescription\n')
        for each_cog_cate in cog_category_list:
            each_cog_cate_depth = 0
            if each_cog_cate in cog_cate_to_gene_member_dict:
                each_cog_cate_gene_member = cog_cate_to_gene_member_dict[each_cog_cate]
                for each_gene in each_cog_cate_gene_member:
                    each_gene_depth = gene_depth_dict[each_gene]
                    each_cog_cate_depth += each_gene_depth
            each_cog_cate_depth = float("{0:.2f}".format(each_cog_cate_depth))
            pwd_func_stats_depth_handle.write('%s\t%s\t%s\n' % (each_cog_cate, each_cog_cate_depth, cog_category_to_description_dict[each_cog_cate]))
        pwd_func_stats_depth_handle.close()

    #################### get pct files ####################

    AnnotateNorm(file_in=pwd_cog_stats_copy,  skip_header=True, value_column=2, Divisor_value=len(genes_with_cog), file_out=pwd_cog_stats_copy_pct,  file_out_header='arCOG\tPercent\tDescription\n')
    AnnotateNorm(file_in=pwd_func_stats_copy, skip_header=True, value_column=2, Divisor_value=len(genes_with_cog), file_out=pwd_func_stats_copy_pct, file_out_header='Category\tPercent\tDescription\n')
    if depth_file is not None:
        AnnotateNorm(file_in=pwd_cog_stats_depth,  skip_header=True, value_column=2, Divisor_value=genes_with_cog_depth, file_out=pwd_cog_stats_depth_pct,  file_out_header='arCOG\tPercent\tDescription\n')
        AnnotateNorm(file_in=pwd_func_stats_depth, skip_header=True, value_column=2, Divisor_value=genes_with_cog_depth, file_out=pwd_func_stats_depth_pct, file_out_header='Category\tPercent\tDescription\n')
    if pct_by_all is True:
        AnnotateNorm(file_in=pwd_cog_stats_copy,  skip_header=True, value_column=2, Divisor_value=len(query_seq_list), file_out=pwd_cog_stats_copy_pct_by_all,  file_out_header='arCOG\tPercent_by_all\tDescription\n')
        AnnotateNorm(file_in=pwd_func_stats_copy, skip_header=True, value_column=2, Divisor_value=len(query_seq_list), file_out=pwd_func_stats_copy_pct_by_all, file_out_header='Category\tPercent_by_all\tDescription\n')
        if depth_file is not None:
            AnnotateNorm(file_in=pwd_cog_stats_depth,  skip_header=True, value_column=2, Divisor_value=total_depth_for_all_query_genes, file_out=pwd_cog_stats_depth_pct_by_all,  file_out_header='arCOG\tPercent_by_all\tDescription\n')
            AnnotateNorm(file_in=pwd_func_stats_depth, skip_header=True, value_column=2, Divisor_value=total_depth_for_all_query_genes, file_out=pwd_func_stats_depth_pct_by_all, file_out_header='Category\tPercent_by_all\tDescription\n')


def get_COG_annot_df(annotation_dir, stats_level, annotation_df_absolute_num, annotation_df_percentage, annotation_df_percentage_by_all, with_depth, pct_by_all):

    annotation_dir_re = '%s/*_arCOG_wd' % annotation_dir
    annotation_folder_list = [os.path.basename(file_name) for file_name in glob.glob(annotation_dir_re)]

    cog_num_dict = {}
    cog_num_pct_dict = {}
    cog_num_pct_by_all_dict = {}
    all_identified_cog = set()
    for annotation_folder in annotation_folder_list:
        annotation_folder_basename = annotation_folder.split('_arCOG_wd')[0]
        pwd_annotation_stats_file = ''
        pwd_annotation_stats_file_pct = ''

        if stats_level == 'cog_id':
            if with_depth is False:
                pwd_annotation_stats_file =             '%s/%s/%s_arcog_stats.txt'                  % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct =         '%s/%s/%s_arcog_stats_pct.txt'              % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct_by_all =  '%s/%s/%s_arcog_stats_pct_by_all.txt'       % (annotation_dir, annotation_folder, annotation_folder_basename)
            else:
                pwd_annotation_stats_file =             '%s/%s/%s_arcog_stats_depth.txt'            % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct =         '%s/%s/%s_arcog_stats_depth_pct.txt'        % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct_by_all =  '%s/%s/%s_arcog_stats_depth_pct_by_all.txt' % (annotation_dir, annotation_folder, annotation_folder_basename)

        if stats_level == 'cog_cate':
            if with_depth is False:
                pwd_annotation_stats_file =             '%s/%s/%s_func_stats.txt'                   % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct =         '%s/%s/%s_func_stats_pct.txt'               % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct_by_all =  '%s/%s/%s_func_stats_pct_by_all.txt'        % (annotation_dir, annotation_folder, annotation_folder_basename)
            else:
                pwd_annotation_stats_file =             '%s/%s/%s_func_stats_depth.txt'             % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct =         '%s/%s/%s_func_stats_depth_pct.txt'         % (annotation_dir, annotation_folder, annotation_folder_basename)
                pwd_annotation_stats_file_pct_by_all =  '%s/%s/%s_func_stats_depth_pct_by_all.txt'  % (annotation_dir, annotation_folder, annotation_folder_basename)

        current_cog_to_num_dict = {}
        for cog in open(pwd_annotation_stats_file):
            if (not cog.startswith('Category	')) and (not cog.startswith('COG	')) and (not cog.startswith('arCOG	')):
                cog_split = cog.strip().split('\t')
                if with_depth is False:
                    current_cog_to_num_dict[cog_split[0]] = int(cog_split[1])
                else:
                    current_cog_to_num_dict[cog_split[0]] = float(cog_split[1])
                all_identified_cog.add(cog_split[0])
        cog_num_dict[annotation_folder_basename] = current_cog_to_num_dict

        current_cog_to_num_pct_dict = {}
        for cog in open(pwd_annotation_stats_file_pct):
            if (not cog.startswith('Category	')) and (not cog.startswith('COG	')) and (not cog.startswith('arCOG	')):
                cog_split = cog.strip().split('\t')
                current_cog_to_num_pct_dict[cog_split[0]] = float(cog_split[1])
                all_identified_cog.add(cog_split[0])
        cog_num_pct_dict[annotation_folder_basename] = current_cog_to_num_pct_dict

        if pct_by_all is True:
            current_cog_to_num_pct_by_all_dict = {}
            for cog in open(pwd_annotation_stats_file_pct_by_all):
                if (not cog.startswith('Category	')) and (not cog.startswith('COG	')) and (not cog.startswith('arCOG	')):
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
        annotation_folder_basename = annotation_folder.split('_arCOG_wd')[0]
        current_cog_num_dict = cog_num_dict[annotation_folder_basename]
        current_cog_num_pct_dict = cog_num_pct_dict[annotation_folder_basename]

        current_cog_num_list = []
        current_cog_num_list_pct = []
        for identified_cog in all_identified_cog_list:
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
            annotation_folder_basename = annotation_folder.split('_arCOG_wd')[0]
            current_cog_num_pct_by_all_dict = cog_num_pct_by_all_dict[annotation_folder_basename]

            current_cog_num_list_pct_by_all = []
            for identified_cog in all_identified_cog_list:
                identified_cog_num_pct_by_all = 0
                if identified_cog in current_cog_num_pct_by_all_dict:
                    identified_cog_num_pct_by_all = current_cog_num_pct_by_all_dict[identified_cog]
                current_cog_num_list_pct_by_all.append(identified_cog_num_pct_by_all)

            # write out
            annotation_df_percentage_by_all_handle.write('%s\t%s\n' % (annotation_folder_basename, '\t'.join([str(i) for i in current_cog_num_list_pct_by_all])))
        annotation_df_percentage_by_all_handle.close()


def arCOG(args):

    file_in             = args['i']
    output_folder       = args['o']
    file_extension      = args['x']
    sequence_type       = args['m']
    depth_file          = args['d']
    pct_by_all          = args['pct_by_all']
    DB_dir              = args['db_dir']
    num_threads         = args['t']
    run_diamond         = args['diamond']
    evalue_cutoff       = args['evalue']
    force_create_op_dir = args['f']

    ar18_fa           = '%s/ar18.fa'                 % DB_dir
    ar18_fa_diamond   = '%s/ar18.fa.dmnd'            % DB_dir
    ar18_ar14_02_csv  = '%s/ar18.ar14.02.csv'        % DB_dir
    cog_des_txt       = '%s/arCOGdef.tab'            % DB_dir
    pwd_fun_20_tab    = '%s/fun-20.tab'              % DB_dir

    ############################################ check whether db file exist ###########################################

    # check whether db file exist
    unfound_inputs = []
    for each_input in [ar18_ar14_02_csv, cog_des_txt, ar18_fa, pwd_fun_20_tab]:
        if (not os.path.isfile(each_input)) and (not os.path.isdir(each_input)):
            unfound_inputs.append(each_input)
    if len(unfound_inputs) > 0:
        for each_unfound in unfound_inputs:
            print('%s not found' % each_unfound)
        exit()

    if run_diamond is True:
        if os.path.isfile(ar18_fa_diamond) is False:
            print(datetime.now().strftime(time_format) + 'DB file for diamond not found, please refers to the help info for diamond db preparation')
            print(datetime.now().strftime(time_format) + 'Program exited!')
            exit()

    ################################################# read db into dict ################################################

    # get cog_id_to_category_dict and cog_id_to_description_dict (arCOGdef.tab)
    cog_category_set = set()
    cog_id_to_category_dict = dict()
    cog_id_to_description_dict = dict()
    for each_cog in open(cog_des_txt, encoding="ISO-8859-1"):
        each_cog_split = each_cog.strip().split('\t')
        cog_id = each_cog_split[0]
        cog_cate_str = each_cog_split[1]
        cog_cate_split = [i for i in cog_cate_str]
        cog_desc = each_cog_split[3]
        cog_id_to_description_dict[cog_id] = cog_desc
        cog_id_to_category_dict[cog_id] = cog_cate_split
        for each_cate in cog_cate_split:
            cog_category_set.add(each_cate)

    cog_category_list = sorted([i for i in cog_category_set])

    # get protein_to_cog_dict (ar18_ar14_02_csv)
    protein_to_cog_dict = {}
    for each_line in open(ar18_ar14_02_csv):
        each_line_split = each_line.strip().split(',')
        arcog_id = each_line_split[6]
        seq_id = each_line_split[2]
        if seq_id not in protein_to_cog_dict:
            protein_to_cog_dict[seq_id] = {arcog_id}
        else:
            protein_to_cog_dict[seq_id].add(arcog_id)

    # get cog_category_to_description_dict (fun-20.tab)
    cog_category_to_description_dict = {}
    for cog_category in open(pwd_fun_20_tab):
        cog_category_split = cog_category.strip().split('\t')
        cog_category_to_description_dict[cog_category_split[0]] = cog_category_split[2]

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
        arCOG_worker([file_in, ar18_fa, protein_to_cog_dict, cog_id_to_category_dict, cog_id_to_description_dict,
                      cog_category_list, cog_category_to_description_dict, sequence_type, file_in_path, num_threads,
                      run_diamond, evalue_cutoff, depth_file, pct_by_all])

    ################################################ if input is folder ################################################

    # if input is folder
    else:
        if os.path.isdir(file_in) is False:
            print(datetime.now().strftime(time_format) + 'input folder not found, program exited')
            exit()
        else:
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

            file_in_folder_name = file_in
            if '/' in file_in:
                file_in_folder_name = file_in.split('/')[-1]

            # create output folder
            if os.path.isdir(output_folder) is True:
                if force_create_op_dir is True:
                    os.system('rm -r %s' % output_folder)
                else:
                    print('Output folder detected, program exited!')
                    exit()
            os.system('mkdir %s' % output_folder)

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

                list_for_multiple_arguments_COG.append([pwd_input_file, ar18_fa, protein_to_cog_dict, cog_id_to_category_dict, cog_id_to_description_dict,
                                                        cog_category_list, cog_category_to_description_dict, sequence_type, output_folder, 1, run_diamond,
                                                        evalue_cutoff, input_file_depth, pct_by_all])

            # run COG annotaion files with multiprocessing
            pool = mp.Pool(processes=num_threads)
            pool.map(arCOG_worker, list_for_multiple_arguments_COG)
            pool.close()
            pool.join()

            ######################################################### get dataframe #########################################################

            annotation_df_cog_cate_copy =             '%s/%s_arCOG_cate.txt'                    % (output_folder, file_in_folder_name)
            annotation_df_cog_cate_copy_pct =         '%s/%s_arCOG_cate_pct.txt'                % (output_folder, file_in_folder_name)
            annotation_df_cog_cate_copy_pct_by_all =  '%s/%s_arCOG_cate_pct_by_all.txt'         % (output_folder, file_in_folder_name)
            annotation_df_cog_cate_depth =            '%s/%s_arCOG_cate_depth.txt'              % (output_folder, file_in_folder_name)
            annotation_df_cog_cate_depth_pct =        '%s/%s_arCOG_cate_depth_pct.txt'          % (output_folder, file_in_folder_name)
            annotation_df_cog_cate_depth_pct_by_all = '%s/%s_arCOG_cate_depth_pct_by_all.txt'   % (output_folder, file_in_folder_name)
            annotation_df_cog_id_copy =               '%s/%s_arCOG_id.txt'                      % (output_folder, file_in_folder_name)
            annotation_df_cog_id_copy_pct =           '%s/%s_arCOG_id_pct.txt'                  % (output_folder, file_in_folder_name)
            annotation_df_cog_id_copy_pct_by_all =    '%s/%s_arCOG_id_pct_by_all.txt'           % (output_folder, file_in_folder_name)
            annotation_df_cog_id_depth =              '%s/%s_arCOG_id_depth.txt'                % (output_folder, file_in_folder_name)
            annotation_df_cog_id_depth_pct =          '%s/%s_arCOG_id_depth_pct.txt'            % (output_folder, file_in_folder_name)
            annotation_df_cog_id_depth_pct_by_all =   '%s/%s_arCOG_id_depth_pct_by_all.txt'     % (output_folder, file_in_folder_name)

            print(datetime.now().strftime(time_format) + 'Data matrix exported to:')

            # get df
            get_COG_annot_df(output_folder, 'cog_cate', annotation_df_cog_cate_copy, annotation_df_cog_cate_copy_pct, annotation_df_cog_cate_copy_pct_by_all, with_depth=False, pct_by_all=False)
            get_COG_annot_df(output_folder, 'cog_id', annotation_df_cog_id_copy, annotation_df_cog_id_copy_pct, annotation_df_cog_id_copy_pct_by_all, with_depth=False, pct_by_all=False)
            if pct_by_all is True:
                get_COG_annot_df(output_folder, 'cog_cate', annotation_df_cog_cate_copy, annotation_df_cog_cate_copy_pct, annotation_df_cog_cate_copy_pct_by_all, with_depth=False, pct_by_all=True)
                get_COG_annot_df(output_folder, 'cog_id', annotation_df_cog_id_copy, annotation_df_cog_id_copy_pct, annotation_df_cog_id_copy_pct_by_all, with_depth=False, pct_by_all=True)

            # report
            if pct_by_all is False:
                print(datetime.now().strftime(time_format) + '%s and %s'        % (annotation_df_cog_id_copy.split('/')[-1], annotation_df_cog_id_copy_pct.split('/')[-1]))
                print(datetime.now().strftime(time_format) + '%s and %s'        % (annotation_df_cog_cate_copy.split('/')[-1], annotation_df_cog_cate_copy_pct.split('/')[-1]))
            else:
                print(datetime.now().strftime(time_format) + '%s, %s and %s'    % (annotation_df_cog_id_copy.split('/')[-1], annotation_df_cog_id_copy_pct.split('/')[-1], annotation_df_cog_id_copy_pct_by_all.split('/')[-1]))
                print(datetime.now().strftime(time_format) + '%s, %s and %s'    % (annotation_df_cog_cate_copy.split('/')[-1], annotation_df_cog_cate_copy_pct.split('/')[-1], annotation_df_cog_cate_copy_pct_by_all.split('/')[-1]))

            if depth_file is not None:
                get_COG_annot_df(output_folder, 'cog_cate', annotation_df_cog_cate_depth, annotation_df_cog_cate_depth_pct, annotation_df_cog_cate_depth_pct_by_all, with_depth=True, pct_by_all=False)
                get_COG_annot_df(output_folder, 'cog_id', annotation_df_cog_id_depth, annotation_df_cog_id_depth_pct, annotation_df_cog_id_depth_pct_by_all, with_depth=True, pct_by_all=False)
                if pct_by_all is True:
                    get_COG_annot_df(output_folder, 'cog_cate', annotation_df_cog_cate_depth, annotation_df_cog_cate_depth_pct, annotation_df_cog_cate_depth_pct_by_all, with_depth=True, pct_by_all=True)
                    get_COG_annot_df(output_folder, 'cog_id', annotation_df_cog_id_depth, annotation_df_cog_id_depth_pct, annotation_df_cog_id_depth_pct_by_all, with_depth=True, pct_by_all=True)

                # report
                if pct_by_all is False:
                    print(datetime.now().strftime(time_format) + '%s and %s' % (annotation_df_cog_id_depth.split('/')[-1], annotation_df_cog_id_depth_pct.split('/')[-1]))
                    print(datetime.now().strftime(time_format) + '%s and %s' % (annotation_df_cog_cate_depth.split('/')[-1], annotation_df_cog_cate_depth_pct.split('/')[-1]))
                else:
                    print(datetime.now().strftime(time_format) + '%s, %s and %s' % (annotation_df_cog_id_depth.split('/')[-1], annotation_df_cog_id_depth_pct.split('/')[-1], annotation_df_cog_id_depth_pct_by_all.split('/')[-1]))
                    print(datetime.now().strftime(time_format) + '%s, %s and %s' % (annotation_df_cog_cate_depth.split('/')[-1], annotation_df_cog_cate_depth_pct.split('/')[-1], annotation_df_cog_cate_depth_pct_by_all.split('/')[-1]))

    ################################################## Final report ####################################################

    print(datetime.now().strftime(time_format) + 'Done!')


if __name__ == '__main__':

    arCOG_parser = argparse.ArgumentParser(usage=arCOG_parser_usage)
    arCOG_parser.add_argument('-i',               required=True,                              help='path to input sequences (in multi-fasta format)')
    arCOG_parser.add_argument('-o',               required=True,                              help='output directory')
    arCOG_parser.add_argument('-x',               required=False,                             help='file extension')
    arCOG_parser.add_argument('-m',               required=True,                              help='sequence type, "N/n" for "nucleotide", "P/p" for "protein"')
    arCOG_parser.add_argument('-d',               required=False, default=None,               help='gene depth file/folder')
    arCOG_parser.add_argument('-db_dir',          required=True,                              help='COG_db_dir')
    arCOG_parser.add_argument('-pct_by_all',      required=False, action='store_true',        help='normalize by all query genes, including those without annotation')
    arCOG_parser.add_argument('-diamond',         required=False, action='store_true',        help='run diamond (for big dataset), default: blastp from NCBI')
    arCOG_parser.add_argument('-evalue',          required=False, default=0.0001,             help='evalue cutoff, accepted format 0.001, 1e-10, 1e-30, default: 0.0001')
    arCOG_parser.add_argument('-t',               required=False, type=int, default=1,        help='number of threads')
    arCOG_parser.add_argument('-f',               required=False, action="store_true",        help='force overwrite')
    args = vars(arCOG_parser.parse_args())
    arCOG(args)
