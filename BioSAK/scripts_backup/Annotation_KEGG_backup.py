import os
import glob
import shutil
import argparse
from Bio import SeqIO
from time import sleep
from datetime import datetime
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing as mp


KEGG_parser_usage = '''
========================================= KEGG_annot example commands ========================================

# module needed
module load python/3.7.3
module load blast+/2.6.0
module load diamond/0.9.24

# run with NCBI blastp (default, for small dataset)
MyBioTools KEGG -t 6 -db_dir /srv/scratch/z5039045/DB/KEGG_DB -seq_in input.faa
MyBioTools KEGG -t 6 -db_dir /srv/scratch/z5039045/DB/KEGG_DB -seq_in faa_files -x faa
MyBioTools KEGG -t 6 -db_dir /srv/scratch/z5039045/DB/KEGG_DB -seq_in faa_files -x faa -plot

# run with Diamond blastp (for big dataset)
MyBioTools KEGG -t 6 -db_dir /srv/scratch/z5039045/DB/KEGG_DB -diamond -seq_in input.faa
MyBioTools KEGG -t 6 -db_dir /srv/scratch/z5039045/DB/KEGG_DB -diamond -seq_in faa_folder -x faa
MyBioTools KEGG -t 6 -db_dir /srv/scratch/z5039045/DB/KEGG_DB -diamond -seq_in faa_folder -x faa -plot

# get summary for BlastKOALA/GhostKOALA produced results (put all [prefix]_user_ko.txt files in input folder)
MyBioTools KEGG -t 6 -db_dir /srv/scratch/z5039045/DB/KEGG_DB -ko_in user_ko.txt
MyBioTools KEGG -t 6 -db_dir /srv/scratch/z5039045/DB/KEGG_DB -ko_in user_ko_folder -x txt
MyBioTools KEGG -t 6 -db_dir /srv/scratch/z5039045/DB/KEGG_DB -ko_in user_ko_folder -x txt -plot

# DB files:
prokaryotes.pep.fasta : refers to https://www.kegg.jp/kegg/download/Readme/README.fasta (no need for "-ko_in" mode)
prokaryotes.dat       : refers to https://www.kegg.jp/kegg/download/Readme/README.fasta
ko00001.keg           : https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir=

==============================================================================================================
'''


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_ext = os.path.splitext(file_name)

    return file_path, file_basename, file_ext


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


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


def write_out_summary(ko_level, list_all, list_uniq, query_gene_all, query_gene_ko_NA, ko_description_dict, summary_file_out):

    summary_file_out_handle = open(summary_file_out, 'w')
    summary_file_out_handle.write('Level\tKO\tNumber\tPercent(%)\tDescription\n')
    for each_uniq_ko in list_uniq:
        each_uniq_ko_ID = each_uniq_ko.split('_')[1]
        each_uniq_ko_count = list_all.count(each_uniq_ko)
        each_uniq_ko_percent = float("{0:.2f}".format(list_all.count(each_uniq_ko) * 100 / len(query_gene_all)))
        each_uniq_ko_description = ko_description_dict[each_uniq_ko_ID]
        summary_file_out_handle.write('%s\t%s\t%s\t%s\t%s\n' % (
        ko_level, each_uniq_ko_ID, each_uniq_ko_count, each_uniq_ko_percent, each_uniq_ko_description))

    summary_file_out_handle.write('%s\tNA\t%s\t%s\tNA\n' % (ko_level, query_gene_ko_NA, float("{0:.2f}".format(query_gene_ko_NA * 100 / len(query_gene_all)))))
    summary_file_out_handle.close()


def barh_plotter(num_list, label_list, query_seq_num, query_ko_NA, fig_width, fig_height, plot_file):

    fig, ax = plt.subplots()
    fig.set_size_inches(fig_width, fig_height)

    y_pos = range(len(num_list))
    ax.barh(y_pos, num_list, height=0.8, align='center', alpha=0.2, linewidth=0)
    ax.set_yticks([])  # not show yticks
    ax.invert_xaxis()  # line up bar on right
    ax.invert_yaxis()  # put first number on top
    ax.axis('tight')   # remove extra spaces at the top and bottom, equal to: ax.margins(0, 0)
    # ax.margins(0, 0.01) # customize space percentage

    ax.set_xlabel('Number of gene')
    ax.set_title('Query genes number: %s, genes without KO: %s' % (query_seq_num, query_ko_NA))

    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(label_list)

    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    plt.close()
    plt.clf()


def Annotation_KEGG_worker(argument_list):

    pwd_input_file =        argument_list[0]
    run_blast =             argument_list[1]
    run_diamond =           argument_list[2]
    KEGG_DB_seq =           argument_list[3]
    KEGG_DB_seq_diamond =   argument_list[4]
    As_description_dict =   argument_list[5]
    Bs_description_dict =   argument_list[6]
    Cs_description_dict =   argument_list[7]
    Ds_description_dict =   argument_list[8]
    D2ABCD_dict =           argument_list[9]
    db_seq_to_KO_dict =     argument_list[10]
    op_dir =                argument_list[11]
    plot_stats =            argument_list[12]

    ################################################### define file name ###################################################

    input_file_path, in_file_basename, input_file_ext = sep_path_basename_ext(pwd_input_file)

    blast_results =           '%s/%s_KEGG_wd/%s_blast.tab'              % (op_dir, in_file_basename, in_file_basename)
    blast_results_best_hit =  '%s/%s_KEGG_wd/%s_blast_best_hits.tab'    % (op_dir, in_file_basename, in_file_basename)

    KO_assignment_file_D =    '%s/%s_KEGG_wd/%s_KO_assignment_D.txt'    % (op_dir, in_file_basename, in_file_basename)
    KO_assignment_file_DCBA = '%s/%s_KEGG_wd/%s_KO_assignment_DCBA.txt' % (op_dir, in_file_basename, in_file_basename)

    stats_file_A =            '%s/%s_KEGG_wd/%s_ko_stats_A.txt'         % (op_dir, in_file_basename, in_file_basename)
    stats_file_B =            '%s/%s_KEGG_wd/%s_ko_stats_B.txt'         % (op_dir, in_file_basename, in_file_basename)
    stats_file_C =            '%s/%s_KEGG_wd/%s_ko_stats_C.txt'         % (op_dir, in_file_basename, in_file_basename)
    stats_file_D =            '%s/%s_KEGG_wd/%s_ko_stats_D.txt'         % (op_dir, in_file_basename, in_file_basename)

    stats_plot_A =            '%s/%s_KEGG_wd/%s_ko_stats_A.png'         % (op_dir, in_file_basename, in_file_basename)
    stats_plot_B =            '%s/%s_KEGG_wd/%s_ko_stats_B.png'         % (op_dir, in_file_basename, in_file_basename)
    stats_plot_C =            '%s/%s_KEGG_wd/%s_ko_stats_C.png'         % (op_dir, in_file_basename, in_file_basename)
    stats_plot_D =            '%s/%s_KEGG_wd/%s_ko_stats_D.png'         % (op_dir, in_file_basename, in_file_basename)


    # create output folder
    force_create_folder('%s/%s_KEGG_wd' % (op_dir, in_file_basename))


    ########################################## blast against KEGG database (Shan) ##########################################

    if run_blast is True:

        if run_diamond is False:
            blastp_cmd = 'blastp -query %s -db %s -out %s -outfmt 6 -evalue 0.001 -num_alignments 10 -num_threads 1' % (pwd_input_file, KEGG_DB_seq, blast_results)
            os.system(blastp_cmd)

        else:
            diamond_cmd = 'diamond blastp -q %s --db %s --out %s --outfmt 6 --evalue 0.001 --threads 1 --quiet' % (pwd_input_file, KEGG_DB_seq_diamond, blast_results)
            os.system(diamond_cmd)

        # only keep the best hit
        keep_blast_hit_with_highest_bit_score(blast_results, blast_results_best_hit)

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

    if run_blast is False:
        KO_assignment_file_D = pwd_input_file

    # get ko id at all levels for all query genes
    ko_assig_ABCD_handle = open(KO_assignment_file_DCBA, 'w')
    ko_assig_ABCD_handle.write('Gene_id\tID_D\tID_C\tID_B\tID_A\tDesc_D\tDesc_C\tDesc_B\tDesc_A\n')
    query_seq_id_all = set()
    for query_gene in open(KO_assignment_file_D):
        query_gene_split = query_gene.strip().split('\t')
        gene_ID = query_gene_split[0]

        if len(query_gene_split) == 1:
            query_seq_id_all.add(query_gene_split[0])
            ko_assig_ABCD_handle.write('%s\n' % gene_ID)

        if len(query_gene_split) == 2:
            query_seq_id_all.add(query_gene_split[0])
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
                    ko_assig_ABCD_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (gene_ID,
                                                                             '\t'.join(KO_DCBA_list),
                                                                             desc_D, desc_C, desc_B, desc_A))

                if len(KO_ID_ABCD) > 1:
                    for each_ABCD in KO_ID_ABCD:
                        each_KO_DCBA_list = each_ABCD.split('|')[::-1]
                        each_KO_DCBA_list_only_id = [i.split('_')[1] for i in each_KO_DCBA_list]
                        each_desc_A = As_description_dict[each_KO_DCBA_list_only_id[3]]
                        each_desc_B = Bs_description_dict[each_KO_DCBA_list_only_id[2]]
                        each_desc_C = Cs_description_dict[each_KO_DCBA_list_only_id[1]]
                        each_desc_D = Ds_description_dict[each_KO_DCBA_list_only_id[0]]
                        ko_assig_ABCD_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (gene_ID,
                                                                                 '\t'.join(each_KO_DCBA_list),
                                                                                 each_desc_D,
                                                                                 each_desc_C,
                                                                                 each_desc_B,
                                                                                 each_desc_A))

    ko_assig_ABCD_handle.close()


    ##################################################### Get summary ######################################################

    query_As_dict = {}
    query_Bs_dict = {}
    query_Cs_dict = {}
    query_Ds_dict = {}
    query_ko_NA = 0
    for each_query in open(KO_assignment_file_DCBA):
        if not each_query.startswith('Gene_id'):
            each_query_split = each_query.strip().split('\t')

            if len(each_query_split) == 1:
                query_ko_NA += 1

            if len(each_query_split) > 1:
                query_id = each_query_split[0]
                query_ko_A = each_query_split[4]
                query_ko_B = each_query_split[3]
                query_ko_C = each_query_split[2]
                query_ko_D = each_query_split[1]

                # get all KOs at level A assigned to current query
                if query_id not in query_As_dict:
                    query_As_dict[query_id] = {query_ko_A}
                else:
                    query_As_dict[query_id].add(query_ko_A)

                # get all KOs at level B assigned to current query
                if query_id not in query_Bs_dict:
                    query_Bs_dict[query_id] = {query_ko_B}
                else:
                    query_Bs_dict[query_id].add(query_ko_B)

                # get all KOs at level C assigned to current query
                if query_id not in query_Cs_dict:
                    query_Cs_dict[query_id] = {query_ko_C}
                else:
                    query_Cs_dict[query_id].add(query_ko_C)

                # get all KOs at level D assigned to current query
                if query_id not in query_Ds_dict:
                    query_Ds_dict[query_id] = {query_ko_D}
                else:
                    query_Ds_dict[query_id].add(query_ko_D)

    As_list = []
    for query_A in query_As_dict:
        query_A_kos = query_As_dict[query_A]
        for each_query_A_ko in query_A_kos:
            As_list.append(each_query_A_ko)

    Bs_list = []
    for query_B in query_Bs_dict:
        query_B_kos = query_Bs_dict[query_B]
        for each_query_B_ko in query_B_kos:
            Bs_list.append(each_query_B_ko)

    Cs_list = []
    for query_C in query_Cs_dict:
        query_C_kos = query_Cs_dict[query_C]
        for each_query_C_ko in query_C_kos:
            Cs_list.append(each_query_C_ko)

    Ds_list = []
    for query_D in query_Ds_dict:
        query_D_kos = query_Ds_dict[query_D]
        for eaDh_query_D_ko in query_D_kos:
            Ds_list.append(eaDh_query_D_ko)

    # uniq list
    As_list_uniq = sorted(unique_list_elements(As_list))
    Bs_list_uniq = sorted(unique_list_elements(Bs_list))
    Cs_list_uniq = sorted(unique_list_elements(Cs_list))
    Ds_list_uniq = sorted(unique_list_elements(Ds_list))

    # write out summary
    write_out_summary('A', As_list, As_list_uniq, query_seq_id_all, query_ko_NA, As_description_dict, stats_file_A)
    write_out_summary('B', Bs_list, Bs_list_uniq, query_seq_id_all, query_ko_NA, Bs_description_dict, stats_file_B)
    write_out_summary('C', Cs_list, Cs_list_uniq, query_seq_id_all, query_ko_NA, Cs_description_dict, stats_file_C)
    write_out_summary('D', Ds_list, Ds_list_uniq, query_seq_id_all, query_ko_NA, Ds_description_dict, stats_file_D)


    ################################################## Plot distribution ###################################################

    if plot_stats is True:
        ko_num_list_A = [As_list.count(i) for i in As_list_uniq]
        ko_num_list_B = [Bs_list.count(i) for i in Bs_list_uniq]
        ko_num_list_C = [Cs_list.count(i) for i in Cs_list_uniq]

        As_list_uniq_desc = [('A %s %s' % (i.split('_')[1], As_description_dict[i.split('_')[1]])) for i in As_list_uniq]
        Bs_list_uniq_desc = [('B %s %s' % (i.split('_')[1], Bs_description_dict[i.split('_')[1]])) for i in Bs_list_uniq]
        Cs_list_uniq_desc = [('C %s %s' % (i.split('_')[1], Cs_description_dict[i.split('_')[1]])) for i in Cs_list_uniq]

        query_num = len(query_seq_id_all)

        fig_width_ABCD = 9

        fig_height_A = round(len(As_list_uniq) / 2)
        fig_height_B = round(len(Bs_list_uniq) / 4)
        fig_height_C = round(len(Cs_list_uniq) / 3)

        barh_plotter(ko_num_list_A, As_list_uniq_desc, query_num, query_ko_NA, fig_width_ABCD, fig_height_A, stats_plot_A)
        barh_plotter(ko_num_list_B, Bs_list_uniq_desc, query_num, query_ko_NA, fig_width_ABCD, fig_height_B, stats_plot_B)
        barh_plotter(ko_num_list_C, Cs_list_uniq_desc, query_num, query_ko_NA, fig_width_ABCD, fig_height_C, stats_plot_C)


def Annotation_KEGG(args):

    input_file_faa =      args['seq_in']
    input_file_user_ko =  args['ko_in']
    file_extension =      args['x']
    KEGG_DB_folder =      args['db_dir']
    run_diamond =         args['diamond']
    num_threads =         args['t']
    plot_stats =          args['plot']

    time_format = '[%Y-%m-%d %H:%M:%S] '

    run_blast = None
    if (input_file_faa is not None) and (input_file_user_ko is None):
        run_blast = True
    elif (input_file_faa is None) and (input_file_user_ko is not None):
        run_blast = False
    else:
        print(datetime.now().strftime(time_format) + 'Please provide input file with either "-seq_in" or "-ko_in", do not provide the two options at the same time')
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
        print(datetime.now().strftime(time_format) + 'Input sequence file detected, will run blastp/diamond against KEGG database at first')
        sleep(0.5)
    else:
        print(datetime.now().strftime(time_format) + 'Annotation results provided, blastp/diamond skipped')
        sleep(0.5)


    ################################################### define file name ###################################################

    KEGG_DB_seq =         '%s/prokaryotes.pep.fasta'      % KEGG_DB_folder
    KEGG_DB_seq_diamond = '%s/prokaryotes.pep.fasta.dmnd' % KEGG_DB_folder
    KEGG_DB_seq2ko =      '%s/prokaryotes.dat'            % KEGG_DB_folder
    KEGG_DB_ko =          '%s/ko00001.keg'                % KEGG_DB_folder


    ############################################ check whether diamond db exist ############################################

    if (run_blast is True) and (run_diamond is True):
        if os.path.isfile(KEGG_DB_seq_diamond) is False:
            print(datetime.now().strftime(time_format) + 'DB file not found, making diamond db with %s' % KEGG_DB_seq)

            if os.path.isfile(KEGG_DB_seq) is True:
                diamond_makedb_cmd = 'diamond makedb --in %s --db %s --quiet' % (KEGG_DB_seq, KEGG_DB_seq_diamond)
                os.system(diamond_makedb_cmd)
            else:
                print(datetime.now().strftime(time_format) + '%s not found, program exited' % KEGG_DB_seq)
                exit()


    ################################################ Read in KEGG DB files #################################################

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

    # get D2A, D2B and D2C dict
    D2A_dict = {}
    D2B_dict = {}
    D2C_dict = {}
    for each_D in D2ABCD_dict:
        ABCD_list = D2ABCD_dict[each_D]
        Current_As = []
        Current_Bs = []
        Current_Cs = []
        for each_ABCD in ABCD_list:
            each_ABCD_split = each_ABCD.split('|')

            if each_ABCD_split[0] not in Current_As:
                Current_As.append(each_ABCD_split[0])

            if each_ABCD_split[1] not in Current_Bs:
                Current_Bs.append(each_ABCD_split[1])

            if each_ABCD_split[2] not in Current_Cs:
                Current_Cs.append(each_ABCD_split[2])

        D2A_dict[each_D] = Current_As
        D2B_dict[each_D] = Current_Bs
        D2C_dict[each_D] = Current_Cs

    # get db_seq_to_KO_dict
    db_seq_to_KO_dict = {}
    for each_hit in open(KEGG_DB_seq2ko):
        each_hit_split = each_hit.strip().split('\t')
        db_seq = each_hit_split[0]
        hit_id_KO = each_hit_split[1]
        if hit_id_KO != '':
            db_seq_to_KO_dict[db_seq] = hit_id_KO


    ########################################################################################################################

    # check whether the input file is a file or fplder
    if os.path.isfile(input_file_folder) is True:

        # create output folder
        input_file_path, input_file_basename, input_file_ext = sep_path_basename_ext(input_file_folder)

        Annotation_KEGG_worker([input_file_folder,
                                run_blast,
                                run_diamond,
                                KEGG_DB_seq,
                                KEGG_DB_seq_diamond,
                                As_description_dict,
                                Bs_description_dict,
                                Cs_description_dict,
                                Ds_description_dict,
                                D2ABCD_dict,
                                db_seq_to_KO_dict,
                                input_file_path,
                                plot_stats])


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
        else:
            print(datetime.now().strftime(time_format) + 'Running KEGG annotation for %s input files with %s cores' % (len(input_file_name_list), num_threads))

        list_for_multiple_arguments_KEGG = []
        for input_file in input_file_name_list:

            pwd_input_file = '%s/%s' % (input_file_folder, input_file)

            list_for_multiple_arguments_KEGG.append([pwd_input_file,
                                                     run_blast,
                                                     run_diamond,
                                                     KEGG_DB_seq,
                                                     KEGG_DB_seq_diamond,
                                                     As_description_dict,
                                                     Bs_description_dict,
                                                     Cs_description_dict,
                                                     Ds_description_dict,
                                                     D2ABCD_dict,
                                                     db_seq_to_KO_dict,
                                                     output_folder,
                                                     plot_stats])

        # run KEGG annotaion files with multiprocessing
        pool = mp.Pool(processes=num_threads)
        pool.map(Annotation_KEGG_worker, list_for_multiple_arguments_KEGG)
        pool.close()
        pool.join()


    ################################################## Final report ####################################################

    print(datetime.now().strftime(time_format) + 'Done!')


# if __name__ == "__main__":
#
#     parser = argparse.ArgumentParser()
#
#     parser.add_argument('-seq_in',  required=False, help='faa file')
#     parser.add_argument('-ko_in',   required=False, help='annotation results from BlastKOALA/GhostKOALA, normally with name user_ko.txt')
#     parser.add_argument('-x',       required=True,  help='file extension')
#     parser.add_argument('-p',       required=True,  help='output prefix')
#     parser.add_argument('-db_dir',  required=True,  help='folder holds prokaryotes.pep.fasta, prokaryotes.dat and ko00001.keg')
#     parser.add_argument('-diamond', required=False, action='store_true', help='run diamond (for big dataset), default is NCBI blastp')
#     parser.add_argument('-t',       required=False, default=1, type=int, help='number of threads, default: 1')
#
#     args = vars(parser.parse_args())
#
#     Annotation_KEGG(args)
