import os
import glob
import argparse


Combine_KEGG_arCOG_usage = '''
================= Combine_KEGG_arCOG example commands =================

BioSAK Combine_KEGG_arCOG -kegg KEGG_wd -arcog arCOG_wd -o op_dir -f

=======================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def Combine_KEGG_arCOG(args):

    output_folder       = args['o']
    kegg_annotation_wd  = args['kegg']
    arcog_annotation_wd = args['arcog']
    force_create_op_dir = args['f']

    ############################################## define output file name #############################################

    fun_assignment_dir  = '%s/fun_assignment'   % output_folder
    stats_dir_num       = '%s/fun_stats_num'    % output_folder
    stats_dir_pct       = '%s/fun_stats_pct'    % output_folder

    # create output folder
    if os.path.isdir(output_folder) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % output_folder)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % output_folder)
    os.system('mkdir %s' % fun_assignment_dir)
    os.system('mkdir %s' % stats_dir_num)
    os.system('mkdir %s' % stats_dir_pct)

    ################################################### get file list ##################################################

    arcog_assignment_re         = '%s/*_arCOG_wd/*_query_to_cog.txt'        %   arcog_annotation_wd
    ko_assignment_re            = '%s/*_KEGG_wd/*_ko_assignment_ABCD.txt'   %   kegg_annotation_wd

    arcog_assignment_file_list  = glob.glob(arcog_assignment_re)
    ko_assignment_file_list     = glob.glob(ko_assignment_re)

    ######################################### read in kegg annotation results ##########################################

    fun_id_to_description_dict = dict()
    fun_id_to_B_dict = dict()
    fun_id_to_C_dict = dict()
    gnm_to_gene_dict = dict()
    ko_assignment_dod = dict()
    for each_file in ko_assignment_file_list:
        f_name, f_path, f_base, f_ext = sep_path_basename_ext(each_file)
        gnm_id = f_base[:-len('_ko_assignment_ABCD')]

        if gnm_id not in gnm_to_gene_dict:
            gnm_to_gene_dict[gnm_id] = set()

        current_genome_ko_assignment_dict = dict()
        for each_line in open(each_file):
            if not each_line.startswith('Gene_id\t'):
                each_line_split = each_line.strip().split('\t')
                gene_id = each_line_split[0]
                gnm_to_gene_dict[gnm_id].add(gene_id)
                if len(each_line_split) > 1:
                    ko_id       = each_line_split[4][2:]
                    ko_b        = each_line_split[2]
                    ko_b_desc   = each_line_split[6]
                    ko_c        = each_line_split[3]
                    ko_c_desc   = each_line_split[7]
                    ko_desc     = each_line_split[8].strip()
                    ko_b_with_desc = '%s__%s' % (ko_b, ko_b_desc)
                    ko_c_with_desc = '%s__%s' % (ko_c, ko_c_desc)
                    fun_id_to_B_dict[ko_id] = ko_b_with_desc
                    fun_id_to_C_dict[ko_id] = ko_c_with_desc
                    fun_id_to_description_dict[ko_id] = ko_desc
                    current_genome_ko_assignment_dict[gene_id] = ko_id
        ko_assignment_dod[gnm_id] = current_genome_ko_assignment_dict

    ######################################### read in arcog annotation results #########################################

    arcog_assignment_dod = dict()
    arcog_id_to_cate_dict = dict()
    arcog_cate_to_desc_dict = dict()
    for each_file in arcog_assignment_file_list:

        f_name, f_path, f_base, f_ext = sep_path_basename_ext(each_file)
        gnm_id = f_base[:-len('_query_to_cog')]

        if gnm_id not in gnm_to_gene_dict:
            gnm_to_gene_dict[gnm_id] = set()

        current_genome_arcog_assignment_dict = dict()
        for each_line in open(each_file):
            if not each_line.startswith('Query\tarCOG'):
                each_line_split = each_line.strip().split('\t')
                gene_id = each_line_split[0]
                gnm_to_gene_dict[gnm_id].add(gene_id)
                if len(each_line_split) > 1:
                    arcog_id = each_line_split[1]
                    arcog_cates = each_line_split[2]
                    arcog_desc = each_line_split[3]

                    if arcog_id not in arcog_id_to_cate_dict:
                        arcog_id_to_cate_dict[arcog_id] = set()

                    for arcog_cate in arcog_cates:
                        arcog_id_to_cate_dict[arcog_id].add(arcog_cate)

                    fun_id_to_description_dict[arcog_id] = arcog_desc
                    current_genome_arcog_assignment_dict[gene_id] = arcog_id
        arcog_assignment_dod[gnm_id] = current_genome_arcog_assignment_dict

        # read in cog cate to desc info
        func_stats_txt = each_file.replace('_query_to_cog.txt', '_func_stats.txt')
        for each in open(func_stats_txt):
            each_split = each.strip().split('\t')
            cate_id = each_split[0]
            cate_desc = each_split[2]
            arcog_cate_to_desc_dict[cate_id] = cate_desc

    ###################################### write out combined annotation results #######################################

    # write out combined annotation results
    combined_annotation_dod = dict()
    for each_gnm in sorted(list(gnm_to_gene_dict.keys())):
        gnm_arcog_dict = arcog_assignment_dod[each_gnm]
        gnm_kegg_dict  = ko_assignment_dod[each_gnm]
        gnm_gene_list  = sorted(list(gnm_to_gene_dict[each_gnm]))
        combined_annotation_dict = dict()
        combined_annotation_txt = '%s/%s.txt' % (fun_assignment_dir, each_gnm)
        combined_annotation_txt_handle = open(combined_annotation_txt, 'w')
        for each_gene in gnm_gene_list:
            if each_gene in gnm_kegg_dict:
                fun_assignment = gnm_kegg_dict[each_gene]
                fun_desc = fun_id_to_description_dict[fun_assignment]
                combined_annotation_txt_handle.write('%s\t%s\t%s\n' % (each_gene, fun_assignment, fun_desc))
                combined_annotation_dict[each_gene] = gnm_kegg_dict[each_gene]
            elif each_gene in gnm_arcog_dict:
                fun_assignment = gnm_arcog_dict[each_gene]
                fun_desc = fun_id_to_description_dict[fun_assignment]
                combined_annotation_txt_handle.write('%s\t%s\t%s\n' % (each_gene, fun_assignment, fun_desc))
                combined_annotation_dict[each_gene] = gnm_arcog_dict[each_gene]
            else:
                combined_annotation_txt_handle.write(each_gene + '\n')
        combined_annotation_txt_handle.close()

        combined_annotation_dod[each_gnm] = combined_annotation_dict

    ############################################ write out annotation stats ############################################

    for each_gnm in sorted(list(gnm_to_gene_dict.keys())):
        annotation_dict = combined_annotation_dod[each_gnm]

        fun_stats_num_txt = '%s/%s.txt' % (stats_dir_num, each_gnm)
        fun_stats_pct_txt = '%s/%s.txt' % (stats_dir_pct, each_gnm)

        stats_dict = dict()
        gene_num = 0
        for each_gene in annotation_dict:
            fun_id = annotation_dict[each_gene]
            if fun_id not in stats_dict:
                stats_dict[fun_id] = 1
            else:
                stats_dict[fun_id] += 1
            gene_num += 1

        stats_dict_pct = dict()
        for each_fun in stats_dict:
            fun_num = stats_dict[each_fun]
            fun_pct = fun_num*100/gene_num
            fun_pct = float("{0:.2f}".format(fun_pct))
            stats_dict_pct[each_fun] = fun_pct

        fun_stats_num_txt_handle = open(fun_stats_num_txt, 'w')
        fun_stats_pct_txt_handle = open(fun_stats_pct_txt, 'w')
        for each_fun in sorted(list(stats_dict.keys())):
            fun_num = stats_dict[each_fun]
            fun_pct = stats_dict_pct[each_fun]
            fun_desc = fun_id_to_description_dict[each_fun]

            fun_note = ''
            if each_fun.startswith('arCOG'):
                fun_cate_list = sorted(list(arcog_id_to_cate_dict[each_fun]))
                for fun_cate in fun_cate_list:
                    fun_cate_desc = arcog_cate_to_desc_dict[fun_cate]
                    fun_note += ('%s__%s' % (fun_cate, fun_cate_desc))
            elif each_fun.startswith('K'):
                print(each_fun)

                fun_b = fun_id_to_B_dict[each_fun]
                fun_c = fun_id_to_C_dict[each_fun]
                fun_note = '%s;%s' % (fun_c, fun_b)


            fun_stats_num_txt_handle.write('%s\t%s\t%s\t%s\n' % (each_fun, fun_num, fun_desc, fun_note))
            fun_stats_pct_txt_handle.write('%s\t%s\t%s\t%s\n' % (each_fun, fun_pct, fun_desc, fun_note))
        fun_stats_num_txt_handle.close()
        fun_stats_pct_txt_handle.close()

    ####################################################################################################################


if __name__ == '__main__':

    Combine_KEGG_arCOG_parser = argparse.ArgumentParser(usage=Combine_KEGG_arCOG_usage)
    Combine_KEGG_arCOG_parser.add_argument('-o',        required=True,                          help='output directory')
    Combine_KEGG_arCOG_parser.add_argument('-kegg',     required=True,                          help='annotation results from KEGG module')
    Combine_KEGG_arCOG_parser.add_argument('-arcog',    required=True,                          help='annotation results from arCOG module')
    Combine_KEGG_arCOG_parser.add_argument('-f',        required=False, action="store_true",    help='force overwrite')
    args = vars(Combine_KEGG_arCOG_parser.parse_args())
    Combine_KEGG_arCOG(args)


'''

cd /Users/songweizhi/Desktop
python3 /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/Combine_KEGG_arCOG.py -kegg 3_combined_genomes_50_5_dRep97_284_KEGG_wd -arcog 3_combined_genomes_50_5_dRep97_284_arCOG_wd_0.0001 -o 3_combined_genomes_50_5_dRep97_284_combined_KEGG_arCOG -f

'''
