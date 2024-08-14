import os
import argparse


stats_arcog_usage = '''
========================== stats_arcog example commands ==========================

BioSAK stats_arcog -cog cog.txt -o op_dir -db /Users/songweizhi/DB/arCOG18

# Required DB files:
arCOGdef.tab: https://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/tmp.ar18/arCOGdef.tab
fun-20.tab:   https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab

==================================================================================
'''


def stats_arcog(args):

    cog_txt             = args['cog']
    op_dir              = args['o']
    db_dir              = args['db']
    force_create_op_dir = args['f']

    ####################################################################################################################

    cog_des_txt     = '%s/arCOGdef.tab'         % db_dir
    pwd_fun_20_tab  = '%s/fun-20.tab'           % db_dir
    op_txt          = '%s/COG_cate.txt'         % op_dir
    op_stats_txt    = '%s/COG_cate_stats.txt'   % op_dir

    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    ############################################## Read in KEGG DB files ###############################################

    print('Read in arCOG DB files')

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

    cog_category_to_description_dict = {}
    for cog_category in open(pwd_fun_20_tab):
        cog_category_split = cog_category.strip().split('\t')
        cog_category_to_description_dict[cog_category_split[0]] = cog_category_split[2]

    ########################################################################################################################

    input_cog_id_set = set()
    for each_cog in open(cog_txt):
        cog_id = each_cog.strip().split('\t')[0]
        if cog_id.startswith('arCOG'):
            input_cog_id_set.add(cog_id)

    cog_cate_stats_dict = dict()
    for cog_id in input_cog_id_set:
        cog_cate_list = cog_id_to_category_dict[cog_id]
        for cog_cate in cog_cate_list:
            if cog_cate not in cog_cate_stats_dict:
                cog_cate_stats_dict[cog_cate] = set()
            cog_cate_stats_dict[cog_cate].add(cog_id)

    op_txt_handle = open(op_txt, 'w')
    op_txt_handle.write('Category\tCategory_description\tCOG\tDescription\n')
    op_stats_txt_handle = open(op_stats_txt, 'w')
    for cog_cate in sorted(list(cog_cate_stats_dict.keys())):
        cog_id_set = cog_cate_stats_dict[cog_cate]
        cog_cate_desc = cog_category_to_description_dict[cog_cate]
        op_stats_txt_handle.write('%s\t%s\t%s\n' % (cog_cate, len(cog_id_set), cog_cate_desc))
        for cog_id in sorted(list(cog_id_set)):
            cog_id_desc = cog_id_to_description_dict[cog_id]
            op_txt_handle.write('%s\t%s\t%s\t%s\t%s\n' % (cog_cate, len(cog_id_set), cog_cate_desc.strip(), cog_id, cog_id_desc.strip()))
    op_txt_handle.close()
    op_stats_txt_handle.close()

    print('Done!')


if __name__ == "__main__":

    stats_arcog_parser = argparse.ArgumentParser(usage=stats_arcog_usage)
    stats_arcog_parser.add_argument('-cog',         required=True,                          help='cog_ids.txt')
    stats_arcog_parser.add_argument('-db',          required=True,                          help='DB folder, contains arCOGdef.tab and fun-20.tab')
    stats_arcog_parser.add_argument('-o',           required=True,                          help='output directory')
    stats_arcog_parser.add_argument('-f',           required=False, action="store_true",    help='force overwrite')
    args = vars(stats_arcog_parser.parse_args())
    stats_arcog(args)
