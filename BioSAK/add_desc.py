import os
import argparse


add_desc_usage = '''
====================== add_desc example commands ======================

BioSAK add_desc -i arCOG_df.txt -o arCOG_df_desc.txt -d arCOGdef.tab
BioSAK add_desc -i KEGG_df.txt -o KEGG_df_desc.txt -d ko00001.keg

=======================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)
    return f_path, f_base, f_ext


def add_desc(args):

    df_in   = args['i']
    df_out  = args['o']
    db_file = args['d']

    f_path, f_base, f_ext = sep_path_basename_ext(db_file)

    # read in cog_des_txt
    description_dict = dict()
    if f_base == 'arCOGdef':
        for each_cog in open(db_file, encoding="ISO-8859-1"):
            each_cog_split = each_cog.strip().split('\t')
            cog_id = each_cog_split[0]
            cog_desc = each_cog_split[3]
            description_dict[cog_id] = cog_desc

    elif f_base == 'ko00001':
        description_dict = {}
        for each_line in open(db_file):
            if each_line[0] in ['A', 'B', 'C', 'D']:
                each_line_split = each_line.strip().split(' ')
                if each_line[0] == 'D':
                    current_D_id = each_line_split[6]
                    current_D_description = ' '.join(each_line_split[8:])
                    description_dict[current_D_id] = current_D_description
    else:
        print('Please provide the correct db file, either arCOGdef.tab or ko00001.keg')
        print('Program exited!')
        exit()

    df_out_handle = open(df_out, 'w')
    line_index = 0
    for each_line in open(df_in):
        each_line_split = each_line.strip().split('\t')
        if line_index == 0:
            list_with_desc = ['']
            for each_id in each_line_split:
                list_with_desc.append('%s__%s' % (each_id, description_dict[each_id]))
            df_out_handle.write('\t'.join(list_with_desc) + '\n')
        else:
            df_out_handle.write(each_line)
        line_index += 1
    df_out_handle.close()


if __name__ == '__main__':

    add_desc_parser = argparse.ArgumentParser()
    add_desc_parser.add_argument('-i',  required=True,  help='input file ')
    add_desc_parser.add_argument('-o',  required=True,  help='output file')
    add_desc_parser.add_argument('-d',  required=True,  help='database file')
    args = vars(add_desc_parser.parse_args())
    add_desc(args)

