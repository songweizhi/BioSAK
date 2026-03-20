import os
import argparse


rstudio_color_usage = '''
======================= rstudio_color example commands =======================

BioSAK rstudio_color -i color_code.txt -o color_code.R
BioSAK rstudio_color -i color_code.txt -o color_code.R -t WoRMS_taxonomy.txt

==============================================================================
'''


def rstudio_color(args):

    color_code_txt  = args['i']
    taxonomy_txt    = args['t']
    fixed_tax_len   = args['l']
    rstudio_txt     = args['o']

    if os.path.isfile(color_code_txt) is False:
        print('%s does not exist, program exited!' % color_code_txt)
        exit()

    if taxonomy_txt is None:
        rstudio_txt_handle = open(rstudio_txt, 'w')
        for each_line in open(color_code_txt):
            each_line_split = each_line.strip().split('\t')
            if len(each_line_split) == 2:
                rstudio_txt_handle.write('"%s"\t%s\n' % (each_line_split[1], each_line_split[0]))
        rstudio_txt_handle.close()
    else:
        if os.path.isfile(taxonomy_txt) is False:
            print('%s does not exist, program exited!' % taxonomy_txt)
            exit()

        color_code_dict = dict()
        for each_line in open(color_code_txt):
            each_line_split = each_line.strip().split('\t')
            color_code_dict[each_line_split[0]] = each_line_split[1]

        rstudio_txt_handle = open(rstudio_txt, 'w')
        for each_line in open(taxonomy_txt):
            each_line_split = each_line.strip().split(';')
            str_for_r = ''
            for each_r in each_line_split:
                each_r_no_space = each_r
                each_r_no_space = each_r_no_space.replace(' ', '_')
                current_color_code = color_code_dict.get(each_r, '')
                space_to_add = ' '*(fixed_tax_len-len(each_r))
                if current_color_code == '':
                    space_to_add = space_to_add + '       '
                each_r_fixed_len = '%s%s' % (each_r, ' '*(fixed_tax_len-len(each_r)))
                str_for_r += '%s("%s")%s' % (each_r_no_space, current_color_code, space_to_add)
            rstudio_txt_handle.write(str_for_r + '\n')
        rstudio_txt_handle.close()


if __name__ == '__main__':

    rstudio_color_parser = argparse.ArgumentParser(usage=rstudio_color_usage)
    rstudio_color_parser.add_argument('-i', required=True,                          help='color code txt')
    rstudio_color_parser.add_argument('-t', required=False, default=None,           help='taxonomy, in GTDB format')
    rstudio_color_parser.add_argument('-l', required=False, default=25, type=int,   help='fix length for each label column, default is 25')
    rstudio_color_parser.add_argument('-o', required=True,                          help='output R file')
    args = vars(rstudio_color_parser.parse_args())
    rstudio_color(args)

