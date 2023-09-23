import os
import argparse


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


parser = argparse.ArgumentParser()
parser.add_argument('-i', required=False, type=int, help='input file, -o2')
args = vars(parser.parse_args())

txt_in      = args['i']

f_path, f_base, f_ext = sep_path_basename_ext(txt_in)

txt_out     = '%s/%s_krona.txt'     % (f_path, f_base)
krona_out   = '%s/%s_krona.html'    % (f_path, f_base)

tax_count_dict = dict()
for each in open(txt_in):
    each_split = each.strip().split(',')
    tax_str = '\t'.join(each_split[1:]).replace(' ', '_')
    print(tax_str)
    if tax_str not in tax_count_dict:
        tax_count_dict[tax_str] = 1
    else:
        tax_count_dict[tax_str] += 1

txt_out_handle = open(txt_out, 'w')
for each in tax_count_dict:
    txt_out_handle.write('%s\t%s\n' % (tax_count_dict[each], each))
txt_out_handle.close()

ktImportText_cmd = 'ktImportText %s -o %s' % (txt_out, krona_out)
os.system(ktImportText_cmd)
