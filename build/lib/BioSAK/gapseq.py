import os
import glob
import argparse


gapseq_usage = '''
============ gapseq example commands ============

BioSAK gapseq -i gapseq_op_dir -o df.txt
BioSAK gapseq -i gapseq_op_dir -o df.txt -name

=================================================
'''


def sep_path_basename_ext(file_in):

    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)
    return f_path, f_base, f_ext


def gapseq(args):

    dir_in       = args['i']
    df_out       = args['o']
    include_desc = args['name']

    file_re   = '%s/*-Pathways.tbl' % dir_in
    file_list = glob.glob(file_re)

    # read in annotation results
    gnm_id_set = set()
    detected_pwy_set = set()
    pwy_to_gnm_dict = dict()
    pwy_id_to_name_dict = dict()
    annotation_results_dict = dict()
    for each_file in sorted(file_list):
        _, f_base, _ = sep_path_basename_ext(each_file)
        gnm_id = '-'.join(f_base.split('-')[:-2])
        gnm_id_set.add(gnm_id)
        encoded_pwy_set = set()
        for each_line in open(each_file):
            if not each_line.startswith('#'):
                if not each_line.startswith('ID\tName'):
                    each_line_split = each_line.strip().split('\t')
                    pwy_id = each_line_split[0]
                    if (pwy_id[0] == '|') and (pwy_id[-1] == '|'):
                        pwy_id = pwy_id[1:-1]
                    pwy_name = each_line_split[1]
                    pwy_id_to_name_dict[pwy_id] = pwy_name
                    pa = each_line_split[2]
                    if pa == 'true':
                        detected_pwy_set.add(pwy_id)
                        encoded_pwy_set.add(pwy_id)
                        if pwy_id not in pwy_to_gnm_dict:
                            pwy_to_gnm_dict[pwy_id] = set()
                        pwy_to_gnm_dict[pwy_id].add(gnm_id)
        annotation_results_dict[gnm_id] = encoded_pwy_set

    pwy_id_list_sorted = sorted(list(detected_pwy_set))
    pwy_id_list_sorted_desc = [('%s__%s' % (i, pwy_id_to_name_dict[i])) for i in pwy_id_list_sorted]

    df_out_handle = open(df_out, 'w')
    if include_desc is False:
        df_out_handle.write('\t'.join(pwy_id_list_sorted) + '\n')
    else:
        df_out_handle.write('\t'.join(pwy_id_list_sorted_desc) + '\n')
    for each_gnm in sorted(list(gnm_id_set)):
        encoded_pwys = annotation_results_dict[each_gnm]
        value_list = [each_gnm]
        for each_pwy in pwy_id_list_sorted:
            if each_pwy in encoded_pwys:
                value_list.append('1')
            else:
                value_list.append('0')
        df_out_handle.write('\t'.join(value_list) + '\n')
    df_out_handle.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',       required=True,                        help='folder holds the *Pathways.tbl files')
    parser.add_argument('-o',       required=True,                        help='output file')
    parser.add_argument('-name',    required=False, action="store_true",  help='include pathway name in the output dataframe')
    args = vars(parser.parse_args())
    gapseq(args)
