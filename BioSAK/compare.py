import os
import argparse


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


compare_usage = '''
========= compare example command =========

BioSAK compare -1 file1.txt -2 file2.txt

===========================================
'''


def compare(args):

    file_1 = args['1']
    file_2 = args['2']

    _, _, f1_base, _ = sep_path_basename_ext(file_1)
    _, _, f2_base, _ = sep_path_basename_ext(file_2)

    file_1_uniq    = '%s_uniq.txt'      % f1_base
    file_2_uniq    = '%s_uniq.txt'      % f2_base
    file_12_shared = '%s_%s_shared.txt' % (f1_base, f2_base)

    set_1 = set()
    for each_line in open(file_1):
        set_1.add(each_line.strip().split()[0])

    set_2 = set()
    for each_line in open(file_2):
        set_2.add(each_line.strip().split()[0])

    shared_id_set = set(set_1).intersection(set_2)

    set_1_uniq = set()
    for e1 in set_1:
        if e1 not in shared_id_set:
            set_1_uniq.add(e1)

    set_2_uniq = set()
    for e2 in set_2:
        if e2 not in shared_id_set:
            set_2_uniq.add(e2)

    if len(shared_id_set) == 0:
        print('No id shared between %s and %s' % (file_1, file_2))
    else:
        if len(set_1_uniq) > 0:
            file_1_uniq_handle = open(file_1_uniq, 'w')
            file_1_uniq_handle.write('\n'.join(sorted(list(set_1_uniq))) + '\n')
            file_1_uniq_handle.close()
        if len(set_2_uniq) > 0:
            file_2_uniq_handle = open(file_2_uniq, 'w')
            file_2_uniq_handle.write('\n'.join(sorted(list(set_2_uniq))) + '\n')
            file_2_uniq_handle.close()
        if len(shared_id_set) > 0:
            file_12_shared_uniq_handle = open(file_12_shared, 'w')
            file_12_shared_uniq_handle.write('\n'.join(sorted(list(shared_id_set))) + '\n')
            file_12_shared_uniq_handle.close()


if __name__ == '__main__':

    compare_parser = argparse.ArgumentParser(usage=compare_usage)
    compare_parser.add_argument('-1', required=True, help='file 1')
    compare_parser.add_argument('-2', required=True, help='file 2')
    args = vars(compare_parser.parse_args())
    compare(args)
