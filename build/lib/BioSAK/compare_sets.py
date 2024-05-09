import argparse


compare_sets_usage = '''
========= compare_sets example command =========

BioSAK compare_sets -1 file1.txt -2 file2.txt

================================================
'''


def compare_sets(args):

    file_1 = args['1']
    file_2 = args['2']

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
        print('Shared IDs (%s)\n%s\n' % (len(shared_id_set), '\n'.join(shared_id_set)))
        print('IDs uniq to %s (%s)\n%s\n' % (file_1, len(set_1_uniq), '\n'.join(sorted(list(set_1_uniq)))))
        print('IDs uniq to %s (%s)\n%s\n' % (file_2, len(set_2_uniq), '\n'.join(sorted(list(set_2_uniq)))))


if __name__ == '__main__':

    compare_sets_parser = argparse.ArgumentParser(usage=compare_sets_usage)
    compare_sets_parser.add_argument('-1', required=True, help='file 1')
    compare_sets_parser.add_argument('-2', required=True, help='file 2')
    args = vars(compare_sets_parser.parse_args())
    compare_sets(args)
