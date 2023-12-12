import argparse


ezaai2mat_usage = '''
============== ezaai2mat example command ==============

BioSAK ezaai2mat -i ezaai.table -o ezaai_matrix.txt

=======================================================
'''


def ezaai2mat(args):

    ezaai_table = args['i']
    op_matrix   = args['o']

    seq_set = set()
    aai_dict = dict()
    for each_line in open(ezaai_table):
        if not each_line.startswith('ID 1	ID 2'):
            each_line_split = each_line.strip().split('\t')
            seq_1 = each_line_split[2]
            seq_2 = each_line_split[3]
            aai_value = each_line_split[4]
            key_str = '%s__|__%s' % (seq_1, seq_2)
            aai_dict[key_str] = aai_value
            seq_set.add(seq_1)
            seq_set.add(seq_2)

    seq_list_sorted = sorted([i for i in seq_set])

    op_matrix_handle = open(op_matrix, 'w')
    op_matrix_handle.write('\t' + '\t'.join(seq_list_sorted) + '\n')
    for seq1 in seq_list_sorted:
        aai_list = [seq1]
        for seq2 in seq_list_sorted:
            key_str = '%s__|__%s' % (seq1, seq2)
            aai = aai_dict.get(key_str, 'na')
            aai_list.append(aai)
        op_matrix_handle.write('\t'.join(aai_list) + '\n')
    op_matrix_handle.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, help='ezaai produced table')
    parser.add_argument('-o', required=True, help='output AAI matrix')
    args = vars(parser.parse_args())
    ezaai2mat(args)
