import argparse


count_num_usage = '''
==== count_num example command ====

BioSAK count_num -i genus.txt

===================================
'''

def count_num(args):

    file_in = args['i']

    count_num_dict = dict()
    for each_line in open(file_in):
        each_line_split = each_line.strip().split()
        if each_line_split[0] not in count_num_dict:
            count_num_dict[each_line_split[0]] = 1
        else:
            count_num_dict[each_line_split[0]] += 1

    for each_key in sorted(list(count_num_dict.keys())):
        print('%s\t%s' % (each_key, count_num_dict[each_key]))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage=count_num_usage)
    parser.add_argument('-i', required=True, help='input file')
    args = vars(parser.parse_args())
    count_num(args)
