import argparse


def get_reads_id_in_sam(args):

    sam_file     = args['i']
    reads_id_txt = args['o']

    read_id_set = set()
    for each_line in open(sam_file):
        if not each_line.startswith('@'):
            each_line_split = each_line.strip().split('\t')
            read_id = each_line_split[0]
            read_id_set.add(read_id)

    reads_id_txt_handle = open(reads_id_txt, 'w')
    for each_read in read_id_set:
        reads_id_txt_handle.write(each_read + '\n')
    reads_id_txt_handle.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, help='input sam file')
    parser.add_argument('-o', required=True, help='output txt file')
    args = vars(parser.parse_args())
    get_reads_id_in_sam(args)
