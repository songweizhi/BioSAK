import argparse


def get_reads_from_sam(args):

    sam_file    = args['sam']
    output_file = args['out']
    ctgs_file   = args['ctgs']
    option      = args['option']

    # get contig list
    ctg_list = []
    for each_ctg in open(ctgs_file):
        each_ctg = each_ctg.strip()
        ctg_list.append(each_ctg)

    # export reads
    output = open(output_file, 'w')
    for each in open(sam_file):
        if not each.startswith('@'):
            each_split = each.strip().split('\t')
            query_name = each_split[0]
            ref_name = each_split[2]
            query_seq = each_split[9]

            if option == 1:
                if ref_name in ctg_list:
                    output.write('>%s\n' % query_name)
                    output.write('%s\n' % query_seq)

            if option == 0:
                if ref_name not in ctg_list:
                    output.write('>%s\n' % query_name)
                    output.write('%s\n' % query_seq)

    output.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-sam',     required=True,              help='Input sam file')
    parser.add_argument('-ctgs',    required=True,              help='Contig list')
    parser.add_argument('-option',  required=True, type=int,    help="Specify '1' to get reads mapped to provided contigs, '0' to get reads not mapped to provided contigs")
    parser.add_argument('-out',     required=True,              help='Output fasta file')
    args = vars(parser.parse_args())
    get_reads_from_sam(args)
