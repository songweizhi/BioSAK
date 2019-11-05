import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-sam',
                    required=True,
                    help='Input sam file')

parser.add_argument('-ctgs',
                    required=True,
                    help='Contig list')

parser.add_argument('-option',
                    required=True,
                    type=int,
                    help="'0' to export unmapped reads, '1' to export mapped reads")

parser.add_argument('-out',
                    required=True,
                    help='Output fasta file')

args = vars(parser.parse_args())
sam_file = args['sam']
output_file = args['out']
ctgs_file = args['ctgs']
option = args['option']


# get contig list
ctg_list = []
for each_ctg in open(ctgs_file):
    each_ctg = each_ctg.strip()
    ctg_list.append(each_ctg)


# export reads
output = open(output_file, 'w')
for each in open(sam_file):
    if each.startswith('@'):
        pass
    else:
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
