import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument('-sam',
                    required=True,
                    help='Input sam file')

parser.add_argument('-ctgs',
                    required=True,
                    help='Contig list')

parser.add_argument('-seq',
                    required=True,
                    help='Reference sequences')

args = vars(parser.parse_args())
sam_file = args['sam']
ctgs_file = args['ctgs']
seq_file = args['seq']


# get contig list
ctg_list = []
for each_ctg in open(ctgs_file):
    each_ctg = each_ctg.strip()
    ctg_list.append(each_ctg)


for each_ref in ctg_list:

    # get sequence length
    ref_len = 0
    for each_seq in SeqIO.parse(seq_file, 'fasta'):
        if each_seq.id == each_ref:
            ref_len = len(each_seq.seq)
            #print('reference lenth %s\t%s' % (each_ref, ref_len))

    # get total length of mapped reads
    total_aligned_length = 0
    for each_align in open(sam_file):
        if each_align.startswith('@'):
            pass
        else:
            each_split = each_align.strip().split('\t')
            query_name = each_split[0]
            ref_name = each_split[2]
            query_seq_len = len(each_split[9])
            if ref_name == each_ref:
                #print('reads length %s\t%s' % (each_ref, query_seq_len))
                total_aligned_length += query_seq_len
    #print('total_aligned_length %s' % total_aligned_length)


    # calculate average coverage of current contig
    mean_coverage = float("{0:.2f}".format(total_aligned_length/ref_len))

    print('%s\t%s' % (each_ref, mean_coverage))
