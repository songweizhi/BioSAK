import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


slice_seq_usage = '''
========================== slice_seq example commands ==========================

BioSAK slice_seq -in idba.fa -id ctg_01 -range 20-400 -out ctg1_20_400bp.fa
BioSAK slice_seq -in idba.fa -id ctg_01 -range 1-50 -out ctg1_1_50bp.fa -rc
BioSAK slice_seq -in idba.fa -id ctg_01 -range 30-end -out ctg1_30bp_to_end.fa

================================================================================
'''


def slice_seq(args):

    seq_file    = args['in']
    ctg_id      = args['id']
    seq_range   = args['range']
    get_rc      = args['rc']
    output_file = args['out']

    seq_range_split = seq_range.split('-')
    seq_range_l = int(seq_range_split[0])
    seq_range_r = seq_range_split[1]

    sequence_found = False
    output_file_handle = open(output_file, 'w')
    for seq_record in SeqIO.parse(seq_file, 'fasta'):
        if seq_record.id == ctg_id:

            sequence_found = True

            seq_subset = ''
            if seq_range_r == 'end':
                seq_subset = seq_record.seq[seq_range_l:]
            else:
                if int(seq_range_r) > len(seq_record.seq):
                    print('Specified end position longer than sequence, program exited!')
                    exit()
                else:
                    seq_subset = seq_record.seq[(seq_range_l - 1):int(seq_range_r)]

            if get_rc is True:
                seq_subset_rc = SeqRecord(seq_subset).reverse_complement().seq
                output_file_handle.write('>%s_%s_%s_rc\n' % (ctg_id, seq_range_l, seq_range_r))
                output_file_handle.write('%s\n' % seq_subset_rc)
            else:
                output_file_handle.write('>%s_%s_%s\n' % (ctg_id, seq_range_l, seq_range_r))
                output_file_handle.write('%s\n' % seq_subset)

    output_file_handle.close()

    if sequence_found is False:
        print('No sequence matched provided id, please double check!')


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()

    # arguments for select_seq
    parser.add_argument('-in',    required=True,                        help='sequence file')
    parser.add_argument('-id',    required=True,                        help='sequence id')
    parser.add_argument('-range', required=True,                        help='sequence range, start-end (in bp). e.g. 200-4000')
    parser.add_argument('-rc',    required=False, action='store_true',  help='write out reverse complement sequence')
    parser.add_argument('-out',   required=True,                        help='output file')

    args = vars(parser.parse_args())
    slice_seq(args)
