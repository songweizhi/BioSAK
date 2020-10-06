import argparse
from Bio import SeqIO
from BioSAK.global_functions import sep_path_basename_ext


rename_seq_usage = '''
========================= rename_seq example commands =========================

# rename "NODE_941_length_17600_cov_52.7123" to "NODE_941"
BioSAK rename_seq -in Contigs.fa -sep_in "_" -n 2

# rename "Seawater|NODE|941|length|17600|cov|52.7123" to "Seawater_NODE_941"
BioSAK rename_seq -in Contigs.fa -sep_in "|" -sep_out "_" -n 3

# add prefix to all sequences in a fasta file
BioSAK rename_seq -in Contigs.fa -prefix seawater

# rename "NODE_941_length_17600_cov_52.7123" to "Seawater_NODE_941"
BioSAK rename_seq -in Contigs.fa -sep_in "_" -n 2 -prefix Seawater

===============================================================================
'''


def rename_seq(args):

    ctg_file_in =    args['in']
    sep_in =         args['sep_in']
    sep_out =        args['sep_out']
    column_to_keep = args['n']
    add_prefix =     args['prefix']
    one_line =       args['oneline']

    ctg_file_path, ctg_file_basename, ctg_file_ext = sep_path_basename_ext(ctg_file_in)
    ctg_file_out = '%s/%s_renamed%s' % (ctg_file_path, ctg_file_basename, ctg_file_ext)

    if sep_out is None:
        sep_out = sep_in

    ctg_file_out_handle = open(ctg_file_out, 'w')
    for Seq_record in SeqIO.parse(ctg_file_in, 'fasta'):

        if (sep_in is not None) and (add_prefix is None):
            Seq_record_id_new = sep_out.join(Seq_record.id.split(sep_in)[:column_to_keep])

        elif (sep_in is None) and (add_prefix is not None):
            Seq_record_id_new = '%s_%s' % (add_prefix, Seq_record.id)

        elif (sep_in is not None) and (add_prefix is not None):
            Seq_record_id_new = '%s_%s' % (add_prefix, sep_out.join(Seq_record.id.split(sep_in)[:column_to_keep]))

        else:
            Seq_record_id_new = ''
            print('Don\'t know what to do, program exited!')
            exit()

        if one_line is True:
            ctg_file_out_handle.write('>%s\n' % Seq_record_id_new)
            ctg_file_out_handle.write('%s\n' % Seq_record.seq)
        else:
            Seq_record.id = Seq_record_id_new
            SeqIO.write(Seq_record, ctg_file_out_handle, 'fasta')

    ctg_file_out_handle.close()


if __name__ == '__main__':

    rename_seq_parser = argparse.ArgumentParser()

    # arguments for rename_seq_parser
    rename_seq_parser.add_argument('-in',         required=True,                          help='input sequence file')
    rename_seq_parser.add_argument('-sep_in',     required=False, default=None,           help='separator for input sequences')
    rename_seq_parser.add_argument('-sep_out',    required=False, default=None,           help='separator for output sequences, default: same as sep_in')
    rename_seq_parser.add_argument('-n',          required=False, default=None, type=int, help='the number of columns to keep')
    rename_seq_parser.add_argument('-prefix',     required=False, default=None,           help='add prefix to sequence')
    rename_seq_parser.add_argument('-oneline',    required=False, action="store_true",    help='put sequence in a single line')

    args = vars(rename_seq_parser.parse_args())

    rename_seq(args)
