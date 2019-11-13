import argparse
from Bio import SeqIO
from BioSAK.global_functions import sep_path_basename_ext


COG2014_parser_usage = '''
========================= rename_ctg example commands =========================

# rename "NODE_941_length_17600_cov_52.7123" to "NODE_941"
BioSAK rename_ctg i Contigs.fa -si "_" -n 2

# rename "Seawater|NODE|941|length|17600|cov|52.7123" to "Seawater_NODE_941"
BioSAK rename_ctg i Contigs.fa -si "|" -so "_" -n 3

===============================================================================
'''


def rename_ctg(args):

    ctg_file_in = args['i']
    sep_in = args['si']
    sep_out = args['so']
    column_to_keep = args['n']

    if sep_out is None:
        sep_out = sep_in

    ctg_file_path, ctg_file_basename, ctg_file_ext = sep_path_basename_ext(ctg_file_in)
    ctg_file_out = '%s/%s_renamed%s' % (ctg_file_path, ctg_file_basename, ctg_file_ext)

    ctg_file_out_handle = open(ctg_file_out, 'w')
    for Seq_record in SeqIO.parse(ctg_file_in, 'fasta'):
        Seq_record_id_new = sep_out.join(Seq_record.id.split(sep_in)[:column_to_keep])
        Seq_record.id = Seq_record_id_new
        SeqIO.write(Seq_record, ctg_file_out_handle, 'fasta')
    ctg_file_out_handle.close()


if __name__ == '__main__':

    rename_ctg_parser = argparse.ArgumentParser()

    # arguments for rename_ctg_parser
    rename_ctg_parser.add_argument('-i',    required=True,                  help='input sequence file')
    rename_ctg_parser.add_argument('-si',   required=True,                  help='separator for input sequences')
    rename_ctg_parser.add_argument('-so',   required=False, default=None,   help='separator for output sequences, default: same as si')
    rename_ctg_parser.add_argument('-n',    required=True, type=int,        help='the number of columns to keep')

    args = vars(rename_ctg_parser.parse_args())

    rename_ctg(args)
