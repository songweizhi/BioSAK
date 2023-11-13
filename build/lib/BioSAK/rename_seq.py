import os
import glob
import argparse
from Bio import SeqIO
import multiprocessing as mp
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

# rename sequences in multiple files
BioSAK rename_seq -in seq_folder -x fa -sep_in "_" -n 2 -t 12
BioSAK rename_seq -in seq_folder -x fa -prefix prefix_file.txt

# prefix file format: file name (with extension) followed by prefix to be added 
# to the sequences in the file, tab separated.
MAG_1.fa mag1
MAG_2.fa mag2
MAG_3.fa mag3

===============================================================================
'''

def rename_seq_worker(argument_list):

    seq_file_in     = argument_list[0]
    seq_file_out    = argument_list[1]
    sep_in          = argument_list[2]
    sep_out         = argument_list[3]
    column_to_keep  = argument_list[4]
    add_prefix      = argument_list[5]
    one_line        = argument_list[6]

    if sep_out is None:
        sep_out = sep_in

    ctg_file_out_handle = open(seq_file_out, 'w')
    for Seq_record in SeqIO.parse(seq_file_in, 'fasta'):

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


def rename_seq(args):

    seq_file_in =    args['in']
    file_extension = args['x']
    sep_in =         args['sep_in']
    sep_out =        args['sep_out']
    column_to_keep = args['n']
    add_prefix =     args['prefix']
    one_line =       args['oneline']
    num_threads =    args['t']

    if os.path.isfile(seq_file_in) is True:
        ctg_file_path, ctg_file_basename, ctg_file_ext = sep_path_basename_ext(seq_file_in)
        seq_file_out = '%s/%s_renamed%s' % (ctg_file_path, ctg_file_basename, ctg_file_ext)
        if os.path.isfile(seq_file_out) is True:
            print('Output file detected, program exited: %s' % seq_file_out)
            exit()
        else:
            rename_seq_worker([seq_file_in, seq_file_out, sep_in, sep_out, column_to_keep, add_prefix, one_line])

    if os.path.isdir(seq_file_in) is True:

        if seq_file_in[-1] == '/':
            seq_file_in = seq_file_in[:-1]

        seq_in_folder_no_path = seq_file_in.split('/')[-1]
        seq_out_folder = '%s_renamed' % seq_in_folder_no_path

        seq_in_re   = '%s/*.%s' % (seq_file_in, file_extension)
        seq_in_list = [os.path.basename(file_name) for file_name in glob.glob(seq_in_re)]
        if len(seq_in_list) == 0:
            print('No sequence file detected, program exited!')
            exit()

        if os.path.isdir(seq_out_folder) is True:
            print('Output folder detected, program exited: %s' % seq_out_folder)
            exit()
        else:
            os.mkdir(seq_out_folder)

        prefix_dict = {}
        if add_prefix is not None:
            if os.path.isfile(add_prefix) is False:
                print('Prefix file not detected, program exited!')
                exit()
            else:
                # read in prefix
                for each_genome in open(add_prefix):
                    each_genome_split = each_genome.strip().split('\t')
                    prefix_dict[each_genome_split[0]] = each_genome_split[1]

                genome_without_prefix = set()
                for each_seq_file in seq_in_list:
                    if each_seq_file not in prefix_dict:
                        genome_without_prefix.add(each_seq_file)
                if len(genome_without_prefix) > 0:
                    print('Prefix for the following files not found, , program exited!')
                    print(','.join(genome_without_prefix))
                    exit()

        argument_lol_for_rename_seq_worker = []
        for each_seq_file in seq_in_list:
            pwd_seq_in  = '%s/%s' % (seq_file_in, each_seq_file)
            pwd_seq_out = '%s/%s' % (seq_out_folder, each_seq_file)
            current_argument_list = [pwd_seq_in, pwd_seq_out, sep_in, sep_out, column_to_keep, prefix_dict.get(each_seq_file, None), one_line]
            argument_lol_for_rename_seq_worker.append(current_argument_list)

        # rename sequence files with multiprocessing
        pool = mp.Pool(processes=num_threads)
        pool.map(rename_seq_worker, argument_lol_for_rename_seq_worker)
        pool.close()
        pool.join()

    print('Done!')


if __name__ == '__main__':

    rename_seq_parser = argparse.ArgumentParser()
    rename_seq_parser.add_argument('-in',         required=True,                          help='input sequence file')
    rename_seq_parser.add_argument('-x',          required=False, default='fasta',        help='file extension, default: fasta')
    rename_seq_parser.add_argument('-sep_in',     required=False, default=None,           help='separator for input sequences')
    rename_seq_parser.add_argument('-sep_out',    required=False, default=None,           help='separator for output sequences, default: same as sep_in')
    rename_seq_parser.add_argument('-n',          required=False, default=None, type=int, help='the number of columns to keep')
    rename_seq_parser.add_argument('-prefix',     required=False, default=None,           help='add prefix to sequence')
    rename_seq_parser.add_argument('-oneline',    required=False, action="store_true",    help='put sequence in single line')
    rename_seq_parser.add_argument('-t',          required=False, type=int, default=1,    help='number of threads')

    args = vars(rename_seq_parser.parse_args())
    rename_seq(args)
