import os
import glob
import shutil
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def remove_dash(aln_in, seq_out, min_len):
    seq_out_handle = open(seq_out, 'w')
    for each_seq in SeqIO.parse(aln_in, 'fasta'):
        each_seq_sequence = str(each_seq.seq)
        each_seq_sequence_no_dash = ''
        for each_base in each_seq_sequence:
            if each_base != '-':
                each_seq_sequence_no_dash += each_base
        if len(each_seq_sequence_no_dash) >= min_len:
            export_dna_record(each_seq_sequence_no_dash, each_seq.id, '', seq_out_handle)
    seq_out_handle.close()


parser = argparse.ArgumentParser()

parser.add_argument('-in',
                    required=True,
                    help='path to folder bac120_msa_individual_genes_r83')

parser.add_argument('-makeblastdb',
                    required=False,
                    default='makeblastdb',
                    help='path to makeblastdb')

args = vars(parser.parse_args())
folder_in = args['in']
pwd_makeblastdb = args['makeblastdb']


folder_path, folder_name = os.path.split(folder_in)
folder_out = '%s_db' % folder_name
pwd_folder_out = '%s/%s' % (folder_path, folder_out)

extension = 'faa'
min_len = 20

# create output folder
if os.path.isdir(folder_out):
    shutil.rmtree(folder_out)
    os.makedirs(folder_out)
else:
    os.makedirs(folder_out)

# get file list
file_re = '%s/*.%s' % (folder_in, extension)
file_list = [os.path.basename(file_name) for file_name in glob.glob(file_re)]

# remove dash for all files
n = 1
for each_file in file_list:
    print('Processing %s/%s: %s' % (n, len(file_list), each_file))

    pwd_each_file = '%s/%s' % (folder_in, each_file)
    pwd_each_file_out = '%s/%s' % (folder_out, each_file)
    remove_dash(pwd_each_file, pwd_each_file_out, min_len)
    os.system('%s -in %s -dbtype prot -parse_seqids >/dev/null ' % (pwd_makeblastdb, pwd_each_file_out))
    n += 1

print('makeblastdb done!')
