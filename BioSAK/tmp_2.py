import os
import glob
from Bio import SeqIO


file_list   = glob.glob('b_rep_set_splitted/*.fa')
blastn_dir  = 'blastn_splitted'

seq_to_file_dict = dict()
for each_file in file_list:
    for each_seq in SeqIO.parse(each_file, 'fasta'):
        seq_id = each_seq.id
        seq_to_file_dict[seq_id] = os.path.basename(each_file)


dict = dict()
for each_line in open('b_rep_set.fasta.blastn'):

    seq_id = each_line.split('\t')[0]
    fa_file_name = seq_to_file_dict[seq_id]
    blastn_file = '%s/%s.blastn' % (blastn_dir, fa_file_name)

    if blastn_file not in dict:
        dict[blastn_file] = []

    dict[blastn_file].append(each_line.strip())


for each_b_file in dict:
    lines = dict[each_b_file]
    each_b_file_handle = open(each_b_file, 'w')
    each_b_file_handle.write('\n'.join(lines))
    each_b_file_handle.close()

