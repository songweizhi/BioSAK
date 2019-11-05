import os
import glob
import shutil
from Bio import SeqIO


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


# define input and output files
faa_file_folder = '/Users/songweizhi/Desktop/test'
output_folder = '/Users/songweizhi/Desktop/test_aa_comp'


force_create_folder(output_folder)

amino_acid_str = 'ARNDBCEQZGHILKMFPSTWYV'
amino_acid_list = [i for i in amino_acid_str]
faa_file_re = '%s/*.faa' % faa_file_folder
faa_file_list = [os.path.basename(file_name) for file_name in glob.glob(faa_file_re)]

for faa_file in faa_file_list:

    pwd_faa = '%s/%s' % (faa_file_folder, faa_file)
    faa_file_path, faa_file_basename, faa_file_extension = sep_path_basename_ext(pwd_faa)

    output_file = '%s/%s_aa_comp.txt' % (output_folder, faa_file_basename)
    output_file_handle = open(output_file, 'w')
    output_file_handle.write('Gene\t%s\n' % '\t'.join(amino_acid_list))
    for seq_record in SeqIO.parse(pwd_faa, 'fasta'):
        seq_record_id = seq_record.id
        seq_record_seq = str(seq_record.seq)

        aa_count_list = []
        for aa in amino_acid_str:
            aa_percent = seq_record_seq.count(aa)*100/len(seq_record_seq)
            aa_percent = float("{0:.2f}".format(aa_percent))
            aa_count_list.append(aa_percent)

        aa_count_list_str = [str(i) for i in aa_count_list]
        output_file_handle.write('%s\t%s\n' % (seq_record_id, '\t'.join(aa_count_list_str)))

    output_file_handle.close()

