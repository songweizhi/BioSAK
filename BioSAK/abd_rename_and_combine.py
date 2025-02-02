import os
import glob
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def rename_and_cat_gnm(gnm_dir, gnm_ext, combined_renamed_fna, rename_txt):

    gnm_file_re = '%s/*.%s' % (gnm_dir, gnm_ext)
    gnm_file_list = glob.glob(gnm_file_re)

    gnm_rename_txt_handle = open(rename_txt, 'w')
    combined_renamed_fa_handle = open(combined_renamed_fna, 'w')
    for each_gnm in sorted(gnm_file_list):
        _, _, gnm_base, _ = sep_path_basename_ext(each_gnm)
        gnm_id_no_underscore = gnm_base.replace('_', '')
        gnm_rename_txt_handle.write('%s\t%s\n' % (gnm_base, gnm_id_no_underscore))

        seq_index = 1
        for each_seq in SeqIO.parse(each_gnm, 'fasta'):
            combined_renamed_fa_handle.write('>%s_%s\n%s\n' % (gnm_id_no_underscore, seq_index, str(each_seq.seq)))
            seq_index += 1

    gnm_rename_txt_handle.close()
    combined_renamed_fa_handle.close()


# file in
gnm_dir = '/Users/songweizhi/Desktop/3_combined_genomes_50_5_614_dRep99_406'
gnm_ext = 'fna'

# file out
combined_renamed_fa = '/Users/songweizhi/Desktop/Sponge_r220/11_host_specificity_by_abundance/dRep99_406_renamed_combined.fna'
gnm_rename_txt      = '/Users/songweizhi/Desktop/Sponge_r220/11_host_specificity_by_abundance/dRep99_406_renamed_combined.txt'

rename_and_cat_gnm(gnm_dir, gnm_ext, combined_renamed_fa, gnm_rename_txt)
