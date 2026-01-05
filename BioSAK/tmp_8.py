import os
import glob
import argparse
from Bio import SeqIO
import multiprocessing as mp


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def gbk2faa(gbk_in, faa_out):

    _, _, faa_base, _ = sep_path_basename_ext(faa_out)

    faa_out_handle = open(faa_out, 'w')
    for seq_record in SeqIO.parse(gbk_in, 'genbank'):
        for feature in seq_record.features:
            if feature.type == 'CDS':
                feature_translation = ''
                if 'translation' in feature.qualifiers:
                    feature_translation = feature.qualifiers['translation'][0]
                if feature_translation != '':
                    faa_out_handle.write('>%s__%s\n' % (faa_base, feature.qualifiers['protein_id'][0]))
                    faa_out_handle.write('%s\n' % feature_translation)
    faa_out_handle.close()



file_dir = '/Users/songweizhi/Desktop/Sponge_r226/sponge_phylogeny/COI/4110_wd/tmp'
file_ext = 'gbk'
file_re = '%s/*.%s' % (file_dir, file_ext)
file_list = glob.glob(file_re)


for file in sorted(file_list):
    print(file)
    faa_file = file.replace('.gbk', '.faa')
    gbk2faa(file, faa_file)









