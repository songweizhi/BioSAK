import os
import glob
import argparse
from Bio import SeqIO
import multiprocessing as mp
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def gbk_to_ffn_faa(gbk_in, ffn_out, faa_out):

    _, _, gbk_base, _ = sep_path_basename_ext(gbk_in)
    faa_out_handle = open(faa_out, 'w')
    ffn_out_handle = open(ffn_out, 'w')
    for seq_record in SeqIO.parse(gbk_in, 'genbank'):
        record_sequence = str(seq_record.seq)
        for feature in seq_record.features:
            if feature.type == 'CDS':
                feature_pos = '%s-%s' % (feature.location.start, feature.location.end)
                feature_seq_nt = record_sequence[feature.location.start:feature.location.end]
                if feature.location.strand in ['-1', -1]:
                    feature_pos = '%s-%s' % (feature.location.end, feature.location.start)
                    feature_seq_nt = str(Seq(feature_seq_nt).reverse_complement())
                feature_seq_aa = ''
                if 'translation' in feature.qualifiers:
                    feature_seq_aa = feature.qualifiers['translation'][0]

                protein_id = ''
                if 'protein_id' in feature.qualifiers:
                    protein_id = feature.qualifiers['protein_id'][0]

                cds_id_to_used = feature_pos.replace('<', '').replace('>', '')
                if protein_id != '':
                    cds_id_to_used = protein_id

                feature_product = ''
                if 'protein_id' in feature.qualifiers:
                    feature_product = feature.qualifiers['product'][0]

                gene_name = ''
                if 'gene' in feature.qualifiers:
                    gene_name = feature.qualifiers['gene'][0]

                feature_note = ''
                if 'note' in feature.qualifiers:
                    feature_note = feature.qualifiers['note'][0]

                desc_to_use = ''
                if feature_product != '':
                    desc_to_use = feature_product
                elif feature_note != '':
                    desc_to_use = feature_note
                elif gene_name != '':
                    desc_to_use = gene_name

                # write out nt sequence
                ffn_out_handle.write('>%s__%s %s\n' % (gbk_base, cds_id_to_used, desc_to_use))
                ffn_out_handle.write('%s\n' % feature_seq_nt)

                # write out aa sequence
                if feature_seq_aa != '':
                    faa_out_handle.write('>%s__%s %s\n' % (gbk_base, cds_id_to_used, desc_to_use))
                    faa_out_handle.write('%s\n' % feature_seq_aa)
    faa_out_handle.close()
    ffn_out_handle.close()


file_dir = '/Users/songweizhi/Desktop/Sponge_r226/sponge_phylogeny/COI/RefSeq_COI_wd/tmp'
file_ext = 'gbk'
file_re = '%s/*.%s' % (file_dir, file_ext)
file_list = glob.glob(file_re)

for file in sorted(file_list):
    faa_file = file.replace('.gbk', '.faa')
    ffn_file = file.replace('.gbk', '.ffn')
    gbk_to_ffn_faa(file, ffn_file, faa_file)


# gbk_file = '/Users/songweizhi/Desktop/Sponge_r226/sponge_phylogeny/to_add/COI_wd/tmp/KR911863.1.gbk'
# fna_file = '/Users/songweizhi/Desktop/KR911863.1.fas'
# faa_file = '/Users/songweizhi/Desktop/KR911863.1.faa'
#
# gbk_to_ffn_faa(gbk_file, fna_file, faa_file)

