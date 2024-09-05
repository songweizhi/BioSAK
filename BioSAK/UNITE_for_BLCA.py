import argparse
from Bio import SeqIO
from BioSAK.global_functions import sep_path_basename_ext


UNITE_for_BLCA_usage = '''
========================= UNITE_for_BLCA example commands =========================

BioSAK UNITE_for_BLCA -i sh_general_release_dynamic_04.04.2024.fasta

# Input file is available from: https://doi.plutof.ut.ee/doi/10.15156/BIO/2959332

# output_file:
sh_general_release_dynamic_04.04.2024_blca_compatible.fasta
sh_general_release_dynamic_04.04.2024_blca_compatible.taxonomy

==================================================================================
'''


def UNITE_for_BLCA(args):

    db_file = args['i']

    db_file_path, db_file_basename, _ = sep_path_basename_ext(db_file)
    file_out_sequence = '%s/%s_BLCA_compatible.fasta'    % (db_file_path, db_file_basename)
    file_out_taxonomy = '%s/%s_BLCA_compatible.taxonomy' % (db_file_path, db_file_basename)

    file_out_sequence_handle = open(file_out_sequence, 'w')
    file_out_taxonomy_handle = open(file_out_taxonomy, 'w')
    for seq_record in SeqIO.parse(db_file, 'fasta'):
        seq_id_split = seq_record.id.split('|k__Fungi;')
        seq_id_new = seq_id_split[0]
        seq_id_new = seq_id_new.replace('|reps', '')
        seq_id_new = seq_id_new.replace('|refs', '')
        seq_id_new = '_'.join(seq_id_new.split('|')[1:])
        tax_str = seq_id_split[1]
        tax_str_split = tax_str.strip().split(';')
        tax_str_split_reverse = tax_str_split[::-1]
        tax_str_split_reverse.append('k__Fungi;')
        tax_str_split_reverse_str = ';'.join(tax_str_split_reverse)
        tax_str_split_reverse_str = tax_str_split_reverse_str.replace('s__', 'species:')
        tax_str_split_reverse_str = tax_str_split_reverse_str.replace('g__', 'genus:')
        tax_str_split_reverse_str = tax_str_split_reverse_str.replace('f__', 'family:')
        tax_str_split_reverse_str = tax_str_split_reverse_str.replace('o__', 'order:')
        tax_str_split_reverse_str = tax_str_split_reverse_str.replace('c__', 'class:')
        tax_str_split_reverse_str = tax_str_split_reverse_str.replace('p__', 'phylum:')
        tax_str_split_reverse_str = tax_str_split_reverse_str.replace('k__', 'superkingdom:')

        # write out to taxonomy file
        file_out_taxonomy_handle.write('%s\t%s\n' % (seq_id_new, tax_str_split_reverse_str))

        # write out to sequence file
        file_out_sequence_handle.write('>%s\n' % seq_id_new)
        file_out_sequence_handle.write('%s\n' % seq_record.seq)

    file_out_sequence_handle.close()
    file_out_taxonomy_handle.close()

    print('You may want to run makeblastdb with:')
    print('makeblastdb -in %s -dbtype nucl -parse_seqids' % file_out_sequence)


if __name__ == "__main__":

    UNITE_for_BLCA_parser = argparse.ArgumentParser()
    UNITE_for_BLCA_parser.add_argument('-i', required=True, help='sh_general_release_dynamic_04.04.2024.fasta')
    args = vars(UNITE_for_BLCA_parser.parse_args())
    UNITE_for_BLCA(args)
