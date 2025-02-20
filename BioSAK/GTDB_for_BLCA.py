#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from BioSAK.global_functions import sep_path_basename_ext


GTDB_for_BLCA_usage = '''
===== GTDB_for_BLCA example commands =====

BioSAK GTDB_for_BLCA -i ssu_all_r220.fna

# output_file:
ssu_all_r220.blca.fa
ssu_all_r220.blca.tax

==========================================
'''


def GTDB_for_BLCA(args):

    GTDB_db_file = args['i']

    db_path, db_base, db_ext = sep_path_basename_ext(GTDB_db_file)
    file_out_sequence = '%s/%s.blca.fa'     % (db_path, db_base)
    file_out_taxonomy = '%s/%s.blca.tax'    % (db_path, db_base)

    rank_list = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']

    file_out_sequence_handle = open(file_out_sequence, 'w')
    file_out_taxonomy_handle = open(file_out_taxonomy, 'w')
    for seq_record in SeqIO.parse(GTDB_db_file, 'fasta'):
        seq_record_taxon_split = ' '.join(seq_record.description.strip().split(' [')[0].split(' ')[1:]).split(';')
        seq_record_taxon_split_no_rank = [i[3:] for i in seq_record_taxon_split]
        seq_record_taxon_split_no_rank_reverse = seq_record_taxon_split_no_rank[::-1]

        GTDB_seq_taxon_str = ''
        n = 0
        for taxon_rank in rank_list:
            GTDB_seq_taxon_str += '%s:%s;' % (taxon_rank, seq_record_taxon_split_no_rank_reverse[n])
            n += 1

        # write out to taxonomy file
        file_out_taxonomy_handle.write('%s\t%s\n' % (seq_record.id, GTDB_seq_taxon_str))

        # write out to sequence file
        file_out_sequence_handle.write('>%s\n' % seq_record.description)
        file_out_sequence_handle.write('%s\n' % seq_record.seq)

    file_out_sequence_handle.close()
    file_out_taxonomy_handle.close()

    print('Done!')
    print('You may need to run:')
    print('makeblastdb -dbtype nucl -parse_seqids -in %s' % file_out_sequence)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, help='GTDB SSU file, e.g. ssu_all_r220.fna')
    args = vars(parser.parse_args())
    GTDB_for_BLCA(args)
