#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from BioSAK.global_functions import sep_path_basename_ext


GTDB_for_BLCA_usage = '''
====================== GTDB_for_BLCA example commands ======================

BioSAK GTDB_for_BLCA -GTDB_ssu bac120_ar122_ssu_r89.fna

# output_file:
GTDB_bac120_ar122_ssu_r89_BLCAparsed.fasta
GTDB_bac120_ar122_ssu_r89_BLCAparsed.taxonomy

============================================================================
'''


def GTDB_for_BLCA(args):

    GTDB_db_file = args['GTDB_ssu']

    GTDB_db_file_path, GTDB_db_file_basename, GTDB_db_file_ext = sep_path_basename_ext(GTDB_db_file)
    file_out_sequence = '%s/%s_BLCAparsed.fasta'    % (GTDB_db_file_path, GTDB_db_file_basename)
    file_out_taxonomy = '%s/%s_BLCAparsed.taxonomy' % (GTDB_db_file_path, GTDB_db_file_basename)

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
        file_out_sequence_handle.write('>%s\n' % seq_record.id)
        file_out_sequence_handle.write('%s\n' % seq_record.seq)

    file_out_sequence_handle.close()
    file_out_taxonomy_handle.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-GTDB_ssu', required=True, help='GTDB SSU file, e.g. bac120_ar122_ssu_r89.fna')

    args = vars(parser.parse_args())
    GTDB_for_BLCA(args)
