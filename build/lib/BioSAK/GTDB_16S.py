import os
import argparse
from datetime import datetime
from BioSAK.global_functions import time_format


gtdb_16s_usage = '''
===================================== GTDB_16S example commands =====================================

module load blast+/2.9.0
BioSAK GTDB_16S -p 16S_sequences -q 16S_sequences.fa -r GTDB_ssu_all_r95.fna -t 12

# prepare database file
wget https://data.gtdb.ecogenomic.org/releases/release95/95.0/genomic_files_all/ssu_all_r95.tar.gz
tar xvzf ssu_all_r95.tar.gz
makeblastdb -in GTDB_ssu_all_r95.fna -dbtype nucl -parse_seqids -logfile /dev/null

=====================================================================================================
'''


def keep_best_blast_hit(file_in, file_out):

    file_out_handle = open(file_out, 'w')
    best_hit_line = ''
    best_hit_query_id = ''
    best_hit_score = 0
    for blast_hit in open(file_in):
        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]
        bit_score = float(blast_hit_split[11])

        if best_hit_query_id == '':
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

        elif (query_id == best_hit_query_id) and (bit_score > best_hit_score):
            best_hit_score = bit_score
            best_hit_line = blast_hit

        elif query_id != best_hit_query_id:
            file_out_handle.write(best_hit_line)
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

    file_out_handle.write(best_hit_line)
    file_out_handle.close()


def gtdb_16s(args):

    output_prefix = args['p']
    query_seq     = args['q']
    gtdb_ssu_file = args['r']
    num_threads   = args['t']
    evalue_cutoff = args['evalue']

    # define output file name
    blast_op            = '%s_blastn.tab'           % output_prefix
    blast_op_best_hit   = '%s_blastn_best_hit.tab'  % output_prefix
    op_classified       = '%s_classified.txt'       % output_prefix
    op_unclassified     = '%s_unclassified.txt'     % output_prefix

    # get_query_seq_list
    print(datetime.now().strftime(time_format) + 'Get the id of query sequences')
    query_seq_set = set()
    for each_query in open(query_seq):
        if each_query.startswith('>'):
            query_id = each_query.strip()[1:].split(' ')[0]
            query_seq_set.add(query_id)

    # get ref_seq_to_taxon_dict
    print(datetime.now().strftime(time_format) + 'Read in database file')
    ref_seq_to_taxon_dict = {}
    for each_line in open(gtdb_ssu_file):
        if each_line.startswith('>'):
            each_line_split = each_line[1:].split(' [')[0].split(' ')
            ref_seq_id = each_line_split[0]
            ref_seq_taxon = '_'.join(each_line_split[1:])
            ref_seq_to_taxon_dict[ref_seq_id] = ref_seq_taxon

    # run blast and keep best hits
    print(datetime.now().strftime(time_format) + 'Running blast with %s cores, be patient!' % num_threads)
    blastn_cmd = 'blastn -query %s -db %s -out %s -outfmt 6 -evalue %s -num_threads %s' % (query_seq, gtdb_ssu_file, blast_op, evalue_cutoff, num_threads)
    os.system(blastn_cmd)
    print(datetime.now().strftime(time_format) + 'Keep only best hits')
    keep_best_blast_hit(blast_op, blast_op_best_hit)

    # get query_to_ref_dict
    print(datetime.now().strftime(time_format) + 'Read in blast results')
    query_to_ref_dict = {}
    for each_hit in open(blast_op_best_hit):
        each_hit_split = each_hit.strip().split('\t')
        query_id = each_hit_split[0]
        ref_id = each_hit_split[1]
        query_to_ref_dict[query_id] = ref_id

    # get taxonomy of query sequences
    print(datetime.now().strftime(time_format) + 'Get taxonomy of query sequences')
    op_classified_handle = open(op_classified, 'w')
    op_unclassified_handle = open(op_unclassified, 'w')
    for each_query_seq in sorted([i for i in query_seq_set]):
        matched_ref = query_to_ref_dict.get(each_query_seq, None)
        matched_ref_taxon = ref_seq_to_taxon_dict.get(matched_ref, None)
        if matched_ref_taxon is None:
            op_unclassified_handle.write('%s\n' % each_query_seq)
        else:
            op_classified_handle.write('%s\t%s\n' % (each_query_seq, matched_ref_taxon))
    op_classified_handle.close()
    op_unclassified_handle.close()

    # final report
    print(datetime.now().strftime(time_format) + 'Classification results exported to %s'        % op_classified)
    print(datetime.now().strftime(time_format) + 'ID of unclassified sequences exported to %s'  % op_unclassified)
    print(datetime.now().strftime(time_format) + 'Done!')


if __name__ == '__main__':

    gtdb_16s_parser = argparse.ArgumentParser()
    gtdb_16s_parser.add_argument('-p',           required=True,                           help='output prefix')
    gtdb_16s_parser.add_argument('-q',           required=True,                           help='query sequences')
    gtdb_16s_parser.add_argument('-r',           required=True,                           help='GTDB SSU sequences file')
    gtdb_16s_parser.add_argument('-evalue',      required=False, default='1e-20',         help='evalue cutoff, default: 1e-20')
    gtdb_16s_parser.add_argument('-t',           required=False, type=int, default=1,     help='number of threads')
    args = vars(gtdb_16s_parser.parse_args())
    gtdb_16s(args)
