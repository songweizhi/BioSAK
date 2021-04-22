import os
import argparse
from datetime import datetime
from BioSAK.global_functions import time_format


top_16S_hits_usage = '''
======================================= top_16S_hits example commands =======================================

module load blast+/2.9.0
BioSAK top_16S_hits -p 16S_vs_GTDB -q 16S_sequences.fa -r GTDB_ssu_all_r95.fna -t 6
BioSAK top_16S_hits -p 16S_vs_SILVA -q 16S_sequences.fa -r SILVA_138_SSURef_NR99_tax_silva.fasta -t 6 -top 5

# prepare GTDB SSU database file
wget https://data.gtdb.ecogenomic.org/releases/release95/95.0/genomic_files_all/ssu_all_r95.tar.gz
tar xvzf ssu_all_r95.tar.gz
makeblastdb -in GTDB_ssu_all_r95.fna -dbtype nucl -parse_seqids -logfile /dev/null

# prepare SILVA SSU database file
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva.fasta.gz
gunzip SILVA_138_SSURef_NR99_tax_silva.fasta.gz
makeblastdb -in SILVA_138_SSURef_NR99_tax_silva.fasta -dbtype nucl -parse_seqids -logfile /dev/null

=============================================================================================================
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


def keep_top_n_blast_hit(file_in, file_out, top_hit_num):

    file_out_handle = open(file_out, 'w')

    current_query_id = ''
    current_query_wrote_hits = 0
    for blast_hit in open(file_in):

        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]

        if current_query_id == '':
            current_query_id = query_id
            file_out_handle.write(blast_hit)
            current_query_wrote_hits += 1

        elif current_query_id != '':
            if query_id == current_query_id:
                if current_query_wrote_hits < top_hit_num:
                    file_out_handle.write(blast_hit)
                    current_query_wrote_hits += 1
            else:
                current_query_id = query_id
                file_out_handle.write(blast_hit)
                current_query_wrote_hits = 1

    file_out_handle.close()


def top_16S_hits(args):

    output_prefix = args['p']
    query_seq     = args['q']
    gtdb_ssu_file = args['r']
    num_threads   = args['t']
    evalue_cutoff = args['evalue']
    top_hit_num   = args['top']

    # define output file name
    blast_op            = '%s_blastn.tab'           % output_prefix
    blast_op_best_hit   = '%s_blastn_best_hit.tab'  % output_prefix
    op_classified       = '%s_classifications.txt'  % output_prefix
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

            each_line_split = each_line.strip()[1:].split(' ')
            ref_seq_id = each_line_split[0]

            # get ref_seq_des
            ref_seq_des = ' '.join(each_line_split[1:])
            if ' [' in ref_seq_des:
                ref_seq_des = ref_seq_des.split(' [')[0]
            if ' ' in ref_seq_des:
                ref_seq_des = '_'.join(ref_seq_des.split(' '))

            ref_seq_to_taxon_dict[ref_seq_id] = ref_seq_des


    # run blast and keep best hits
    print(datetime.now().strftime(time_format) + 'Running blast with %s cores, be patient!' % num_threads)
    blastn_cmd = 'blastn -query %s -db %s -out %s -outfmt 6 -evalue %s -num_threads %s' % (query_seq, gtdb_ssu_file, blast_op, evalue_cutoff, num_threads)
    os.system(blastn_cmd)

    print(datetime.now().strftime(time_format) + 'Keep the top %s hits for each query' % top_hit_num)
    if top_hit_num == 1:
        keep_best_blast_hit(blast_op, blast_op_best_hit)
    else:
        keep_top_n_blast_hit(blast_op, blast_op_best_hit, top_hit_num)

    # Parse blast results
    print(datetime.now().strftime(time_format) + 'Parse blast results')
    op_classified_handle = open(op_classified, 'w')
    op_classified_handle.write('Query\tReference\tIdentity\tAligned_length\tTaxonomy\n')
    classified_queries = set()
    for each_hit in open(blast_op_best_hit):
        each_hit_split  = each_hit.strip().split('\t')
        query_id        = each_hit_split[0]
        ref_id          = each_hit_split[1]
        identity        = float(each_hit_split[2])
        aln_len         = int(each_hit_split[3])
        ref_taxon       = ref_seq_to_taxon_dict.get(ref_id, 'NA')
        if ref_taxon != 'NA':
            op_classified_handle.write('%s\t%s\t%s\t%s\t%s\n' % (query_id, ref_id, identity, aln_len, ref_taxon))
            classified_queries.add(query_id)
    op_classified_handle.close()

    # write out id of unclassified queries
    unclassified_queries = query_seq_set - classified_queries
    if len(unclassified_queries) > 0:
        op_unclassified_handle = open(op_unclassified, 'w')
        op_unclassified_handle.write('\n'.join(unclassified_queries))
        op_unclassified_handle.close()

    # final report
    print(datetime.now().strftime(time_format) + 'Classifications exported to %s'        % op_classified)
    if len(unclassified_queries) > 0:
        print(datetime.now().strftime(time_format) + 'Unclassified sequences exported to %s'  % op_unclassified)
    print(datetime.now().strftime(time_format) + 'Done!')


if __name__ == '__main__':

    top_16S_hits_parser = argparse.ArgumentParser()
    top_16S_hits_parser.add_argument('-p',           required=True,                           help='output prefix')
    top_16S_hits_parser.add_argument('-q',           required=True,                           help='query sequence file')
    top_16S_hits_parser.add_argument('-r',           required=True,                           help='SILVA or GTDB SSU sequence file')
    top_16S_hits_parser.add_argument('-evalue',      required=False, default='1e-20',         help='evalue cutoff, default: 1e-20')
    top_16S_hits_parser.add_argument('-top',         required=False, type=int, default=1,     help='Number of top hits to report, default: 1')
    top_16S_hits_parser.add_argument('-t',           required=False, type=int, default=1,     help='number of threads')
    args = vars(top_16S_hits_parser.parse_args())
    top_16S_hits(args)
