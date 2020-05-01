import os
import glob
import argparse
from Bio import SeqIO


def get_ref_taxonomy(args):

    fsa_folder =    args['i']
    file_out =      args['o']
    taxonomy_file = args['t']

    fsa_file_re = '%s/*.hitdb.fsa' % fsa_folder
    fsa_file_list = [os.path.basename(file_name) for file_name in glob.glob(fsa_file_re)]

    # store taxonomy info in dict
    taxonomy_dict = {}
    for ref in open(taxonomy_file):
        ref_split = ref.strip().split('\t')
        ref_id = ref_split[0]
        ref_taxon = ref_split[1]
        ref_taxon_split_no_rank = [i.split(':')[1] for i in ref_taxon.split(';')[:-1]]
        ref_taxon_no_rank = ';'.join(ref_taxon_split_no_rank[::-1])
        taxonomy_dict[ref_id] = ref_taxon_no_rank

    file_out_handle = open(file_out, 'w')
    for fsa_file in fsa_file_list:
        pwd_fas_file = '%s/%s' % (fsa_folder, fsa_file)
        for seq in SeqIO.parse(pwd_fas_file, 'fasta'):
            if seq.id in taxonomy_dict:
                file_out_handle.write('%s\t%s\t%s\n' % (fsa_file[:-10], seq.id, taxonomy_dict[seq.id]))


if __name__ == '__main__':

    get_best_hit_parser = argparse.ArgumentParser()

    # arguments for COG_parser
    get_best_hit_parser.add_argument('-i',   required=True,  help='folder holds .hitdb.fsa files')
    get_best_hit_parser.add_argument('-t',   required=True,  help='taxonomy file')
    get_best_hit_parser.add_argument('-o',   required=True,  help='output file')

    args = vars(get_best_hit_parser.parse_args())

    get_ref_taxonomy(args)




fsa_folder = '/Users/songweizhi/Desktop/OTU_table_GTDB_BLCA'
taxonomy_file = '/Users/songweizhi/Desktop/OTU_table_GTDB_BLCA/00_DataNeeded/BLCA/db_GTDB_SSU/GTDB_bac120_ar122_ssu_r89_BLCAparsed.taxonomy'
file_out = '/Users/songweizhi/Desktop/OTU_table_GTDB_BLCA.txt'
