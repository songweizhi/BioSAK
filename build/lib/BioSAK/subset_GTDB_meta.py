import os
import argparse


subset_GTDB_meta_usage = '''
================================ subset_GTDB_meta example commands ================================

BioSAK subset_GTDB_meta -meta ar122_metadata_r202.tsv -id gnm_id.txt -out metadata_subset.txt
BioSAK subset_GTDB_meta -meta bac120_metadata_r202.tsv -id gnm_id.txt -out metadata_subset.txt

# format of gnm_id.txt (one id per line)
GCA_002726655.1
GCA_002728675.1

===================================================================================================
'''


def subset_GTDB_meta(args):

    gnm_metadata   = args['meta']
    interested_gnm = args['id']
    output_txt     = args['out']

    interested_gnm_set = set()
    for each_gnm in open(interested_gnm):
        interested_gnm_set.add(each_gnm.strip())

    found_gnm_set = set()
    output_txt_handle = open(output_txt, 'w')
    for each_ref in open(gnm_metadata):
        each_ref_split = each_ref.strip().split('\t')
        if each_ref.startswith('accession	ambiguous_bases'):
            output_txt_handle.write(each_ref)
        else:
            ref_accession = each_ref_split[0][3:]
            if ref_accession in interested_gnm_set:
                output_txt_handle.write(each_ref)
                found_gnm_set.add(ref_accession)
    output_txt_handle.close()

    if len(found_gnm_set) < len(interested_gnm_set):
        not_found_gnm_list = []
        for each_gnm in interested_gnm_set:
            if each_gnm not in found_gnm_set:
                not_found_gnm_list.append(each_gnm)
        print('The following genomes were not found:')
        print(','.join(not_found_gnm_list))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-meta', required=True, help='GTDB reference genome metadata')
    parser.add_argument('-id',   required=True, help='id of genomes')
    parser.add_argument('-out',  required=True, help='output file')
    args = vars(parser.parse_args())
    subset_GTDB_meta(args)
