import argparse


get_GTDB_taxon_gnm_parser_usage = '''
================================== get_GTDB_taxon_gnm example commands ==================================

# Example commands
BioSAK get_GTDB_taxon_gnm -meta bac120_metadata_r202.tsv -path genome_paths.tsv -taxon taxons.txt
BioSAK get_GTDB_taxon_gnm -meta ar122_metadata_r202.tsv -path genome_paths.tsv -taxon taxons.txt

# genome metadata: bac120_metadata_r202.tsv or ar122_metadata_r202.tsv
# genome paths: in the fastani folder from auxillary_files
# example of interested taxon file (one taxon per line):
p__Thermoproteota
c__Methanosarcinia
f__Thalassarchaeaceae

=====================================================================================================
'''


def get_GTDB_taxon_gnm(args):

    gtdb_gnm_metadata =         args['meta']
    gtdb_ref_gnm_path_txt =     args['path']
    taxons_to_retrieve_txt =    args['taxon']
    output_prefix =             args['p']

    gnms_to_retrieve_gtdb_txt = '%s_in_GTDB.txt' % output_prefix
    gnms_to_retrieve_ncbi_txt = '%s_in_NCBI.txt' % output_prefix

    # read in taxons_to_retrieve_txt
    taxons_to_retrieve_set = set()
    for each_taxon in open(taxons_to_retrieve_txt):
        taxons_to_retrieve_set.add(each_taxon.strip())

    # read in gtdb_ref_gnm_path_txt
    gtdb_ref_gnm_set = set()
    for each_gnm_path in open(gtdb_ref_gnm_path_txt):
        gnm_path_split = each_gnm_path.strip().split(' ')
        gnm_accession = gnm_path_split[0].split('_genomic')[0]
        gtdb_ref_gnm_set.add(gnm_accession)

    # get genomes from interested taxons
    gnm_set_to_retrieve = set()
    col_index = {}
    for each_ref in open(gtdb_gnm_metadata):
        each_ref_split = each_ref.strip().split('\t')
        if each_ref.startswith('accession	ambiguous_bases'):
            col_index = {key: i for i, key in enumerate(each_ref_split)}
        else:
            ref_accession = each_ref_split[0][3:]
            gtdb_taxonomy = each_ref_split[col_index['gtdb_taxonomy']]
            for each_interested_taxon in taxons_to_retrieve_set:
                if each_interested_taxon in gtdb_taxonomy:
                    gnm_set_to_retrieve.add(ref_accession)

    # write out
    from_gtdb = 0
    from_ncbi = 0
    gnms_to_retrieve_gtdb_txt_handle = open(gnms_to_retrieve_gtdb_txt, 'w')
    gnms_to_retrieve_ncbi_txt_handle = open(gnms_to_retrieve_ncbi_txt, 'w')
    for each_gnm_to_retrieve in sorted([i for i in gnm_set_to_retrieve]):
        if each_gnm_to_retrieve in gtdb_ref_gnm_set:
            gnms_to_retrieve_gtdb_txt_handle.write('%s\n' % each_gnm_to_retrieve)
            from_gtdb += 1
        else:
            gnms_to_retrieve_ncbi_txt_handle.write('%s\n' % each_gnm_to_retrieve)
            from_ncbi += 1
    gnms_to_retrieve_gtdb_txt_handle.close()
    gnms_to_retrieve_ncbi_txt_handle.close()

    print('Total number of genomes from interested taxons: %s' % len(gnm_set_to_retrieve))
    print('Genomes in GTDB: %s, see details in %s'   % (from_gtdb, gnms_to_retrieve_gtdb_txt))
    print('Genomes in NCBI: %s, see details in %s' % (from_ncbi, gnms_to_retrieve_ncbi_txt))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-meta',        required=True,              help='GTDB reference genome metadata')
    parser.add_argument('-path',        required=True,              help='GTDB genome paths')
    parser.add_argument('-taxon',       required=True,              help='interested taxons, one taxon per line')
    parser.add_argument('-p',           required=False, default='', help='output prefix')
    args = vars(parser.parse_args())
    get_GTDB_taxon_gnm(args)
