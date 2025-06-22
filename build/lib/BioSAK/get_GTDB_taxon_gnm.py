import os
import argparse


get_GTDB_taxon_gnm_parser_usage = '''
====================================== get_GTDB_taxon_gnm example commands ======================================

# Example commands
BioSAK get_GTDB_taxon_gnm -p Bac -meta bac120_metadata_r202.tsv -path genome_paths.tsv -taxon taxons.txt -cpl 85 -ctm 5
BioSAK get_GTDB_taxon_gnm -p Arc -meta ar122_metadata_r202.tsv -path genome_paths.tsv -taxon taxons.txt -cpl 50 -ctm 10

# -meta: bac120_metadata_r202.tsv or ar122_metadata_r202.tsv
# -path: genome_paths.tsv, in fastani folder from auxillary_files
# -taxon: one taxon per line, example below:
p__Thermoproteota
c__Methanosarcinia
f__Thalassarchaeaceae

=================================================================================================================
'''


def get_GTDB_taxon_gnm(args):

    output_prefix =             args['p']
    gtdb_gnm_metadata =         args['meta']
    gtdb_ref_gnm_path_txt =     args['path']
    taxons_to_retrieve_txt =    args['taxon']
    cpl_cutoff =                args['cpl']
    ctm_cutoff =                args['ctm']

    try:
        cpl_cutoff = int(cpl_cutoff)
    except:
        cpl_cutoff = cpl_cutoff

    try:
        ctm_cutoff = int(ctm_cutoff)
    except:
        ctm_cutoff = ctm_cutoff

    op_name_cpl_part = ''
    if cpl_cutoff is not None:
        op_name_cpl_part = '_%s' % cpl_cutoff

    op_name_ctm_part = ''
    if ctm_cutoff is not None:
        op_name_ctm_part = '_%s' % ctm_cutoff

    # define output file name
    gnms_to_retrieve_gtdb_txt = '%s%s%s_in_GTDB.txt'    % (output_prefix, op_name_cpl_part, op_name_ctm_part)
    gnms_to_retrieve_ncbi_txt = '%s%s%s_in_NCBI.txt'    % (output_prefix, op_name_cpl_part, op_name_ctm_part)
    gnms_taxonomy_txt         = '%s%s%s_taxonomy.txt'   % (output_prefix, op_name_cpl_part, op_name_ctm_part)

    already_exist_files = []
    for each_op in [gnms_to_retrieve_gtdb_txt, gnms_to_retrieve_ncbi_txt, gnms_taxonomy_txt]:
        if os.path.isfile(each_op) is True:
            already_exist_files.append(each_op)

    if len(already_exist_files) > 0:
        print('The following files already exist, use a different prefix with -p')
        print('\n'.join(already_exist_files))
        print('Program exited!')
        exit()

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
    gnm_cpl_dict = {}
    gnm_ctm_dict = {}
    gnm_size_dict = {}
    gnm_set_to_retrieve = set()
    gnm_to_taxon_dict = {}
    col_index = {}
    for each_ref in open(gtdb_gnm_metadata):
        each_ref_split = each_ref.strip().split('\t')
        if each_ref.startswith('accession	ambiguous_bases'):
            col_index = {key: i for i, key in enumerate(each_ref_split)}
        else:
            ref_accession = each_ref_split[0][3:]
            gnm_completeness  = float(each_ref_split[2])
            gnm_contamination = float(each_ref_split[3])
            gnm_size = each_ref_split[col_index['genome_size']]
            gtdb_taxonomy = each_ref_split[col_index['gtdb_taxonomy']]
            for each_interested_taxon in taxons_to_retrieve_set:
                with_semicolon = '%s;' % each_interested_taxon
                if with_semicolon in gtdb_taxonomy:

                    # check completeness
                    cpl_passed = False
                    if cpl_cutoff is None:
                        cpl_passed = True
                    else:
                        if gnm_completeness >= float(cpl_cutoff):
                            cpl_passed = True

                    # check contamination
                    ctm_passed = False
                    if ctm_cutoff is None:
                        ctm_passed = True
                    else:
                        if gnm_contamination <= float(ctm_cutoff):
                            ctm_passed = True

                    if (cpl_passed is True) and (ctm_passed is True):
                        gnm_set_to_retrieve.add(ref_accession)
                        gnm_to_taxon_dict[ref_accession] = gtdb_taxonomy
                        gnm_cpl_dict[ref_accession] = gnm_completeness
                        gnm_ctm_dict[ref_accession] = gnm_contamination
                        gnm_size_dict[ref_accession] = gnm_size

    # write out
    from_gtdb = 0
    from_ncbi = 0
    gnms_taxonomy_txt_handle = open(gnms_taxonomy_txt, 'w')
    gnms_taxonomy_txt_handle.write('Genome\tCompleteness\tContamination\tSize\tTaxon\n')
    gnms_to_retrieve_gtdb_txt_handle = open(gnms_to_retrieve_gtdb_txt, 'w')
    gnms_to_retrieve_ncbi_txt_handle = open(gnms_to_retrieve_ncbi_txt, 'w')
    for each_gnm_to_retrieve in sorted([i for i in gnm_set_to_retrieve]):
        gnm_taxon = gnm_to_taxon_dict.get(each_gnm_to_retrieve, 'NA')
        gnms_taxonomy_txt_handle.write('%s\t%s\t%s\t%s\t%s\n' % (each_gnm_to_retrieve, gnm_cpl_dict[each_gnm_to_retrieve], gnm_ctm_dict[each_gnm_to_retrieve], gnm_size_dict[each_gnm_to_retrieve], gnm_taxon))
        if each_gnm_to_retrieve in gtdb_ref_gnm_set:
            gnms_to_retrieve_gtdb_txt_handle.write('%s\n' % each_gnm_to_retrieve)
            from_gtdb += 1
        else:
            gnms_to_retrieve_ncbi_txt_handle.write('%s\n' % each_gnm_to_retrieve)
            from_ncbi += 1
    gnms_taxonomy_txt_handle.close()
    gnms_to_retrieve_gtdb_txt_handle.close()
    gnms_to_retrieve_ncbi_txt_handle.close()

    gnms_to_retrieve_gtdb_txt_with_num = '%s%s%s_in_GTDB_%s.txt'        % (output_prefix, op_name_cpl_part, op_name_ctm_part, from_gtdb)
    gnms_to_retrieve_ncbi_txt_with_num = '%s%s%s_in_NCBI_%s.txt'        % (output_prefix, op_name_cpl_part, op_name_ctm_part, from_ncbi)

    os.system('mv %s %s' % (gnms_to_retrieve_gtdb_txt, gnms_to_retrieve_gtdb_txt_with_num))
    os.system('mv %s %s' % (gnms_to_retrieve_ncbi_txt, gnms_to_retrieve_ncbi_txt_with_num))

    print('Total number of genomes from interested taxons: %s'  % len(gnm_set_to_retrieve))
    print('Genomes in GTDB: %s, see details in %s'              % (from_gtdb, gnms_to_retrieve_gtdb_txt_with_num))
    print('Genomes in NCBI: %s, see details in %s'              % (from_ncbi, gnms_to_retrieve_ncbi_txt_with_num))
    print('Genome taxon information is in %s'                   % (gnms_taxonomy_txt))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p',           required=False, default='Genomes',          help='output prefix')
    parser.add_argument('-meta',        required=True,                              help='GTDB reference genome metadata')
    parser.add_argument('-path',        required=True,                              help='GTDB genome_paths.tsv')
    parser.add_argument('-taxon',       required=True,                              help='interested taxons, one taxon per line')
    parser.add_argument('-cpl',         required=False, default=None, type=float,   help='completeness cutoff (0-100), default: None')
    parser.add_argument('-ctm',         required=False, default=None, type=float,   help='contamination cutoff, default: None')
    args = vars(parser.parse_args())
    get_GTDB_taxon_gnm(args)
