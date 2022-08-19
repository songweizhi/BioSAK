import os
import argparse


sampling_GTDB_gnms_usage = '''
====================================== sampling_GTDB_gnms example commands ======================================

# Example commands
BioSAK sampling_GTDB_gnms -p Firm_o -meta bac120_metadata_r202.tsv -taxon p__Firmicutes -r o
BioSAK sampling_GTDB_gnms -p Firm_f -meta bac120_metadata_r202.tsv -taxon p__Firmicutes -r f -cpl 85 -ctm 5
BioSAK sampling_GTDB_gnms -p Firm_g -meta bac120_metadata_r202.tsv -taxon p__Firmicutes -r g -cpl 85 -ctm 5 -ts

# -meta: bac120_metadata_r202.tsv or ar122_metadata_r202.tsv
# Sampling rank need to be lower than the rank of taxon specified with "-taxon".

=================================================================================================================
'''


def sampling_GTDB_gnms(args):

    output_prefix     = args['p']
    gtdb_gnm_metadata = args['meta']
    interested_taxon  = args['taxon']
    sampling_rank     = args['r']
    gnm_num_per_taxon = args['n']
    cpl_cutoff        = args['cpl']
    ctm_cutoff        = args['ctm']
    type_strain_only  = args['ts']

    # define output file name
    selected_gnm_id_txt   = '%s_selected_gnm_id.txt'       % output_prefix
    selected_gnm_meta_txt = '%s_selected_gnm_metadata.txt' % output_prefix

    to_exit = False
    if os.path.isfile(selected_gnm_id_txt) is True:
        print('%s detected, please use a different prefix' % selected_gnm_id_txt)
        to_exit = True

    if os.path.isfile(selected_gnm_meta_txt) is True:
        print('%s detected, please use a different prefix' % selected_gnm_meta_txt)
        to_exit = True

    if to_exit is True:
        print('Program exited!')
        exit()

    dod = {}
    gnm_cpl_dict = {}
    gnm_ctm_dict = {}
    gnm_size_dict = {}
    gnm_metric_dict = {}
    gnm_to_taxon_dict = {}
    col_index = {}
    for each_ref in open(gtdb_gnm_metadata):
        each_ref_split = each_ref.strip().split('\t')
        if each_ref.startswith('accession	ambiguous_bases'):
            col_index = {key: i for i, key in enumerate(each_ref_split)}
        else:
            ref_accession = each_ref_split[0][3:]
            gnm_completeness = float(each_ref_split[2])
            gnm_contamination = float(each_ref_split[3])
            gnm_size = each_ref_split[col_index['genome_size']]
            gtdb_taxonomy = each_ref_split[col_index['gtdb_taxonomy']]
            is_type = each_ref_split[col_index['gtdb_type_designation']]

            if interested_taxon in gtdb_taxonomy:

                # check type
                type_pass = False
                if type_strain_only is False:
                    type_pass = True
                else:
                    if is_type in ['type strain of species', 'type strain of species']:
                        type_pass = True

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

                if (type_pass is True) and (cpl_passed is True) and (ctm_passed is True):

                    # add to dict
                    gnm_to_taxon_dict[ref_accession] = gtdb_taxonomy
                    gnm_cpl_dict[ref_accession] = gnm_completeness
                    gnm_ctm_dict[ref_accession] = gnm_contamination
                    gnm_size_dict[ref_accession] = gnm_size
                    gnm_quality_metric = gnm_completeness - (5 * gnm_contamination)
                    gnm_metric_dict[ref_accession] = float("{0:.3f}".format(gnm_quality_metric))

                    # get taxon at sampling rank
                    sampling_rank_taxon = ''
                    for each_rank_taxon in gtdb_taxonomy.split(';'):
                        if each_rank_taxon[:3] == '%s__' % sampling_rank:
                            sampling_rank_taxon = each_rank_taxon

                    if sampling_rank_taxon not in dod:
                        dod[sampling_rank_taxon] = dict()

                    dod[sampling_rank_taxon][ref_accession] = gnm_quality_metric

    selected_gnm_id_txt_handle = open(selected_gnm_id_txt, 'w')
    selected_gnm_meta_txt_handle = open(selected_gnm_meta_txt, 'w')
    selected_gnm_meta_txt_handle.write('Rank\tGenome\tCompleteness\tContamination\tSelection_criterion(completeness-5*contamination)\tSize(bp)\tTaxonomy\n')
    selected_gnm_list = []
    for each_sampling_taxon in dod:
        current_taxon_gnm_quality_dict = dod[each_sampling_taxon]
        gnms_sorted_by_quality = sorted(current_taxon_gnm_quality_dict.items(), key=lambda x: x[1], reverse=True)
        selected_gnm = 0
        for each_gnm in gnms_sorted_by_quality:
            if selected_gnm < gnm_num_per_taxon:
                gnm_id = each_gnm[0]
                selected_gnm_list.append(gnm_id)
                str_to_write = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (each_sampling_taxon, gnm_id,
                                                               gnm_cpl_dict[gnm_id],
                                                               gnm_ctm_dict[gnm_id],
                                                               gnm_metric_dict[gnm_id],
                                                               gnm_size_dict[gnm_id],
                                                               gnm_to_taxon_dict[gnm_id])
                selected_gnm_meta_txt_handle.write(str_to_write + '\n')
                selected_gnm_id_txt_handle.write(gnm_id + '\n')
                selected_gnm += 1
    selected_gnm_id_txt_handle.close()
    selected_gnm_meta_txt_handle.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p',           required=False, default='Genomes',          help='output prefix')
    parser.add_argument('-meta',        required=True,                              help='GTDB reference genome metadata')
    parser.add_argument('-taxon',       required=True,                              help='interested taxon')
    parser.add_argument('-r',           required=True,                              help='sampling rank, select from p, c, o, f, g and s')
    parser.add_argument('-n',           required=False, default=1, type=int,        help='numer of genome to retain per taxon')
    parser.add_argument('-cpl',         required=False, default=None, type=float,   help='completeness cutoff (0-100), default: None')
    parser.add_argument('-ctm',         required=False, default=None, type=float,   help='contamination cutoff, default: None')
    parser.add_argument('-ts',          required=False, action='store_true',        help='only consider type strain')
    args = vars(parser.parse_args())
    sampling_GTDB_gnms(args)
