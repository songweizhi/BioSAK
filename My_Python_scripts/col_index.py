
gtdb_gnm_metadata = 'gtdb_gnm_metadata.txt'

col_index = {}
for each_ref in open(gtdb_gnm_metadata):
    each_ref_split = each_ref.strip().split('\t')
    if each_ref.startswith('accession	ambiguous_bases'):
        col_index = {key: i for i, key in enumerate(each_ref_split)}
    else:
        gnm_size = each_ref_split[col_index['genome_size']]
        gtdb_taxonomy = each_ref_split[col_index['gtdb_taxonomy']]
