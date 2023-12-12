
gtdb_gnm_metadata = 'gtdb_gnm_metadata.txt'

col_index = {}
line_num_index = 0
for each_ref in open(gtdb_gnm_metadata):
    each_ref_split = each_ref.strip().split('\t')
    if line_num_index == 0:
        col_index = {key: i for i, key in enumerate(each_ref_split)}
        line_num_index += 1
    else:
        gnm_size = each_ref_split[col_index['genome_size']]
        gtdb_taxonomy = each_ref_split[col_index['gtdb_taxonomy']]

