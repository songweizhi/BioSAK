import os

assembly_id_txt = '/Users/songweizhi/Desktop/Mussel/ncbi_dataset.tsv'
op_dir          = '/Users/songweizhi/Desktop/test'

for assembly_id in open(assembly_id_txt):
    assembly_id = assembly_id.strip()
    print(assembly_id)
    op_metadata_tsv   = '%s/%s_data_report.tsv'   % (op_dir, assembly_id)
    ncbi_datasets_cmd = 'datasets summary genome accession %s --as-json-lines | dataformat tsv genome --fields assminfo-biosample-accession > %s' % (assembly_id, op_metadata_tsv)
    ncbi_datasets_cmd = 'datasets summary genome accession %s --as-json-lines | dataformat tsv genome > %s' % (assembly_id, op_metadata_tsv)
    os.system(ncbi_datasets_cmd)

    # biosample_accession_set = set()
    # line_num_index = 0
    # for each_line in open(op_metadata_tsv):
    #     if
    #     line_num_index += 0




'''

cd /Users/songweizhi/Desktop/test
dataformat tsv genome --inputfile GCA_000297895.2_metadata.txt




cd /Users/songweizhi/Desktop/test
dataformat tsv genome --inputfile GCA_007844125.1_data_report.jsonl

'''