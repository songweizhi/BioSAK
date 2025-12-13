import argparse


metadata_usage = '''
======================================= metadata example commands =======================================

BioSAK metadata -id_col 1 -i meta.txt -o meta_new.txt -a checkm_quality.tsv -c 'Completeness,Contamination'
BioSAK metadata -id_col 2 -i meta.txt -o meta_new.txt -a gtdb_taxonomy.tsv -c classification

=========================================================================================================
'''


def metadata(args):

    metadata_txt_in       = args['i']
    new_attribute_file    = args['a']
    interested_col_header = args['c']
    metadata_txt_out      = args['o']
    metadata_txt_sep      = args['ms']
    new_attribute_sep     = args['as']
    id_column             = args['id_col']

    interested_col_list = interested_col_header.split(',')
    default_key_str = metadata_txt_sep.join(['na']*len(interested_col_list))

    # read in new attribute
    new_attribute_dict = dict()
    col_index = dict()
    line_num_index = 0
    for each_line in open(new_attribute_file):
        line_split = each_line.strip().split(new_attribute_sep)
        if line_num_index == 0:
            col_index = {key: i for i, key in enumerate(line_split)}
        else:
            sample_id = line_split[0]
            value_list = []
            for interested_col in interested_col_list:
                value_list.append(line_split[col_index[interested_col]])
            value_str = metadata_txt_sep.join(value_list)
            new_attribute_dict[sample_id] = value_str
        line_num_index += 1

    # add new attribute to metadata file
    metadata_txt_out_handle = open(metadata_txt_out, 'w')
    line_num_index = 0
    for each_line in open(metadata_txt_in):
        each_line_split = each_line.strip().split(metadata_txt_sep)
        sample_id = each_line_split[id_column - 1]
        if line_num_index == 0:
            metadata_txt_out_handle.write('%s%s%s\n' % (each_line.strip(), metadata_txt_sep, metadata_txt_sep.join(interested_col_list)))
        else:
            metadata_txt_out_handle.write('%s%s%s\n' % (each_line.strip(), metadata_txt_sep, new_attribute_dict.get(sample_id, default_key_str)))
        line_num_index += 1
    metadata_txt_out_handle.close()


if __name__ == '__main__':

    metadata_parser = argparse.ArgumentParser(usage=metadata_usage)
    metadata_parser.add_argument('-i',  required=True,                   help='input metadata file')
    metadata_parser.add_argument('-o',  required=True,                   help='output metadata file')
    metadata_parser.add_argument('-a',  required=True,                   help='new attribute file')
    metadata_parser.add_argument('-c',  required=True,                   help='interested column header')
    metadata_parser.add_argument('-ms', required=False, default='\t',    help='metadata column separator')
    metadata_parser.add_argument('-as', required=False, default='\t',    help='new attribute file column separator')
    metadata_parser.add_argument('-id_col', required=False, type=int, default=1,    help='id column in input metadata file, default is 1 (the first column)')
    args = vars(metadata_parser.parse_args())
    metadata(args)
