import os
import argparse
import certifi
os.environ["SSL_CERT_FILE"] = certifi.where()
os.environ["REQUESTS_CA_BUNDLE"] = certifi.where()
from ete3 import NCBITaxa


parse_mmseqs_tsv_usage = '''
============ parse_mmseqs_tsv example command ============

BioSAK parse_mmseqs_tsv -i taxonomyResult.tsv -o out.txt

==========================================================
'''


def parse_mmseqs_tsv(args):

    mmseqs_tsv_in   = args['i']
    table_out       = args['o']

    ncbi = NCBITaxa()

    table_out_handle = open(table_out, 'w')
    for each_line in open(mmseqs_tsv_in):
        each_line_split = each_line.strip().split('\t')
        seq_id = each_line_split[0]
        tax_id = each_line_split[1]

        name_str = 'na'
        if tax_id not in ['0']:
            lineage = ncbi.get_lineage(int(tax_id))
            name_dict = ncbi.get_taxid_translator(lineage)
            name_list = [name_dict[i] for i in lineage]
            name_str = ','.join(name_list)

        table_out_handle.write('%s\t%s\n' % (seq_id, name_str))
    table_out_handle.close()


if __name__ == '__main__':

    parse_mmseqs_tsv_parser = argparse.ArgumentParser(usage=parse_mmseqs_tsv_usage)
    parse_mmseqs_tsv_parser.add_argument('-i',  required=True, help='tsv file produced by "mmseqs createtsv"')
    parse_mmseqs_tsv_parser.add_argument('-o',  required=True, help='output table')
    args = vars(parse_mmseqs_tsv_parser.parse_args())
    parse_mmseqs_tsv(args)
