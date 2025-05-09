import os
import glob
import argparse


GenBank2_usage = '''
========================= GenBank2 example commands =========================

BioSAK GenBank2 -i accession.txt -o op_dir -seq
BioSAK GenBank2 -i accession.txt -o op_dir -voucher
BioSAK GenBank2 -i accession.txt -o op_dir -organism -tax ncbi_taxonomy.txt

=============================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def combine_esearch_op_organism(op_dir, tax_lineage_dict, op_txt):

    file_re = '%s/*_organism.txt' % op_dir
    file_list = glob.glob(file_re)

    op_txt_handle = open(op_txt, 'w')
    for each_file in file_list:
        f_name, f_path, f_base, f_ext = sep_path_basename_ext(each_file)
        accession_id = f_base.split('_organism')[0]

        organism_info = ''
        with open(each_file) as f:
            organism_info = f.readline().replace('ORGANISM', '').strip()
            organism_g = organism_info.split()[0]
            g_full_lineage = tax_lineage_dict.get(organism_g, ('g__' + organism_g))
            full_lineage_str = '%s;s__%s' % (g_full_lineage, organism_info)
        if len(tax_lineage_dict) == 0:
            op_txt_handle.write('%s\t%s\n' % (accession_id, organism_info))
        else:
            op_txt_handle.write('%s\t%s\n' % (accession_id, full_lineage_str))
    op_txt_handle.close()


def combine_esearch_op_voucher(op_dir, op_txt):

    file_re = '%s/*_voucher.txt' % op_dir
    file_list = glob.glob(file_re)

    op_txt_handle = open(op_txt, 'w')
    for each_file in file_list:
        f_name, f_path, f_base, f_ext = sep_path_basename_ext(each_file)
        accession_id = f_base.split('_voucher')[0]

        voucher_info = ''
        with open(each_file) as f:
            #voucher_info = f.readline().replace('ORGANISM', '').strip()
            voucher_info = f.readline().strip()
            if 'specimen_voucher=' in voucher_info:
                voucher_info = voucher_info.replace('/specimen_voucher="', '')[:-1]
        if voucher_info != '':
            op_txt_handle.write('%s\t%s\n' % (accession_id, voucher_info))
    op_txt_handle.close()


def GenBank2(args):

    accession_txt       = args['i']
    op_dir              = args['o']
    get_sequence        = args['seq']
    get_organism_info   = args['organism']
    get_voucher_info    = args['voucher']
    tax_file            = args['tax']

    tmp_dir                     = '%s/tmp'                      % op_dir
    combined_fa                 = '%s/accession_sequence.fasta' % op_dir
    combined_organism_info_txt  = '%s/accession_organism.txt'   % op_dir
    combined_voucher_info_txt   = '%s/accession_voucher.txt'    % op_dir

    tax_lineage_dict = dict()
    if tax_file is not None:
        for each_line in open(tax_file):
            each_line_split = each_line.strip().split('\t')
            tax_lineage_dict[each_line_split[0]] = each_line_split[1]

    if get_sequence is True:
        os.system('cat %s/*.fasta > %s' % (tmp_dir, combined_fa))

    if get_organism_info is True:
        combine_esearch_op_organism(tmp_dir, tax_lineage_dict, combined_organism_info_txt)

    if get_voucher_info is True:
        combine_esearch_op_voucher(tmp_dir, combined_voucher_info_txt)


if __name__ == '__main__':

    GenBank2_parser = argparse.ArgumentParser(usage=GenBank2_usage)
    GenBank2_parser.add_argument('-i',          required=True,                       help='input txt containing accession id')
    GenBank2_parser.add_argument('-o',          required=True,                       help='output dir')
    GenBank2_parser.add_argument('-seq',        required=False, action="store_true", help='get sequences')
    GenBank2_parser.add_argument('-organism',   required=False, action="store_true", help='get organism info')
    GenBank2_parser.add_argument('-voucher',    required=False, action="store_true", help='get voucher info')
    GenBank2_parser.add_argument('-tax',        required=False, default=None,        help='NCBI taxonomy file, specify to get the full taxon lineage')
    args = vars(GenBank2_parser.parse_args())
    GenBank2(args)
