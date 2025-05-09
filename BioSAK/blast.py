import argparse

blast_usage = '''
=============================== blast example commands ===============================

BioSAK blast -i Alignment.txt -o Alignment_formatted.txt -n 20
BioSAK blast -i Alignment.txt -o Alignment_formatted.txt -n 20 -iden 90
BioSAK blast -i Alignment.txt -o Alignment_formatted.txt -n 20 -iden 90 -tax tax.txt

# The input file
On the top of the blast result page (the RID row), click "Download All", 
then choose "Text" from the top of the pop-up.

# Use "BioSAK taxdump" module to generate the "-tax" file.
-tax /Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt

======================================================================================
'''


def blast(args):

    blast_op_txt     = args['i']
    blast_op_txt_fmt = args['o']
    top_hit_num      = args['n']
    iden_cutoff      = args['iden']
    cov_cutoff       = args['cov']
    tax_file         = args['tax']

    blast_op_txt_fmt_id_txt = '.'.join(blast_op_txt_fmt.split('.')[:-1]) + '_ref_accession.txt'
    ref_to_taxon_txt        = '.'.join(blast_op_txt_fmt.split('.')[:-1]) + '_ref_taxon.txt'

    tax_lineage_dict = dict()
    if tax_file is not None:
        for each_line in open(tax_file):
            each_line_split = each_line.strip().split('\t')
            tax_lineage_dict[each_line_split[0]] = each_line_split[1]

    accession_set = set()
    accession_tax_dict = dict()
    current_query = ''
    current_query_len = 0
    keep_line = 0
    wrote_line_num = 0
    blast_op_txt_fmt_handle = open(blast_op_txt_fmt, 'w')
    for each_line in open(blast_op_txt):
        each_line_split = each_line.strip().split()

        # get current_query and current_query_len
        if (each_line.startswith('Query #')) and ('Query ID:' in each_line):
            current_query = each_line.strip().split(' Query ID: ')[0].split(':')[1].strip()
            current_query_len = each_line.strip().split(' Length: ')[1]
            wrote_line_num = 0
            blast_op_txt_fmt_handle.write('\n')

        # decide to keep current line or not
        if each_line_split == ['Description', 'Name', 'Name', 'Taxid', 'Score', 'Score', 'cover', 'Value', 'Ident', 'Len', 'Accession']:
            keep_line = 1
        if len(each_line.strip()) == 0:
            keep_line = 0

        # write out
        if keep_line == 1:
            if wrote_line_num < top_hit_num:
                if each_line_split == ['Description', 'Name', 'Name', 'Taxid', 'Score', 'Score', 'cover', 'Value', 'Ident', 'Len', 'Accession']:
                    blast_op_txt_fmt_handle.write('%s\t%s\t%s\n' % (current_query, each_line.strip(), 'Label'))
                else:
                    iden      = float(each_line_split[-3])
                    cov       = float(each_line_split[-5][:-1])
                    hit_genus = each_line_split[0]
                    accession = each_line_split[-1]
                    hit_tax_lineage= tax_lineage_dict.get(hit_genus, hit_genus)
                    if (iden >= iden_cutoff) and (cov >= cov_cutoff):
                        blast_op_txt_fmt_handle.write('%s\t%s\t%s__%s(%s)\n' % (current_query, each_line.strip(), current_query, hit_tax_lineage, iden))
                        accession_set.add(accession)
                        accession_tax_dict[accession] = hit_tax_lineage
                        wrote_line_num += 1
    blast_op_txt_fmt_handle.close()

    blast_op_txt_fmt_id_txt_handle = open(blast_op_txt_fmt_id_txt, 'w')
    ref_to_taxon_txt_handle = open(ref_to_taxon_txt, 'w')
    for each_accession in sorted(list(accession_set)):
        accession_tax = accession_tax_dict.get(each_accession, 'na')
        blast_op_txt_fmt_id_txt_handle.write('%s\n' % each_accession)
        ref_to_taxon_txt_handle.write('%s\t%s\n' % (each_accession, accession_tax))
    blast_op_txt_fmt_id_txt_handle.close()
    ref_to_taxon_txt_handle.close()

    # # to download the sequence, run:
    # download_seq_cmd = 'cat %s | paste -sd "," - | xargs -I{} esearch -db nucleotide -query {} | efetch -format fasta > %s' % (blast_op_txt_fmt_id_txt, blast_op_txt_fmt_id_txt.replace('.txt', '.fasta'))
    # print(download_seq_cmd)


if __name__ == '__main__':

    blast_parser = argparse.ArgumentParser(usage=blast_usage)
    blast_parser.add_argument('-i',      required=True,                          help='input file')
    blast_parser.add_argument('-o',      required=True,                          help='output file')
    blast_parser.add_argument('-n',      required=False,type=int, default=9999,  help='top hits to keep, default is all')
    blast_parser.add_argument('-iden',   required=False,type=float, default=80,  help='minimum identity, default is 80')
    blast_parser.add_argument('-cov',    required=False,type=float, default=50,  help='minimum coverage, default is 50')
    blast_parser.add_argument('-tax',    required=False, default=None,           help='NCBI taxonomy file')
    args = vars(blast_parser.parse_args())
    blast(args)


'''

cat 28S_ZPVK91YC016-Alignment_reformatted_ref_accession.txt | paste -sd "," - | xargs -I{} esearch -db nucleotide -query {} | efetch -format fasta > 28S_ZPVK91YC016-Alignment_reformatted_ref_accession.fasta


'''