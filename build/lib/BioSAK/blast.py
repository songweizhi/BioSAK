import argparse

blast_usage = '''
========================= blast example commands =========================

BioSAK blast -i Alignment.txt -o Alignment_formatted.txt -n 20
BioSAK blast -i Alignment.txt -o Alignment_formatted.txt -n 20 -iden 90

# The input file
On the top of the blast result page (the RID row), click "Download All", 
then choose "Text" from the top of the pop-up.

==========================================================================
'''


def blast(args):

    blast_op_txt     = args['i']
    blast_op_txt_fmt = args['o']
    top_hit_num      = args['n']
    iden_cutoff      = args['iden']
    cov_cutoff       = args['cov']

    query_to_hits_dict = dict()
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
            if wrote_line_num <= top_hit_num:
                if each_line_split == ['Description', 'Name', 'Name', 'Taxid', 'Score', 'Score', 'cover', 'Value', 'Ident', 'Len', 'Accession']:
                    blast_op_txt_fmt_handle.write('%s\t%s\t%s\n' % (current_query, each_line.strip(), 'Label'))
                else:
                    iden      = float(each_line_split[-3])
                    cov       = float(each_line_split[-5][:-1])
                    hit_genus = each_line_split[0]
                    if (iden >= iden_cutoff) and (cov >= cov_cutoff):
                        blast_op_txt_fmt_handle.write('%s\t%s\t%s__%s(%s)\n' % (current_query, each_line.strip(), current_query, hit_genus, iden))
                        wrote_line_num += 1

            # add to dict
            if each_line_split != ['Description', 'Name', 'Name', 'Taxid', 'Score', 'Score', 'cover', 'Value', 'Ident', 'Len', 'Accession']:
                if current_query not in query_to_hits_dict:
                    query_to_hits_dict[current_query] = []
                if len(query_to_hits_dict[current_query]) < top_hit_num:
                    query_to_hits_dict[current_query].append(each_line.strip())
    blast_op_txt_fmt_handle.close()


if __name__ == '__main__':

    blast_parser = argparse.ArgumentParser(usage=blast_usage)
    blast_parser.add_argument('-i',      required=True,                          help='input file')
    blast_parser.add_argument('-o',      required=True,                          help='output file')
    blast_parser.add_argument('-n',      required=False,type=int, default=9999,  help='top hits to keep, default is all')
    blast_parser.add_argument('-iden',   required=False,type=float, default=80,  help='minimum identity, default is 80')
    blast_parser.add_argument('-cov',    required=False,type=float, default=50,  help='minimum coverage, default is 50')
    args = vars(blast_parser.parse_args())
    blast(args)
