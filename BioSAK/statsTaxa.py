import argparse


statsTaxa_usage = '''
============ statsTaxa example commands ============

BioSAK statsTaxa -i taxa.txt -o taxa_stats.txt

# Example of input file (no header, GTDB format):
d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__JAHRAL01;s__JAHRAL01 sp003724275
d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__WTHC01;s__WTHC01 sp026708305
d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__TA-20;s__TA-20 sp013287585
d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__Nitrosopumilus;s__

====================================================
'''


def statsTaxa(args):

    taxa_txt        = args['i']
    taxa_stats_txt  = args['o']

    tax_rank_order_list = []
    stats_dod = dict()
    for each_line in open(taxa_txt):
        each_line_split = each_line.strip().split(';')
        for each_r in each_line_split:
            tax_rank = each_r.split('__')[0]

            if tax_rank not in tax_rank_order_list:
                tax_rank_order_list.append(tax_rank)

            if tax_rank not in stats_dod:
                stats_dod[tax_rank] = dict()

            if each_r not in stats_dod[tax_rank]:
                stats_dod[tax_rank][each_r] = 1
            else:
                stats_dod[tax_rank][each_r] += 1

    taxa_stats_txt_handle = open(taxa_stats_txt, 'w')
    for each_rank in tax_rank_order_list:
        current_rank_stats_dict = stats_dod[each_rank]
        tax_sorted = [i[0] for i in sorted(current_rank_stats_dict.items(), key=lambda x: x[1])][::-1]
        for each_tax in tax_sorted:
            taxa_stats_txt_handle.write('%s\t%s\t%s\n' % (each_rank, each_tax, current_rank_stats_dict[each_tax]))
    taxa_stats_txt_handle.close()

    print('Done!')


if __name__ == '__main__':

    statsTaxa_parser = argparse.ArgumentParser(usage=statsTaxa_usage)
    statsTaxa_parser.add_argument('-i', required=True, help='input txt')
    statsTaxa_parser.add_argument('-o', required=True, help='output txt')
    args = vars(statsTaxa_parser.parse_args())
    statsTaxa(args)
