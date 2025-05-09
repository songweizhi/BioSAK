import os
import argparse


stats_blast_hit_usage = '''
================================= stats_blast_hit example commands =================================

BioSAK stats_blast_hit -i otu_vs_nt.txt -o stats.txt -node nodes.dmp -db RefSeq_taxonomy.txt -n 10

python3 /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/stats_blast_hit.py -i /Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/Unclassified_3696_vs_nt.txt -o /Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/Unclassified_3696_vs_nt_stats.txt -n 20 -node /Users/songweizhi/DB/taxdump_20250321/nodes.dmp -db /Users/songweizhi/DB/NCBI/RefSeq_taxonomy.txt

====================================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def group_tax_by_domain(nodes_dmp):
    child_to_parent_dict = dict()
    for each_node in open(nodes_dmp):
        each_node_split = each_node.strip().split('|')
        tax_id = int(each_node_split[0].strip())
        parent_tax_id = int(each_node_split[1].strip())
        child_to_parent_dict[tax_id] = parent_tax_id

    tax_id_set_bac = set()
    tax_id_set_ar = set()
    tax_id_set_eu = set()
    tax_id_set_virus = set()
    for tax_id in child_to_parent_dict.keys():
        # get tax_lineage_list
        tax_lineage_list = [tax_id]
        while tax_id != 1:
            tax_id = child_to_parent_dict[tax_id]
            tax_lineage_list.append(tax_id)
        # Bacteria
        if 2 in tax_lineage_list:
            for each_id in tax_lineage_list[:-3]:
                tax_id_set_bac.add(each_id)
        # Archaea
        elif 2157 in tax_lineage_list:
            for each_id in tax_lineage_list[:-3]:
                tax_id_set_ar.add(each_id)
        # Eukaryota
        elif 2759 in tax_lineage_list:
            for each_id in tax_lineage_list[:-3]:
                tax_id_set_eu.add(each_id)
        # Virus
        elif 10239 in tax_lineage_list:
            for each_id in tax_lineage_list[:-3]:
                tax_id_set_virus.add(each_id)
    return tax_id_set_bac, tax_id_set_ar, tax_id_set_eu, tax_id_set_virus


def stats_blast_hit(args):

    blast_op_txt    = args['i']
    hit_stats_txt   = args['o']
    nodes_dmp       = args['node']
    refseq_tax_txt  = args['db']
    hit_num_to_keep = args['n']
    pct_cutoff      = args['c']

    f_name, f_path, f_base, f_ext = sep_path_basename_ext(hit_stats_txt)
    cmd_txt = '%s/%s_RefSeq_classification.sh' % (f_path, f_base)
    tax_set_bac, tax_set_ar, tax_set_eu, tax_set_virus = group_tax_by_domain(nodes_dmp)

    refseq_tax_dict = dict()
    for each_ref in open(refseq_tax_txt):
        each_ref_split = each_ref.strip().split('\t')
        if len(each_ref_split) == 2:
            refseq_tax_dict[each_ref_split[0]] = each_ref_split[1]

    query_to_hit_dict = dict()
    for each_line in open(blast_op_txt):
        each_line_split = each_line.strip().split('\t')
        query_id = each_line_split[0]
        subject_id = each_line_split[1]
        if query_id not in query_to_hit_dict:
            query_to_hit_dict[query_id] = set()
        if len(query_to_hit_dict[query_id]) < hit_num_to_keep:
            query_to_hit_dict[query_id].add(subject_id)

    hit_stats_txt_handle = open(hit_stats_txt, 'w')
    hit_stats_txt_handle.write('Query\tArchaea\tBacteria\tEukaryote\tVirus\tRest\tUnknown\tClassification\n')
    unclassified_refseq_set = set()
    for query in query_to_hit_dict:
        hit_tax_list = []
        for each_hit in query_to_hit_dict[query]:
            hit_tax = refseq_tax_dict.get(each_hit, 'na')
            if hit_tax == 'na':
                hit_tax_list.append(hit_tax)
                unclassified_refseq_set.add(each_hit)
            else:
                hit_tax_list.append(int(hit_tax))

        hit_num_bac = 0
        hit_num_ar = 0
        hit_num_eu = 0
        hit_num_virus = 0
        hit_num_na = 0
        hit_num_rest = 0
        for hit_tax in hit_tax_list:
            if hit_tax in tax_set_bac:
                hit_num_bac += 1
            elif hit_tax in tax_set_ar:
                hit_num_ar += 1
            elif hit_tax in tax_set_eu:
                hit_num_eu += 1
            elif hit_tax in tax_set_virus:
                hit_num_virus += 1
            elif hit_tax == 'na':
                hit_num_na += 1
            else:
                hit_num_rest += 1

        total_hit_num = len(hit_tax_list)
        try:
            hit_pct_bac = int(hit_num_bac*100/total_hit_num)
        except:
            hit_pct_bac = float("{0:.2f}".format(hit_num_bac*100/total_hit_num))
        try:
            hit_pct_ar = int(hit_num_ar*100/total_hit_num)
        except:
            hit_pct_ar = float("{0:.2f}".format(hit_num_ar*100/total_hit_num))
        try:
            hit_pct_eu = int(hit_num_eu*100/total_hit_num)
        except:
            hit_pct_eu = float("{0:.2f}".format(hit_num_eu*100/total_hit_num))
        try:
            hit_pct_virus = int(hit_num_virus*100/total_hit_num)
        except:
            hit_pct_virus = float("{0:.2f}".format(hit_num_virus*100/total_hit_num))
        try:
            hit_pct_na = int(hit_num_na*100/total_hit_num)
        except:
            hit_pct_na = float("{0:.2f}".format(hit_num_na*100/total_hit_num))
        try:
            hit_pct_rest = int(hit_num_rest*100/total_hit_num)
        except:
            hit_pct_rest = float("{0:.2f}".format(hit_num_rest*100/total_hit_num))

        classification = 'Unclassified'
        if hit_pct_ar >= pct_cutoff:
            classification = 'd__Archaea;p__;c__;o__;f__;g__;s__'
        elif hit_pct_bac >= pct_cutoff:
            classification = 'd__Bacteria;p__;c__;o__;f__;g__;s__'
        elif hit_pct_eu >= pct_cutoff:
            classification = 'Eukaryote'
        elif hit_pct_virus >= pct_cutoff:
            classification = 'Virus'
        elif hit_pct_rest >= pct_cutoff:
            classification = 'Rest'
        elif hit_pct_na >= pct_cutoff:
            classification = 'Unknown'

        hit_stats_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (query, hit_pct_ar, hit_pct_bac, hit_pct_eu, hit_pct_virus, hit_pct_rest, hit_pct_na, classification))
        #hit_stats_txt_handle.write('%s\t%s\n' % (query, classification))
    hit_stats_txt_handle.close()
    print('Done!')


if __name__ == '__main__':

    stats_blast_hit_parser = argparse.ArgumentParser(usage=stats_blast_hit_usage)
    stats_blast_hit_parser.add_argument('-i',    required=True,                         help='blast result, outfmt need to be 6')
    stats_blast_hit_parser.add_argument('-o',    required=True,                         help='output stats')
    stats_blast_hit_parser.add_argument('-node', required=True,                         help='nodes.dmp ')
    stats_blast_hit_parser.add_argument('-db',   required=True,                         help='RefSeq taxonomy file')
    stats_blast_hit_parser.add_argument('-c',    required=False, type=int, default=90,  help='cutoff, default is 90')
    stats_blast_hit_parser.add_argument('-n',    required=False, type=int, default=20,  help='number of top hits to consider, default is 20')
    args = vars(stats_blast_hit_parser.parse_args())
    stats_blast_hit(args)
