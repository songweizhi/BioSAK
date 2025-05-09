import os
import argparse


get_eu_otu_usage = '''
================================= get_eu_otu example commands =================================

BioSAK get_eu_otu -i otu_vs_nt.txt -o EU_OTU.txt -node nodes.dmp -db RefSeq_taxonomy.txt -n 10

# An OTU was considered as eukaryotic if more than half of its hits are from Eukaryotes.

===============================================================================================
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


def get_eu_otu(args):

    blast_op_txt            = args['i']
    eu_otu_txt              = args['o']
    nodes_dmp               = args['node']
    refseq_tax_txt          = args['db']
    hit_num_to_keep         = args['n']
    consider_as_eu_cutoff   = args['c']

    f_name, f_path, f_base, f_ext = sep_path_basename_ext(eu_otu_txt)
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

    unclassified_refseq_set = set()
    eu_otu_set = set()
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

        hit_num_without_na = hit_num_bac + hit_num_ar + hit_num_eu + hit_num_virus

        # An OTU was considered as eukaryotic if more than half of the hits are from Eukaryotes
        eu_hit_pct = 0
        if hit_num_without_na > 0:
            eu_hit_pct = hit_num_eu / hit_num_without_na
        if eu_hit_pct > consider_as_eu_cutoff:
            eu_otu_set.add(query)

    # write out EU OTU
    eu_otu_txt_handle = open(eu_otu_txt, 'w')
    eu_otu_txt_handle.write('\n'.join(sorted(list(eu_otu_set))))
    eu_otu_txt_handle.close()

    # report
    print('Note:')
    if len(unclassified_refseq_set) > 0:
        cmd_txt_handle = open(cmd_txt, 'w')
        for seq_accession in unclassified_refseq_set:
            cmd = 'esearch -db nucleotide -query %s | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId > %s.txt' % (seq_accession, seq_accession)
            cmd_txt_handle.write(cmd + '\n')
        cmd_txt_handle.close()
        print('Taxon ID for %s RefSeqs is missing, to increase accuracy, add their classifications to %s and rerun this script.' % (len(unclassified_refseq_set), refseq_tax_txt))
        print('RefSeq classification commands exported to: %s' % cmd_txt)
    print('The ID of %s eukaryotic OTUs has exported to: %s' % (len(eu_otu_set), eu_otu_txt))
    print('Done!')


if __name__ == '__main__':

    get_eu_otu_parser = argparse.ArgumentParser(usage=get_eu_otu_usage)
    get_eu_otu_parser.add_argument('-i',    required=True,                              help='blast result, outfmt need to be 6')
    get_eu_otu_parser.add_argument('-o',    required=True,                              help='output of identified eukaryotic OTUs')
    get_eu_otu_parser.add_argument('-node', required=True,                              help='nodes.dmp ')
    get_eu_otu_parser.add_argument('-db',   required=True,                              help='RefSeq taxonomy file')
    get_eu_otu_parser.add_argument('-c',    required=False, type=float, default=0.3,    help='consider as eukaryotic if the proportion of eukaryotic hits higher than this cutoff, default is 0.3')
    get_eu_otu_parser.add_argument('-n',    required=False, type=int, default=10,       help='number of top hits to consider, default is 10')
    args = vars(get_eu_otu_parser.parse_args())
    get_eu_otu(args)
