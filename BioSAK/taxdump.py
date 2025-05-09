import argparse

taxdump_usage = '''
====================== taxdump example commands ======================

BioSAK taxdump -node nodes.dmp -name names.dmp -o ncbi_taxonomy.txt

# Input files available at:
https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

======================================================================
'''


def taxdump(args):

    nodes_dmp   = args['node']
    names_dmp   = args['name']
    op_txt      = args['o']

    tax_rank_abbrev_dict = {'domain' : 'd',
                            'phylum' : 'p', 'subphylum' : 'sp', 'superphylum' : 'up',
                            'class'  : 'c', 'subclass'  : 'sc', 'superclass'  : 'uc',
                            'order'  : 'o', 'suborder'  : 'so', 'superorder'  : 'uo',
                            'family' : 'f', 'subfamily' : 'sf', 'superfamily' : 'uf',
                            'genus'  : 'g', 'subgenus'  : 'sg', 'supergenus'  : 'ug',
                            'species': 's', 'subspecies': 'ss', 'superspecies': 'us',
                            'kingdom': 'k', 'subkingdom': 'sk', 'superkingdom': 'uk'}

    id_to_name_dict = dict()
    for each_node in open(names_dmp):
        each_node_split = each_node.strip().split('|')
        node_id    = int(each_node_split[0].strip())
        node_name  = each_node_split[1].strip()
        name_class = each_node_split[3].strip()
        if name_class == 'scientific name':
            id_to_name_dict[node_id] = node_name

    tax_set = set()
    tax_rank_dict = dict()
    child_to_parent_dict = dict()
    for each_node in open(nodes_dmp):
        each_node_split = each_node.strip().split('|')
        tax_id = int(each_node_split[0].strip())
        parent_tax_id = int(each_node_split[1].strip())
        tax_rank = each_node_split[2].strip()
        child_to_parent_dict[tax_id] = parent_tax_id
        tax_rank_dict[tax_id] = tax_rank
        if tax_rank not in ['species', 'species group', 'subspecies', 'no rank', 'strain']:
            tax_set.add(tax_id)

    op_txt_handle = open(op_txt, 'w')
    for tax_id in sorted(list(tax_set)):
        current_tax_name = id_to_name_dict[tax_id]
        # get tax_lineage_list
        tax_lineage_list = [tax_id]
        while tax_id != 1:
            tax_id = child_to_parent_dict[tax_id]
            tax_lineage_list.append(tax_id)

        # get tax_lineage_list_by_name
        tax_lineage_list_by_name = []
        for each_id in tax_lineage_list:
            tax_name = id_to_name_dict[each_id]
            tax_rank = tax_rank_dict[each_id]
            tax_rank_abbrev = tax_rank_abbrev_dict.get(tax_rank, tax_rank)
            if tax_rank in tax_rank_abbrev_dict:
                tax_lineage_list_by_name.append('%s__%s' % (tax_rank_abbrev, tax_name))

        tax_lineage_str = ';'.join(tax_lineage_list_by_name[::-1])

        # add g__ if missing
        if 'g__' not in tax_lineage_str:
            tax_lineage_str = tax_lineage_str + ';g__'

        # add f__ if missing
        if 'f__' not in tax_lineage_str:
            if 'ug__' in tax_lineage_str:
                tax_lineage_str = tax_lineage_str.replace(';ug__', ';f__;ug__')
            else:
                tax_lineage_str = tax_lineage_str.replace(';g__', ';f__;g__')

        # add o__ if missing
        if 'o__' not in tax_lineage_str:
            if 'uf__' in tax_lineage_str:
                tax_lineage_str = tax_lineage_str.replace(';uf__', ';o__;uf__')
            else:
                tax_lineage_str = tax_lineage_str.replace(';f__', ';o__;f__')

        # add c__ if missing
        if 'c__' not in tax_lineage_str:
            if 'uo__' in tax_lineage_str:
                tax_lineage_str = tax_lineage_str.replace(';uo__', ';c__;uo__')
            else:
                tax_lineage_str = tax_lineage_str.replace(';o__', ';c__;o__')

        # add p__ if missing
        if 'p__' not in tax_lineage_str:
            if 'uc__' in tax_lineage_str:
                tax_lineage_str = tax_lineage_str.replace(';uc__', ';p__;uc__')
            else:
                tax_lineage_str = tax_lineage_str.replace(';c__', ';p__;c__')

        # write out
        op_txt_handle.write('%s\t%s\n' % (current_tax_name, tax_lineage_str))
    op_txt_handle.close()


if __name__ == '__main__':

    taxdump_parser = argparse.ArgumentParser(usage=taxdump_usage)
    taxdump_parser.add_argument('-node',  required=True,    help='nodes.dmp')
    taxdump_parser.add_argument('-name',  required=True,    help='names.dmp')
    taxdump_parser.add_argument('-o',     required=True,    help='output file')
    args = vars(taxdump_parser.parse_args())
    taxdump(args)
