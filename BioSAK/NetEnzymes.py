import argparse
import networkx as nx
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime
from BioSAK.BioSAK_config import config_dict
from BioSAK.global_functions import is_number
from BioSAK.global_functions import time_format
from BioSAK.global_functions import sep_path_basename_ext


NetEnzymes_parser_usage = '''
=========================== NetEnzymes example commands ===========================

# get network of all enzymes in Ecoli_ec.txt
BioSAK NetEnzymes -ec Ecoli_ec.txt -to_skip skip.txt -NoHyphen

# get network of enzymes belong to ko 00010 (Glycolysis) in Ecoli_ec.txt
BioSAK NetEnzymes -ec Ecoli_ec.txt -ko 00010 -to_skip skip.txt -NoHyphen -plot 

# EC file format (one EC per line):
2.7.7.23
3.6.3.28

# to_skip file format (one substrate per line):
H2O
H+

===================================================================================
'''

'''
cd /Users/songweizhi/MetaCyc_demo
python3 ~/PycharmProjects/BioSAK/BioSAK/NetEnzymes.py -ec Ecoli_ec.txt -to_skip skip.txt -plot -NoHyphen
python3 ~/PycharmProjects/BioSAK/BioSAK/NetEnzymes.py -ec Ecoli_ec.txt -ko 00010 -to_skip skip.txt -plot -NoHyphen
python3 ~/PycharmProjects/BioSAK/BioSAK/NetEnzymes.py -ec Ecoli_ec.txt -ko 00020 -to_skip skip.txt -plot -NoHyphen
python3 ~/PycharmProjects/BioSAK/BioSAK/NetEnzymes.py -ec Glycolysis_ECs.txt -to_skip skip.txt -NoHyphen -plot 

python3 ~/PycharmProjects/BioSAK/BioSAK/NetEnzymes.py -enzymes TCA_enzymes.txt -compounds TCA_compounds.txt -plot
enzymes 

BioSAK NetEnzymes -ec Ecoli_ec.txt -ko 00020 -to_skip skip.txt -plot -NoHyphen
BioSAK NetEnzymes -ec Glycolysis_ECs.txt -to_skip skip.txt -NoHyphen -plot

To-do:
1. Provide a network for a list of enzymes and substrates
2. '=' symbol

'''


def get_ko2description_dict(ko00001_keg):

    As_description_dict = {}
    Bs_description_dict = {}
    Cs_description_dict = {}
    Ds_description_dict = {}
    D2ABCD_dict = {}
    ko2level_dict = {}

    current_A = ''
    current_B = ''
    current_C = ''
    for each_line in open(ko00001_keg):
        if each_line[0] in ['A', 'B', 'C', 'D']:
            each_line_split = each_line.strip().split(' ')

            if each_line[0] == 'A':
                current_A_id = each_line_split[0]
                current_A_description = ' '.join(each_line_split[1:-1])
                ko2level_dict[current_A_id] = 'A'
                current_A = current_A_id
                As_description_dict[current_A_id] = current_A_description

            elif each_line[0] == 'B':
                if len(each_line_split) > 1:
                    current_B_id = each_line_split[2]
                    current_B_description = ' '.join(each_line_split[3:-1])
                    ko2level_dict[current_B_id] = 'B'
                    current_B = current_B_id
                    Bs_description_dict[current_B_id] = current_B_description

            elif each_line[0] == 'C':
                current_C_id = each_line_split[4]
                current_C_description = ' '.join(each_line_split[5:-1])
                ko2level_dict[current_C_id] = 'C'
                current_C = current_C_id
                Cs_description_dict[current_C_id] = current_C_description

            elif each_line[0] == 'D':
                current_D_id = each_line_split[6]
                current_D_description = ' '.join(each_line_split[7:])
                ko2level_dict[current_D_id] = 'D'
                Ds_description_dict[current_D_id] = current_D_description
                ABCD_value = '%s|%s|%s|%s' % (current_A, current_B, current_C, current_D_id)
                if current_D_id not in D2ABCD_dict:
                    D2ABCD_dict[current_D_id] = [ABCD_value]
                elif (current_D_id in D2ABCD_dict) and (ABCD_value not in D2ABCD_dict[current_D_id]):
                    D2ABCD_dict[current_D_id].append(ABCD_value)

    return As_description_dict, Bs_description_dict, Cs_description_dict, Ds_description_dict, D2ABCD_dict, ko2level_dict


def get_ec_of_interested_ko(D2ABCD_dict, KO_description_D_dict, ko_level, ko_id):

    # get ec list
    interested_ko_ec_list = set()
    for each_ko_D in D2ABCD_dict:

        each_ko_D_description = KO_description_D_dict[each_ko_D]
        if '[EC:' in each_ko_D_description:
            ec_str = each_ko_D_description.strip().split('[EC:')[1][:-1]

            # get ec list
            ec_list = [ec_str]
            if ' ' in ec_str:
                ec_list = ec_str.split(' ')

            for each_group in D2ABCD_dict[each_ko_D]:
                each_ko_D_split = each_group.split('|')
                if (ko_level == 'A') and (each_ko_D_split[0] == ko_id):
                    for ec in ec_list:
                        interested_ko_ec_list.add(ec)

                elif (ko_level == 'B') and (each_ko_D_split[1] == ko_id):
                    for ec in ec_list:
                        interested_ko_ec_list.add(ec)

                elif (ko_level == 'C') and (each_ko_D_split[2] == ko_id):
                    for ec in ec_list:
                        interested_ko_ec_list.add(ec)

                elif (ko_level == 'D') and (each_ko_D_split[3] == ko_id):
                    for ec in ec_list:
                        interested_ko_ec_list.add(ec)

    return interested_ko_ec_list


def parse_biological_raction(G, reaction, compounds_to_exclude, compounds_to_include, node_color_dict):

    reaction_split = reaction.strip().split('\t')
    enzyme = reaction_split[0]
    reaction_equation = reaction_split[1]

    G.add_node(enzyme, data=True, shape='s', color_map=node_color_dict['enzyme'])

    substrates = ''
    products = ''
    if ' → ' in reaction_equation:
        reaction_equation_split = reaction_equation.split(' → ')
        substrates = reaction_equation_split[0]
        products = reaction_equation_split[1]
    elif ' ↔ ' in reaction_equation:
        reaction_equation_split = reaction_equation.split(' ↔ ')
        substrates = reaction_equation_split[0]
        products = reaction_equation_split[1]
    if ' ← ' in reaction_equation:
        reaction_equation_split = reaction_equation.split(' ← ')
        substrates = reaction_equation_split[1]
        products = reaction_equation_split[0]

    substrate_list = [substrates]
    if '+' in substrates:
        substrate_list = substrates.split(' + ')

    product_list = [products]
    if '+' in products:
        product_list = products.split(' + ')

    substrate_list_no_num = []
    for substrate in substrate_list:
       if (' ' in substrate):
           substrate_split = substrate.split(' ')
           if (is_number(substrate_split[0]) is True) or (substrate_split[0] in ['a', 'an', 'n']):
               substrate_list_no_num.append(' '.join(substrate_split[1:]))
           else:
               substrate_list_no_num.append(substrate)
       else:
           substrate_list_no_num.append(substrate)


    product_list_no_num = []
    for product in product_list:
        if (' ' in product):
            product_split = product.split(' ')
            if (is_number(product_split[0]) is True) or (product_split[0] in ['a', 'an', 'n']):
                product_list_no_num.append(' '.join(product_split[1:]))
            else:
                product_list_no_num.append(product)
        else:
            product_list_no_num.append(product)

    # remove brackets
    substrate_list_no_num_brackets = []
    for substrate_no_num in substrate_list_no_num:
        if '[' in substrate_no_num:

            substrate_no_num_split = substrate_no_num.split('[')
            substrate_no_num_brackets = substrate_no_num_split[0]
            substrate_list_no_num_brackets.append(substrate_no_num_brackets)
        else:
            substrate_list_no_num_brackets.append(substrate_no_num)

    #print(substrate_list_no_num)
    #print(substrate_list_no_num_brackets)

    product_list_no_num_brackets = []
    for product_no_num in product_list_no_num:
        if '[' in product_no_num:
            product_no_num_split = product_no_num.split('[')
            product_no_num_brackets = product_no_num_split[0]
            product_list_no_num_brackets.append(product_no_num_brackets)
        else:
            product_list_no_num_brackets.append(product_no_num)

    for substrate in substrate_list_no_num_brackets:
        if substrate not in compounds_to_exclude:

            if compounds_to_include == 'all':
                G.add_node(substrate, data=True, shape='o', color_map=node_color_dict['substrate'])
                G.add_edge(substrate, enzyme)
                if '↔' in reaction_equation:
                    G.add_edge(enzyme, substrate)
            else:
                if substrate in compounds_to_include:
                    G.add_node(substrate, data=True, shape='o', color_map=node_color_dict['substrate'])
                    G.add_edge(substrate, enzyme)
                    if '↔' in reaction_equation:
                        G.add_edge(enzyme, substrate)


    for product in product_list_no_num_brackets:
        if product not in compounds_to_exclude:

            if compounds_to_include == 'all':

                G.add_node(product, data=True, shape='o', color_map=node_color_dict['product'])
                G.add_edge(enzyme, product)
                if '↔' in reaction_equation:
                    G.add_edge(product, enzyme)
            else:
                if product in compounds_to_include:
                    G.add_node(product, data=True, shape='o', color_map=node_color_dict['product'])
                    G.add_edge(enzyme, product)
                    if '↔' in reaction_equation:
                        G.add_edge(product, enzyme)


def NetEnzymes(args, config_dict):

    enzymes_list_file           = args['enzymes']
    compound_list_file          = args['compounds']
    interested_ko_id            = args['ko']
    ignore_ec_with_hyphen       = args['NoHyphen']
    to_skip_file                = args['to_skip']
    plot_network                = args['plot']
    label_font_size             = args['lfs']
    node_size                   = args['ns']
    ko00001_keg                 = config_dict['ko00001_keg']
    db_file_with_ec             = config_dict['MetaCyc_rxns_with_ec']

    ########################################################################################################################

    node_color_dict = {'enzyme': 'lightgreen', 'substrate': 'grey', 'product': 'grey'}

    skip_list = set()
    if to_skip_file is not None:
        for each_to_skip in open(to_skip_file):
            skip_list.add(each_to_skip.strip())

    compounds_to_include_list = 'all'
    if compound_list_file is not None:
        compounds_to_include_list = set()
        for compound_to_include in open(compound_list_file):
            compounds_to_include_list.add(compound_to_include.strip())

    # define output file name
    ec_file_no_path, ec_file_no_ext, ec_file_ext = sep_path_basename_ext(enzymes_list_file)

    if interested_ko_id is None:
        if ignore_ec_with_hyphen is True:
            output_graphml  = '%s/%s_NoHyphen.graphml'      % (ec_file_no_path, ec_file_no_ext)
            output_plot     = '%s/%s_NoHyphen.png'          % (ec_file_no_path, ec_file_no_ext)
        else:
            output_graphml  = '%s/%s.graphml'               % (ec_file_no_path, ec_file_no_ext)
            output_plot     = '%s/%s.png'                   % (ec_file_no_path, ec_file_no_ext)
    else:
        if ignore_ec_with_hyphen is True:
            output_graphml  = '%s/%s_ko%s_NoHyphen.graphml' % (ec_file_no_path, ec_file_no_ext, interested_ko_id)
            output_plot     = '%s/%s_ko%s_NoHyphen.png'     % (ec_file_no_path, ec_file_no_ext, interested_ko_id)
        else:
            output_graphml  = '%s/%s_ko%s.graphml'          % (ec_file_no_path, ec_file_no_ext, interested_ko_id)
            output_plot     = '%s/%s_ko%s.png'              % (ec_file_no_path, ec_file_no_ext, interested_ko_id)

    ########################################################################################################################

    interested_ec_list = []
    if interested_ko_id is not None:

        print(datetime.now().strftime(time_format) + 'get ECs from interested KO')

        # read in KEGG db file
        KO_description_A_dict, KO_description_B_dict, KO_description_C_dict, KO_description_D_dict, D2ABCD_dict, ko2level_dict = get_ko2description_dict(ko00001_keg)

        # get ec list from interested KO category
        interested_ec_list = get_ec_of_interested_ko(D2ABCD_dict, KO_description_D_dict, ko2level_dict[interested_ko_id], interested_ko_id)

    # get identified_ec_list
    print(datetime.now().strftime(time_format) + 'read in provided ECs')

    identified_ec_list = set()
    for ec in open(enzymes_list_file):
        ec = ec.strip()
        if interested_ko_id is not None:
            if ec in interested_ec_list:
                if ignore_ec_with_hyphen is False:
                    identified_ec_list.add(ec)
                else:
                    if '-' not in ec:
                        identified_ec_list.add(ec)
        else:
            if ignore_ec_with_hyphen is False:
                identified_ec_list.add(ec)
            else:
                if '-' not in ec:
                    identified_ec_list.add(ec)

    # initialize a graph
    G = nx.DiGraph()

    print(datetime.now().strftime(time_format) + 'add nodes and edges to network')

    # add node and edge
    for reaction in open(db_file_with_ec):
        ec_id = reaction.strip().split('\t')[0]

        if ec_id in identified_ec_list:
            parse_biological_raction(G, reaction, skip_list, compounds_to_include_list, node_color_dict)

    print(datetime.now().strftime(time_format) + 'write out network to graphml file')

    # write out graphml
    nx.write_graphml(G, output_graphml)

    if plot_network is True:

        print(datetime.now().strftime(time_format) + 'plot network')

        # specify
        graph_layout = nx.layout.kamada_kawai_layout(G)  # kamada_kawai_layout, planar_layout, fruchterman_reingold_layout

        # turn node attributes into dict
        node_attributes_dict = {}
        for node in G.nodes(data=True):
            node_attributes_dict[node[0]] = node[1]

        print(datetime.now().strftime(time_format) + 'plot nodes')

        # plot node
        for node in G:
            nx.draw_networkx_nodes(G, graph_layout,
                                   nodelist=[node],
                                   node_size=node_size,
                                   node_color=node_attributes_dict[node]['color_map'],
                                   node_shape=node_attributes_dict[node]['shape'])

            # add customized node label
            # nx.draw_networkx_labels(g, graph_layout, nodelist=[node], font_size=8, font_color='black')

        #  all nodes label together
        nx.draw_networkx_labels(G, graph_layout, nodelist=G.nodes, font_size=label_font_size, font_color='black')

        print(datetime.now().strftime(time_format) + 'plot edges')

        # plot edges
        nx.draw_networkx_edges(G, graph_layout, width=0.5, arrows=True, arrowsize=6)

        # save plot
        plt.savefig(output_plot, dpi=300)
        plt.close()

    ########################################################################################################################

    # G_in_cytoscape_data = json_graph.cytoscape_data(G)
    # print(G_in_cytoscape_data)
    # G_in_cytoscape_graph = json_graph.cytoscape_graph(G_in_cytoscape_data)
    # print(G_in_cytoscape_data)

    print(datetime.now().strftime(time_format) + 'Done!')


if __name__ == '__main__':

    NetEnzymes_parser = argparse.ArgumentParser()

    # arguments for NetMetaCyc_parser
    NetEnzymes_parser.add_argument('-enzymes',   required=True,                           help='Enzyme list file')
    NetEnzymes_parser.add_argument('-compounds', required=False, default=None,            help='compound list file')
    NetEnzymes_parser.add_argument('-ko',        required=False, default=None,            help='get network of enzymes from specified ko')
    NetEnzymes_parser.add_argument('-to_skip',   required=False, default=None,            help='substrates/products to ignore (e.g. H2O, CO2, H+, ATP, ADP)')
    NetEnzymes_parser.add_argument('-NoHyphen',  required=False, action='store_true',     help='ignore enzymes with "-" in EC')
    NetEnzymes_parser.add_argument('-plot',      required=False, action='store_true',     help='plot network, slow and messy layout for complicated network')
    NetEnzymes_parser.add_argument('-lfs',       required=False,  default=3, type=float,  help='Font size of node labels, default is 3')
    NetEnzymes_parser.add_argument('-ns',        required=False, default=20, type=float,  help='Node size, default is 20')

    args = vars(NetEnzymes_parser.parse_args())

    NetEnzymes(args, config_dict)
