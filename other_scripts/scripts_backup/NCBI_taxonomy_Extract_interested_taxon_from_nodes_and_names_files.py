__author__ = 'weizhisong'

#################################################################################
#                                                                               #
#             This script was writen to extract the information of              #
#          (a list of) interested taxon(s) from NCBI Taxonomy database          #
#                                                                               #
#################################################################################


file1 = open('/Users/weizhisong/Desktop/all_nodes_sorted_uniq.txt')
nodes_file = open('/Users/weizhisong/Research/NCBI_Taxonomy/taxdump/nodes.dmp')
names_file = open('/Users/weizhisong/Research/NCBI_Taxonomy/taxdump/names.dmp')
nodes_output = open('/Users/weizhisong/Desktop/customized_nodes.dmp', 'w')
names_output = open('/Users/weizhisong/Desktop/customized_names.dmp', 'w')


# Generate a list which contains all interested taxons
id_list = []

for file1_each_id in file1:
    file1_each_id = file1_each_id.rstrip('\n')
    id_list.append(file1_each_id)


# Extract information from name.dmp and export to a file named customized_names.dmp
for names_each_id in names_file:
    name_type = names_each_id.split('|')[3].rstrip('\t')[1:]
    names_each_id_id = names_each_id.split('|')[0].rstrip('\t')
    if name_type == 'scientific name':
        if names_each_id_id in id_list:
            names_output.write(names_each_id)
        else:
            pass
    else:
        pass

names_output.close()


# Extract information from nodes.dmp and export to a file named customized_nodes.dmp
for nodes_each_id in nodes_file:
    nodes_each_id_id = nodes_each_id.split('|')[0].rstrip('\t')
    if nodes_each_id_id in id_list:
        nodes_output.write(nodes_each_id)
    else:
        pass

nodes_output.close()
