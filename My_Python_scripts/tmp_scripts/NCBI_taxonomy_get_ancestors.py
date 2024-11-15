__author__ = 'weizhisong'

############################################################################################
#                                                                                          #
#                     This script was writen to get all the parents of                     #
#       (a list of) given taxon(s), and export them in rows and columns respectively       #
#                                                                                          #
############################################################################################


nodes_file = open("/Users/weizhisong/Desktop/nodes.dmp")
bin_id_file = open('/Users/weizhisong/Desktop/Only_ID_sorted_uniq.txt')
output_file1 = open('/Users/weizhisong/Desktop/output_in_rows.txt', 'w')
output_file2 = open('/Users/weizhisong/Desktop/output_in_columns.txt', 'w')


# Generate a dic which contains all "child to parents" information in the nodes.dmp file
def get_dic():
    node_dic = {}
    for each_node in nodes_file:
        child_node = each_node.split('|')[0].rstrip('\t')
        parent_node = each_node.split('|')[1].rstrip('\t')[1:]
        node_dic[child_node] = parent_node
    return node_dic

node_dic = get_dic()


# Define a function which list all the parents of a given taxon
def get_ancestors(taxon):
    result = [taxon]
    while taxon != '1':
        result.append(node_dic.get(taxon))
        taxon = node_dic.get(taxon)
    return result


# Get all the parents of (a list of) given taxon(s), and export them in rows and columns respectively
for each_id in bin_id_file:
    each_id = each_id.rstrip('\n')
    line = get_ancestors(each_id)
    for each_node in line:
        if each_node == '1':
            output_file1.write(each_node + '\n')
        else:
            output_file1.write(each_node + '|')
        output_file2.write(each_node + '\n')

output_file1.close()
output_file2.close()

