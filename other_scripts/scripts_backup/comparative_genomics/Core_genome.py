__author__ = 'weizhisong'

blast_result = open('/Users/weizhisong/Desktop/blast_result_273-4.txt')

def get_gene_list():
    gene_list = []
    for each_line in blast_result:
        query_name = each_line.split("\t")[0]
        if query_name in gene_list:
            pass
        else:
            gene_list.append(query_name)
    return gene_list

gene_list = get_gene_list()

for gene in gene_list:
    print gene
    for each_line in blast_result:
        query_name = each_line.split("\t")[0].strip('\n')
        subject_name = each_line.split("\t")[1].strip('\n')
        if query_name == gene:
            print "a"




 #   gene_cluster = []
  #  gene_cluster.append(gene)
   # b = gene_cluster
#    print gene_cluster
#    for each_line in blast_result:
#        print each_line

        #query_name = each_line.split("\t")[0].strip('\n')
        #subject_name = each_line.split("\t")[1].strip('\n')
        #print subject_name









#    subject_name = each_line.split("\t")[1]
#    e_value = each_line.split("\t")[10]
#    gene_list.append(query_name)
#    print query_name + "\t" + "\t" + subject_name






#    if len(gene_clusters) == 0:
#       gene_clusters.append(query_name)
#
#   else:
#      gene_clusters.append(subject_name)
# print gene_clusters


