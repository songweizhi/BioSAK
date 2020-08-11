__author__ = 'weizhisong'

#################################################################################
#                                                                               #
#     This script was writen to get the .function file for PGAP pipeline,       #
#              which get input from the WebMGA output.2 file                    #
#                                                                               #
#################################################################################

cogs = open('/Users/weizhisong/Psychrobacter/analysis_cog_annotation/output.2/PRwf-1_output.2')
genes = open('/Users/weizhisong/Psychrobacter/analysis_cog_annotation/gene_list/PRwf-1_gene_list.txt')
out = open('/Users/weizhisong/Psychrobacter/analysis_cog_annotation/gene_list/PRwf-1_function_file.txt', 'w')


# Define two dictionaries which contain the name-COG and gene-function correlation respectively.
def get_dic():
    gene_cog_dic = {}
    gene_function_dic = {}
    for cog in cogs:
        cog_split = cog.split('\t')
        gene_name = cog_split[0].strip()
        cog_id = cog_split[1].strip()
        gene_function = cog_split[10].strip()
        gene_cog_dic[gene_name] = cog_id
        gene_function_dic[gene_name] = gene_function
    return gene_cog_dic, gene_function_dic


gene_cog_dic, gene_function_dic = get_dic()

# output in the following format:
# There are three columns, the first column is the gene name, the second is COG classification and the last one is
# function description. If the gene has no COG classification information, put "-" on the corresponding position.
# These three columns are separated by <tab> and there is NO header in this table. The data are like the follows:
# 16758994 - thr operon leader peptide
# 16762656 COG4453S hypothetical protein

for gene in genes:
    gene = gene.strip()[1:]
    cog_assign = gene_cog_dic.get(gene, '-')
    function_assign = gene_function_dic.get(gene, '-')
    out.write(gene + '\t' + str(cog_assign) + '\t' + str(function_assign) + '\n')

out.close()

