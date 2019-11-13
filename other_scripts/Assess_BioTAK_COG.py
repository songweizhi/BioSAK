

MAG_id = '13'  # 11, 13, 18
GhostKOALA_ko =      '/Users/songweizhi/Desktop/Assess_BioTAK/COG/BH_ER_050417_Refined_%s_query_to_cog_blast.txt' % MAG_id
local_ko= '/Users/songweizhi/Desktop/Assess_BioTAK/COG/BH_ER_050417_Refined_%s_query_to_cog_diamond.txt' % MAG_id


query_gene_list = set()
local_ko_dict = {}
for ko in open(local_ko):
    ko_split = ko.strip().split('\t')
    query_gene_list.add(ko_split[0])
    if len(ko_split) > 1:
        local_ko_dict[ko_split[0]] = ko_split[1]

GhostKOALA_ko_dict = {}
for ko in open(GhostKOALA_ko):
    ko_split = ko.strip().split('\t')
    query_gene_list.add(ko_split[0])
    if len(ko_split) > 1:
        GhostKOALA_ko_dict[ko_split[0]] = ko_split[1]


total = 0
both_NA_num = 0
same_ko_num = 0
only_local_num = 0
only_ghost_num = 0
different_ko_num = 0
for query_gene in query_gene_list:

    query_gene_ko_local = 'NA'
    if query_gene in local_ko_dict:
        query_gene_ko_local = local_ko_dict[query_gene]

    query_gene_ko_ghost = 'NA'
    if query_gene in GhostKOALA_ko_dict:
        query_gene_ko_ghost = GhostKOALA_ko_dict[query_gene]

    if (query_gene_ko_local == query_gene_ko_ghost) and (query_gene_ko_local == 'NA'):
        both_NA_num += 1

    elif (query_gene_ko_local == query_gene_ko_ghost) and (query_gene_ko_local != 'NA'):
        same_ko_num += 1

    elif (query_gene_ko_local != 'NA') and (query_gene_ko_ghost == 'NA'):
        only_local_num += 1

    elif (query_gene_ko_local == 'NA') and (query_gene_ko_ghost != 'NA'):
        only_ghost_num += 1

    elif (query_gene_ko_local != 'NA') and (query_gene_ko_ghost != 'NA') and (query_gene_ko_local != query_gene_ko_ghost):
        different_ko_num += 1

    total += 1


print('query gene: %s' % total)
print('both_NA_num: %s' % both_NA_num)
print('same_ko_num: %s' % same_ko_num)
print('different_ko_num: %s' % different_ko_num)
print('only_local_num: %s' % only_local_num)
print('only_ghost_num: %s' % only_ghost_num)
