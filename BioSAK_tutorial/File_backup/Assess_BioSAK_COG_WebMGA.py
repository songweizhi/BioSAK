
MAG_id = '18'  # 11, 13, 18
WebMGA_COG =         '/Users/songweizhi/Desktop/Assess_BioTAK/COG/BH_ER_050417_Refined_%s_WebMGA.txt'                % MAG_id
BioSAK_Diamond_COG = '/Users/songweizhi/Desktop/Assess_BioTAK/COG/BH_ER_050417_Refined_%s_query_to_cog_diamond.txt'  % MAG_id


query_gene_list = set()
BioSAK_Diamond_COG_dict = {}
for ko in open(BioSAK_Diamond_COG):
    ko_split = ko.strip().split('\t')
    query_gene_list.add(ko_split[0])
    if len(ko_split) > 1:
        if ko_split[0] not in BioSAK_Diamond_COG_dict:
            BioSAK_Diamond_COG_dict[ko_split[0]] = [ko_split[1]]
        else:
            BioSAK_Diamond_COG_dict[ko_split[0]].append(ko_split[1])


WebMGA_COG_dict = {}
for ko in open(WebMGA_COG):
    ko_split = ko.strip().split('\t')
    query_gene_list.add(ko_split[0])
    if len(ko_split) > 1:
        if ko_split[0] not in WebMGA_COG_dict:
            WebMGA_COG_dict[ko_split[0]] = [ko_split[1]]
        else:
            WebMGA_COG_dict[ko_split[0]].append(ko_split[1])


total = 0
both_NA_num = 0
same_cog_num = 0
only_BioSAK_Diamond_num = 0
only_WebMGA_cog = 0
different_cog_num = 0
for query_gene in query_gene_list:

    BioSAK_Diamond_cog = 'NA'
    if query_gene in BioSAK_Diamond_COG_dict:
        BioSAK_Diamond_cog = BioSAK_Diamond_COG_dict[query_gene]

    WebMGA_cog = 'NA'
    if query_gene in WebMGA_COG_dict:
        WebMGA_cog = WebMGA_COG_dict[query_gene]

    if (BioSAK_Diamond_cog == WebMGA_cog) and (BioSAK_Diamond_cog == 'NA'):
        both_NA_num += 1

    elif (BioSAK_Diamond_cog == WebMGA_cog) and (BioSAK_Diamond_cog != 'NA'):
        same_cog_num += 1

    elif (BioSAK_Diamond_cog != 'NA') and (WebMGA_cog == 'NA'):
        only_BioSAK_Diamond_num += 1

    elif (BioSAK_Diamond_cog == 'NA') and (WebMGA_cog != 'NA'):
        only_WebMGA_cog += 1

    elif (BioSAK_Diamond_cog != 'NA') and (WebMGA_cog != 'NA') and (BioSAK_Diamond_cog != WebMGA_cog):

        if (len(BioSAK_Diamond_cog) == 1) and (len(WebMGA_cog) == 1):
            different_cog_num += 1
        else:
            if BioSAK_Diamond_cog[0] in WebMGA_cog:
                same_cog_num += 1
            else:
                different_cog_num += 1

    total += 1


print('Query gene\t%s' % total)
print('Both NA\t%s' % both_NA_num)
print('Same COG\t%s' % same_cog_num)
print('Different COG\t%s' % different_cog_num)
print('WebMGA only\t%s' % only_WebMGA_cog)
print('BioSAK (Diamond) only\t%s' % only_BioSAK_Diamond_num)
