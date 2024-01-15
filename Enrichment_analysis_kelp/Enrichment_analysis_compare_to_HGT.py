
# file in
COGs_enriched_in_MAGs =     '/Users/songweizhi/Desktop/enrichment_analysis/COG_enrichment_test_emriched_in_kelp_no_uncharacterized_protein.txt'
COGs_detected_in_HGT =      '/Users/songweizhi/Desktop/enrichment_analysis/zKelp_recipient_genes_COG2014_wd/zKelp_recipient_genes_cog_stats_GeneNumber.txt'

KOs_enriched_in_MAGs =      '/Users/songweizhi/Desktop/enrichment_analysis/KEGG_D_enrichment_test_emriched_in_kelp_no_uncharacterized_protein.txt'
KOs_detected_in_HGT =       '/Users/songweizhi/Desktop/enrichment_analysis/zKelp_recipient_genes_KEGG_wd/zKelp_recipient_genes_ko_stats_D_GeneNumber.txt'

CAZys_enriched_in_MAGs =    '/Users/songweizhi/Desktop/enrichment_analysis/CAZy_enrichment_test_emriched_in_kelp.txt'
CAZys_detected_in_HGT =     '/Users/songweizhi/Desktop/enrichment_analysis/zKelp_recipient_genes_dbCAN_wd/zKelp_recipient_genes_dbCAN_stats_GeneNumber.txt'


enriched_COG_list = []
for enriched_COG in open(COGs_enriched_in_MAGs):
    if not enriched_COG.startswith('COG'):
        enriched_COG_list.append(enriched_COG.strip().split('\t')[1])

enriched_KO_list = []
for enriched_KO in open(KOs_enriched_in_MAGs):
    if not enriched_KO.startswith('KO'):
        enriched_KO_list.append(enriched_KO.strip().split('\t')[0])

enriched_CAZy_list = []
for enriched_CAZy in open(CAZys_enriched_in_MAGs):
    if not enriched_CAZy.startswith('CAZy'):
        enriched_CAZy_list.append(enriched_CAZy.strip().split('\t')[0])

transferred_COG_list = []
for transferred_COG in open(COGs_detected_in_HGT):
    if not transferred_COG.startswith('COG	GeneNumber'):
        transferred_COG_list.append(transferred_COG.strip().split('\t')[0])

transferred_KO_list = []
for transferred_KO in open(KOs_detected_in_HGT):
    if not transferred_KO.startswith('KO	GeneNumber'):
        transferred_KO_list.append(transferred_KO.strip().split('\t')[0])

transferred_CAZy_list = []
for transferred_CAZy in open(CAZys_detected_in_HGT):
    if not transferred_CAZy.startswith('Family'):
        transferred_CAZy_list.append(transferred_CAZy.strip().split('\t')[0])


print(enriched_COG_list)
print(enriched_KO_list)
print(enriched_CAZy_list)
print(transferred_COG_list)
print(transferred_KO_list)
print(transferred_CAZy_list)




for enriched_COG in enriched_COG_list:

    transferred = 'No'
    if enriched_COG in transferred_COG_list:
        transferred = 'Yes'

    print('%s\t%s' % (enriched_COG, transferred))








