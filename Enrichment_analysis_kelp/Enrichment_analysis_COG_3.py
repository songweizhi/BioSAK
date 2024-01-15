
output_test =               '/Users/songweizhi/Desktop/Kelp_NM/enrichment_analysis/COG_enrichment_test_results_Mann_Whitney_U.tab'
COGs_detected_in_HGT =      '/Users/songweizhi/Desktop/Kelp_NM/enrichment_analysis/zKelp_recipient_genes_COG2014_wd/zKelp_recipient_genes_cog_stats_GeneNumber.txt'
COG_fun_db_file =           '/Users/songweizhi/Desktop/Kelp_NM/enrichment_analysis/cognames2003-2014.tab'


# file out
emriched_in_kelp =                              '/Users/songweizhi/Desktop/Kelp_NM/enrichment_analysis/COG_enrichment_test_emriched_in_kelp.txt'
emriched_in_tara =                              '/Users/songweizhi/Desktop/Kelp_NM/enrichment_analysis/COG_enrichment_test_emriched_in_tara.txt'
emriched_in_kelp_no_uncharacterized_protein =   '/Users/songweizhi/Desktop/Kelp_NM/enrichment_analysis/COG_enrichment_test_emriched_in_kelp_no_uncharacterized_protein.txt'
emriched_in_tara_no_uncharacterized_protein =   '/Users/songweizhi/Desktop/Kelp_NM/enrichment_analysis/COG_enrichment_test_emriched_in_tara_no_uncharacterized_protein.txt'


# get cog id to function dict
cog_id_to_function_dict = {}
cog_id_to_cate_dict = {}
for each_cog in open(COG_fun_db_file):
    each_cog_split = each_cog.strip().split('\t')
    cog_id_to_function_dict[each_cog_split[0]] = each_cog_split[2]
    cog_id_to_cate_dict[each_cog_split[0]] = each_cog_split[1]

transferred_COG_list = []
for transferred_COG in open(COGs_detected_in_HGT):
    if not transferred_COG.startswith('COG	GeneNumber'):
        transferred_COG_list.append(transferred_COG.strip().split('\t')[0])

emriched_in_sample_1_handle = open(emriched_in_kelp, 'w')
emriched_in_sample_2_handle = open(emriched_in_tara, 'w')
emriched_in_sample_1_handle.write('Category\tCOG\tP_value\tKelp\tPlanktonic\tMean_diff\tEnriched_in\tHGT\tFunction\n')
emriched_in_sample_2_handle.write('Category\tCOG\tP_value\tKelp\tPlanktonic\tMean_diff\tEnriched_in\tHGT\tFunction\n')
emriched_in_sample_1_no_uncharacterized_protein_handle = open(emriched_in_kelp_no_uncharacterized_protein, 'w')
emriched_in_sample_2_no_uncharacterized_protein_handle = open(emriched_in_tara_no_uncharacterized_protein, 'w')
emriched_in_sample_1_no_uncharacterized_protein_handle.write('Category\tCOG\tP_value\tKelp\tPlanktonic\tMean_diff\tEnriched_in\tHGT\tFunction\n')
emriched_in_sample_2_no_uncharacterized_protein_handle.write('Category\tCOG\tP_value\tKelp\tPlanktonic\tMean_diff\tEnriched_in\tHGT\tFunction\n')
n = 0
enriched_transferred_COG_list = []
for cog in open(output_test):

    if cog.startswith('COG	Kelp'):
        sample_1_name = cog.strip().split('\t')[1]
        sample_2_name = cog.strip().split('\t')[3]
    else:
        cog_split = cog.strip().split('\t')

        cog_id = cog_split[0]
        sample_1_mean = float(cog_split[1])
        sample_1_dectected_pct = float(cog_split[2])
        sample_2_mean = float(cog_split[3])
        sample_2_dectected_pct = float(cog_split[4])
        P_value_adjusted = float(cog_split[6])

        transferred = 'No'
        if cog_id in transferred_COG_list:
            transferred = 'Yes'

        if P_value_adjusted <= 0.05:

            enriched_in = ''
            if (sample_1_mean > 0) and (sample_2_mean == 0):
                if sample_1_dectected_pct >= 50:
                    enriched_in = sample_1_name
                    mean_diff = 'NA'
            elif (sample_1_mean == 0) and (sample_2_mean > 0):
                if sample_2_dectected_pct >= 50:
                    enriched_in = sample_2_name
                    mean_diff = 'NA'
            elif (sample_1_mean > 0) and (sample_2_mean > 0):
                mean_diff = float("{0:.3f}".format(sample_1_mean/sample_2_mean))
                if mean_diff >= 2:
                    enriched_in = sample_1_name
                elif mean_diff <= 0.5:
                    enriched_in = sample_2_name

            if enriched_in == sample_1_name:
                emriched_in_sample_1_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (cog_id_to_cate_dict[cog_id], cog_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, enriched_in, transferred, cog_id_to_function_dict[cog_id]))
                if 'Uncharacterized' not in cog_id_to_function_dict[cog_split[0]]:
                    emriched_in_sample_1_no_uncharacterized_protein_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (cog_id_to_cate_dict[cog_id], cog_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, enriched_in, transferred, cog_id_to_function_dict[cog_id]))
                if (transferred == 'Yes') and (cog_id not in enriched_transferred_COG_list):
                    enriched_transferred_COG_list.append(cog_id)

            elif enriched_in == sample_2_name:
                emriched_in_sample_2_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (cog_id_to_cate_dict[cog_id], cog_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, enriched_in, transferred, cog_id_to_function_dict[cog_id]))
                if 'Uncharacterized' not in cog_id_to_function_dict[cog_split[0]]:
                    emriched_in_sample_2_no_uncharacterized_protein_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (cog_id_to_cate_dict[cog_id], cog_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, enriched_in, transferred, cog_id_to_function_dict[cog_id]))

emriched_in_sample_1_handle.close()
emriched_in_sample_2_handle.close()
emriched_in_sample_1_no_uncharacterized_protein_handle.close()
emriched_in_sample_2_no_uncharacterized_protein_handle.close()


print(enriched_transferred_COG_list)
print(len(enriched_transferred_COG_list))


unenriched_transferred_COG_list = []



