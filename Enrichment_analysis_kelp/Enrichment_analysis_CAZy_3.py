
output_test =                   '/Users/songweizhi/Desktop/enrichment_analysis/CAZy_enrichment_analysis_Mann_Whitney_U.tab'
CAZy_family_activities_file =   '/Users/songweizhi/Desktop/enrichment_analysis/CAZyDB.fam-activities.txt'
CAZys_detected_in_HGT =         '/Users/songweizhi/Desktop/enrichment_analysis/zKelp_recipient_genes_dbCAN_wd/zKelp_recipient_genes_dbCAN_stats_GeneNumber.txt'

# file out
emriched_in_kelp =              '/Users/songweizhi/Desktop/enrichment_analysis/CAZy_enrichment_test_emriched_in_kelp.txt'
emriched_in_tara =              '/Users/songweizhi/Desktop/enrichment_analysis/CAZy_enrichment_test_emriched_in_tara.txt'
emriched_in_na =                '/Users/songweizhi/Desktop/enrichment_analysis/CAZy_enrichment_test_emriched_in_na.txt'

# get cog id to function dict
cazy_id_to_function_dict = {}
for each_cazy in open(CAZy_family_activities_file):
    if not each_cazy.startswith('#'):
        each_cazy_split = each_cazy.strip().split('	  ')

        if len(each_cazy_split) == 1:
           cazy_id_to_function_dict[each_cazy_split[0]] = 'NA'
        elif len(each_cazy_split) > 1:
           cazy_id_to_function_dict[each_cazy_split[0]] = each_cazy_split[1]

transferred_CAZy_list = []
for transferred_CAZy in open(CAZys_detected_in_HGT):
    if not transferred_CAZy.startswith('Family'):
        transferred_CAZy_list.append(transferred_CAZy.strip().split('\t')[0])

emriched_in_sample_1_handle = open(emriched_in_kelp, 'w')
emriched_in_sample_2_handle = open(emriched_in_tara, 'w')
emriched_in_sample_1_handle.write('CAZy\tP_value\tKelp\tPlanktonic\tMean_diff\tEnriched_in\tHGT\tFunction\n')
emriched_in_sample_2_handle.write('CAZy\tP_value\tKelp\tPlanktonic\tMean_diff\tEnriched_in\tHGT\tFunction\n')
n = 0
for cazy in open(output_test):

    if cazy.startswith('CAZy	Kelp'):
        sample_1_name = cazy.strip().split('\t')[1]
        sample_2_name = cazy.strip().split('\t')[3]

    else:
        cazy_split = cazy.strip().split('\t')

        cazy_id = cazy_split[0]
        sample_1_mean = float(cazy_split[1])
        sample_1_dectected_pct = float(cazy_split[2])
        sample_2_mean = float(cazy_split[3])
        sample_2_dectected_pct = float(cazy_split[4])
        P_value_adjusted = float(cazy_split[6])

        transferred = 'No'
        if cazy_id in transferred_CAZy_list:
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
                emriched_in_sample_1_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (cazy_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, enriched_in, transferred, cazy_id_to_function_dict[cazy_id]))
            elif enriched_in == sample_2_name:
                emriched_in_sample_2_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (cazy_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, enriched_in, transferred, cazy_id_to_function_dict[cazy_id]))

emriched_in_sample_1_handle.close()
emriched_in_sample_2_handle.close()
