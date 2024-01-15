
# file in
ko_level = 'D'  # B, C, D
output_test =               '/Users/songweizhi/Desktop/enrichment_analysis/KEGG_%s_enrichment_analysis_Mann_Whitney_U.tab'  % ko_level
KOs_detected_in_HGT =       '/Users/songweizhi/Desktop/enrichment_analysis/zKelp_recipient_genes_KEGG_wd/zKelp_recipient_genes_ko_stats_D_GeneNumber.txt'
KEGG_DB_ko =                '/Users/songweizhi/Desktop/enrichment_analysis/ko00001.keg'

# file out
emriched_in_kelp =                              '/Users/songweizhi/Desktop/enrichment_analysis/KEGG_%s_enrichment_test_emriched_in_kelp.txt'                                % ko_level
emriched_in_tara =                              '/Users/songweizhi/Desktop/enrichment_analysis/KEGG_%s_enrichment_test_emriched_in_tara.txt'                                % ko_level
emriched_in_kelp_no_uncharacterized_protein =   '/Users/songweizhi/Desktop/enrichment_analysis/KEGG_%s_enrichment_test_emriched_in_kelp_no_uncharacterized_protein.txt'     % ko_level
emriched_in_tara_no_uncharacterized_protein =   '/Users/songweizhi/Desktop/enrichment_analysis/KEGG_%s_enrichment_test_emriched_in_tara_no_uncharacterized_protein.txt'     % ko_level


transferred_KO_list = []
for transferred_KO in open(KOs_detected_in_HGT):
    if not transferred_KO.startswith('KO	GeneNumber'):
        transferred_KO_list.append(transferred_KO.strip().split('\t')[0])

# store ko functions in dict
As_description_dict = {}
Bs_description_dict = {}
Cs_description_dict = {}
Ds_description_dict = {}
D2ABCD_dict = {}
current_A = ''
current_B = ''
current_C = ''
for each_line in open(KEGG_DB_ko):
    if each_line[0] in ['A', 'B', 'C', 'D']:
        each_line_split = each_line.strip().split(' ')

        if each_line[0] == 'A':
            current_A_id = each_line_split[0]
            current_A_description = ' '.join(each_line_split[1:])
            current_A = current_A_id
            As_description_dict[current_A_id] = current_A_description

        elif each_line[0] == 'B':
            if len(each_line_split) > 1:
                current_B_id = each_line_split[2]
                current_B_description = ' '.join(each_line_split[3:])
                current_B = current_B_id
                Bs_description_dict[current_B_id] = current_B_description

        elif each_line[0] == 'C':
            current_C_id = each_line_split[4]
            current_C_description = ' '.join(each_line_split[5:])
            current_C = current_C_id
            Cs_description_dict[current_C_id] = current_C_description

        elif each_line[0] == 'D':
            current_D_id = each_line_split[6]
            current_D_description = ' '.join(each_line_split[7:])
            Ds_description_dict[current_D_id] = current_D_description
            ABCD_value = 'A_%s|B_%s|C_%s|D_%s' % (current_A, current_B, current_C, current_D_id)
            if current_D_id not in D2ABCD_dict:
                D2ABCD_dict[current_D_id] = [ABCD_value]
            elif (current_D_id in D2ABCD_dict) and (ABCD_value not in D2ABCD_dict[current_D_id]):
                D2ABCD_dict[current_D_id].append(ABCD_value)


emriched_in_sample_1_handle = open(emriched_in_kelp, 'w')
emriched_in_sample_2_handle = open(emriched_in_tara, 'w')
emriched_in_sample_1_handle.write('KO\tP_value\tKelp\tPlanktonic\tMean_diff\tEnriched_in\tHGT\tFunction\n')
emriched_in_sample_2_handle.write('KO\tP_value\tKelp\tPlanktonic\tMean_diff\tEnriched_in\tHGT\tFunction\n')
emriched_in_sample_1_no_uncharacterized_protein_handle = open(emriched_in_kelp_no_uncharacterized_protein, 'w')
emriched_in_sample_2_no_uncharacterized_protein_handle = open(emriched_in_tara_no_uncharacterized_protein, 'w')
emriched_in_sample_1_no_uncharacterized_protein_handle.write('KO\tP_value\tKelp\tPlanktonic\tMean_diff\tEnriched_in\tHGT\tFunction\n')
emriched_in_sample_2_no_uncharacterized_protein_handle.write('KO\tP_value\tKelp\tPlanktonic\tMean_diff\tEnriched_in\tHGT\tFunction\n')
n = 0
for ko in open(output_test):

    if ko.startswith('KO	Kelp'):
        sample_1_name = ko.strip().split('\t')[1]
        sample_2_name = ko.strip().split('\t')[3]

    else:
        ko_split = ko.strip().split('\t')

        ko_id = ko_split[0]
        sample_1_mean = float(ko_split[1])
        sample_1_dectected_pct = float(ko_split[2])
        sample_2_mean = float(ko_split[3])
        sample_2_dectected_pct = float(ko_split[4])
        P_value_adjusted = float(ko_split[6])

        transferred = 'No'
        if ko_id in transferred_KO_list:
            transferred = 'Yes'

        ko_function = ''
        if ko_level == 'B':
            ko_function = Bs_description_dict[ko_id]
        if ko_level == 'C':
            ko_function = Cs_description_dict[ko_id]
        if ko_level == 'D':
            ko_function = Ds_description_dict[ko_id]

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
                emriched_in_sample_1_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, enriched_in, transferred, ko_function))
                if 'uncharacterized' not in ko_function:
                    emriched_in_sample_1_no_uncharacterized_protein_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, enriched_in, transferred, ko_function))

            elif enriched_in == sample_2_name:
                emriched_in_sample_2_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, enriched_in, transferred, ko_function))
                if 'uncharacterized' not in ko_function:
                    emriched_in_sample_2_no_uncharacterized_protein_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, enriched_in, transferred, ko_function))

emriched_in_sample_1_handle.close()
emriched_in_sample_2_handle.close()
