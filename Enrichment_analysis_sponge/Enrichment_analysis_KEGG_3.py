

def summarize_stats(output_test, summary_txt):

    #output_test = '/Users/songweizhi/PycharmProjects/BioSAK/demo_data/enrich/MAGs_dRep99_KEGG_D_GeneNumber_Mann_Whitney_U.txt'
    #summary_txt = '/Users/songweizhi/PycharmProjects/BioSAK/demo_data/enrich/summary.txt'

    summary_txt_handle = open(summary_txt, 'w')
    summary_txt_handle.write('KO\tP_value\tsponge\tseawater\tMean_diff\tEnriched\n')
    line_num_index = 0
    for ko in open(output_test):
        if line_num_index == 0:
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

                if (enriched_in == sample_1_name) or (enriched_in == sample_2_name):
                    summary_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, enriched_in))
        line_num_index += 1
    summary_txt_handle.close()


