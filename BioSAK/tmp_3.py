
output_df_1                                 = '/Users/songweizhi/Desktop/1.txt'
output_df_2                                 = '/Users/songweizhi/Desktop/2.txt'

intersect_file_list                         = {'mag_1', 'mag_4', 'mag_5', 'APA_bin_1', 'mag_2', 'mag_3'}
pwys_to_detect                              = ['|3HP|', '|4HB3HP|', '|CALVIN-PWY|', '|CODHa-PWY|', '|CODHb-PWY|', '|DC4HB|', '|rTCA1|', '|rTCA2|']
pwy_to_key_enzyme_dict                      = {'|rTCA1|': [['citas', 'citbs']], '|rTCA2|': [['CCSs', 'CCSL', 'CCL']], '|4HB3HP|': [['HBD_alignments']], '|3HP|': [['TIGR04253']], '|DC4HB|': [['HBD_alignments']], '|CALVIN-PWY|': [['PRK'], ['PRK2']], '|CODHb-PWY|': [['MeTr_alignments']], '|CODHa-PWY|': [['MeTr_alignments']]}
pwy_completeness_cutoff                     = 80

detected_key_enzymes_dict                   = {'mag_1': {'TIGR04253'},
                                               'mag_4': {'PRK', 'TIGR04253', 'PRK2'},
                                               'mag_5': {'citas', 'citbs'},
                                               'APA_bin_1': set(),
                                               'mag_2': {'HBD_alignments'},
                                               'mag_3': {'PRK'}}

pwy_completeness_dict                       = {'mag_1': {'|3HP|': 80.0},
                                               'mag_4': {'|3HP|': 80.0, '|CALVIN-PWY|': 84.0},
                                               'mag_5': {'|rTCA1|': 80.0}, 'APA_bin_1': {},
                                               'mag_2': {'|4HB3HP|': 31.0, '|DC4HB|': 53.0},
                                               'mag_3': {'|CALVIN-PWY|': 92.0}}

mag_to_paths_with_qualified_key_enzyme_pct  = {'mag_1': ['|3HP|'],
                                               'mag_4': ['|3HP|', '|CALVIN-PWY|'],
                                               'mag_5': ['|rTCA1|'],
                                               'mag_2': ['|4HB3HP|', '|DC4HB|'],
                                               'mag_3': ['|CALVIN-PWY|']}

# write out as data matrix
output_df_1_handle = open(output_df_1, 'w')
output_df_2_handle = open(output_df_2, 'w')
header_printed = False
for genome in sorted(intersect_file_list):
    current_genome_header_list = []
    current_genome_value_list = []
    current_genome_header_list_2 = []
    current_genome_value_list_2 = []
    for pathway in pwys_to_detect:
        current_pathway_key_enzymes = pwy_to_key_enzyme_dict[pathway]
        group_header_top = []
        group_value_top = []
        for key_enzyme_g in current_pathway_key_enzymes:
            group_header = '_n_'.join(key_enzyme_g)
            group_value = []
            for key_enzyme in key_enzyme_g:
                if key_enzyme in detected_key_enzymes_dict[genome]:
                    group_value.append('1')
                else:
                    group_value.append('0')

            group_value_top.append('_n_'.join(group_value))
            group_header_top.append(group_header)

        current_genome_header_list.append('%s__%s' % (pathway, '_v_'.join(group_header_top)))
        current_genome_value_list.append('_v_'.join(group_value_top))

        current_genome_header_list.append('%s_completeness' % pathway)
        current_genome_header_list.append('%s_found' % pathway)
        current_genome_header_list_2.append(pathway[1:-1])

        if pathway in pwy_completeness_dict[genome]:
            current_genome_value_list.append(str(pwy_completeness_dict[genome][pathway]))
        else:
            current_genome_value_list.append('na')

        # if (pathway in mag_to_paths_with_qualified_key_enzyme_pct[genome]) and (pwy_completeness_dict[genome][pathway] >= pwy_completeness_cutoff):
        #     current_genome_value_list.append('1')
        #     current_genome_value_list_2.append('1')
        # else:
        #     current_genome_value_list.append('0')
        #     current_genome_value_list_2.append('0')

        value_to_append = '0'
        if genome in mag_to_paths_with_qualified_key_enzyme_pct:
            if pathway in mag_to_paths_with_qualified_key_enzyme_pct[genome]:
                if genome in pwy_completeness_dict:
                    if pathway in pwy_completeness_dict[genome]:
                        if pwy_completeness_dict[genome][pathway] >= pwy_completeness_cutoff:
                            value_to_append = '1'
        current_genome_value_list.append(value_to_append)
        current_genome_value_list_2.append(value_to_append)

    if header_printed is False:
        output_df_1_handle.write('Genome\t%s\n' % '\t'.join(current_genome_header_list))
        output_df_2_handle.write('Genome\t%s\n' % '\t'.join(current_genome_header_list_2))
        header_printed = True
    output_df_1_handle.write('%s\t%s\n' % (genome, '\t'.join(current_genome_value_list)))
    output_df_2_handle.write('%s\t%s\n' % (genome, '\t'.join(current_genome_value_list_2)))

output_df_1_handle.close()
output_df_2_handle.close()
