

def get_pwys_with_qualified_key_enzymes(mag_to_key_enzymes_dict, pwy_to_key_enzyme_dict, key_enzyme_pct_cutoff):

    mag_to_qualified_pwys_dict = {}
    for mag in mag_to_key_enzymes_dict:
        current_mag_detected_key_enzymes = mag_to_key_enzymes_dict[mag]
        for pathway in pwy_to_key_enzyme_dict:
            current_pathway_key_enzyme = pwy_to_key_enzyme_dict[pathway]
            current_pathway_key_enzyme_qualified = []
            for key_enzyme_group in current_pathway_key_enzyme:
                key_enzyme_group_detected = set()
                for key_enzyme in key_enzyme_group:
                    if key_enzyme in current_mag_detected_key_enzymes:
                        key_enzyme_group_detected.add(key_enzyme)

                if (len(key_enzyme_group_detected)*100/len(key_enzyme_group)) >= key_enzyme_pct_cutoff:
                    current_pathway_key_enzyme_qualified.append('yes')
                else:
                    current_pathway_key_enzyme_qualified.append('no')

            if 'yes' in current_pathway_key_enzyme_qualified:
                if mag not in mag_to_qualified_pwys_dict:
                    mag_to_qualified_pwys_dict[mag] = [pathway]
                else:
                    mag_to_qualified_pwys_dict[mag].append(pathway)

    return mag_to_qualified_pwys_dict

detected_key_enzymes_dict = {'mag_4': {'TIGR04253', 'PRK2', 'PRK'}, 'mag_6': {'HBD_alignments'}, 'mag_3': {'PRK'}, 'mag_5': {'citbs', 'citas'}, 'mag_2': {'HBD_alignments'}, 'mag_1': {'TIGR04253'}}
pwy_to_key_enzyme_dict = {'|rTCA1|': [['citas', 'citbs']], '|rTCA2|': [['CCSs', 'CCSL', 'CCL']], '|4HB3HP|': [['HBD_alignments']], '|3HP|': [['TIGR04253']], '|DC4HB|': [['HBD_alignments']], '|CALVIN-PWY|': [['PRK'], ['PRK2']]}
key_enzyme_percentage_cutoff = 66
gapseq_exe = 'gapseq'


pwys_with_qualified_key_enzymes = get_pwys_with_qualified_key_enzymes(detected_key_enzymes_dict, pwy_to_key_enzyme_dict, key_enzyme_percentage_cutoff)
print(pwys_with_qualified_key_enzymes)



