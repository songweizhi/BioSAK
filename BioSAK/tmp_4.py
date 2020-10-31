
def get_pwy_ec_cutoff_by_mag_cpl(mag_cpl, mag_cpl_to_cutoff_dict, mag_cpl_threshold_list):

    mag_cpl_threshold_list_high_to_low = sorted(mag_cpl_threshold_list)[::-1]

    cutoffs_determined = False
    pathway_cutoff = 0
    key_enzyme_cutoff = 0
    for cpl_threshold in mag_cpl_threshold_list_high_to_low:

        if cutoffs_determined is False:
            if mag_cpl >= cpl_threshold:
                key_enzyme_cutoff = mag_cpl_to_cutoff_dict[cpl_threshold][0]
                pathway_cutoff = mag_cpl_to_cutoff_dict[cpl_threshold][1]
                cutoffs_determined = True

    return key_enzyme_cutoff, pathway_cutoff


# cut-offs in the last two columns will be used if MAG completeness no less than the value specified in the first column.
mag_completeness = '/Users/songweizhi/Desktop/mag_completeness.txt'
cutoff_table     = '/Users/songweizhi/Desktop/cutoff_table.txt'


# read in cutoff table
mag_cpl_to_cutoff_dict = {}
mag_cpl_threshold_list = []
for line in open(cutoff_table):
    if not line.startswith('MAG	PWY	Enzyme'):
        line_split = line.strip().split()
        mag_c = int(line_split[0])
        mag_e = int(line_split[1])
        mag_p = int(line_split[2])
        mag_cpl_to_cutoff_dict[mag_c] = [mag_e, mag_p]
        mag_cpl_threshold_list.append(mag_c)


mag_completeness_dict = {}
mag_to_pathway_cutoff_dict = {}
mag_to_key_enzyme_cutoff_dict = {}
for mag in open(mag_completeness):
    mag_split = mag.strip().split('\t')
    mag_id = mag_split[0]
    mag_cpl = float(mag_split[1])
    mag_completeness_dict[mag_id] = mag_cpl
    current_mag_pathway_cutoff, current_mag_key_enzyme_cutoff = get_pwy_ec_cutoff_by_mag_cpl(mag_cpl, mag_cpl_to_cutoff_dict, mag_cpl_threshold_list)
    mag_to_pathway_cutoff_dict[mag_id] = current_mag_pathway_cutoff
    mag_to_key_enzyme_cutoff_dict[mag_id] = current_mag_key_enzyme_cutoff


print(mag_completeness_dict)
print(mag_to_key_enzyme_cutoff_dict)
print(mag_to_pathway_cutoff_dict)

