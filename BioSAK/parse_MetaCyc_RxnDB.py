import os
import glob


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


MetaCyc_rxn_db_folder = '/Users/songweizhi/Desktop/MetaCyc_rxn_db'
combined_rxn_db = '/Users/songweizhi/Desktop/MetaCyc_rxn_db/MetaCyc_rxn_db.txt'


MetaCyc_rxn_db_file_re = '%s/*.dat' % MetaCyc_rxn_db_folder
MetaCyc_rxn_db_file_list = [os.path.basename(file_name) for file_name in glob.glob(MetaCyc_rxn_db_file_re)]


# get Attribute_list for all db files
Attribute_list = []
for MetaCyc_rxn_db_file in MetaCyc_rxn_db_file_list:

    pwd_MetaCyc_rxn_db_file = '%s/%s' % (MetaCyc_rxn_db_folder, MetaCyc_rxn_db_file)

    Attribute_line = False
    for rxn_attribute in open(pwd_MetaCyc_rxn_db_file):
        if rxn_attribute.startswith('#'):
            if '# Attributes:' in rxn_attribute:
                Attribute_line = True
            elif (Attribute_line is True) and (len(rxn_attribute) > 2):
                attribute_line_split = rxn_attribute.strip().split('    ')
                if attribute_line_split[1] not in Attribute_list:
                    Attribute_list.append(attribute_line_split[1])


combined_rxn_db_handle = open(combined_rxn_db, 'w')
combined_rxn_db_handle.write('Species\t%s\n' % '\t'.join(Attribute_list))
for MetaCyc_rxn_db_file in MetaCyc_rxn_db_file_list:

    pwd_MetaCyc_rxn_db_file = '%s/%s' % (MetaCyc_rxn_db_folder, MetaCyc_rxn_db_file)
    db_file_path, db_file_basename, db_file_extension = sep_path_basename_ext(pwd_MetaCyc_rxn_db_file)

    Attribute_line = False
    current_attribute_dict = {}
    for rxn in open(pwd_MetaCyc_rxn_db_file):
        if not rxn.startswith('#'):
            if rxn.strip() == '//':
                write_out_list = []
                for attribute in Attribute_list:
                    attribute_value = 'NA'
                    if attribute in current_attribute_dict:
                        attribute_value = '__'.join(current_attribute_dict[attribute])
                    write_out_list.append(attribute_value)
                combined_rxn_db_handle.write('%s\t%s\n' % (db_file_basename, '\t'.join(write_out_list)))
                current_attribute_dict = {}
            else:
                if ' - ' in rxn:
                    rxn_split = rxn.strip().split(' - ')
                    if rxn_split[0] not in current_attribute_dict:
                        current_attribute_dict[rxn_split[0]] = [rxn_split[1]]
                    else:
                        current_attribute_dict[rxn_split[0]].append(rxn_split[1])

combined_rxn_db_handle.close()




print(len(Attribute_list))
