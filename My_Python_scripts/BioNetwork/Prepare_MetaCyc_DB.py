import os


db_file_raw =           '/Users/songweizhi/DB/MetaCyc/all_reactions_with_identified_enzyme.txt'
db_file_with_ec_tmp =   '/Users/songweizhi/DB/MetaCyc/all_reactions_with_ec_tmp.txt'
db_file_with_ec =       '/Users/songweizhi/DB/MetaCyc/all_reactions_with_ec.txt'

db_file_with_ec_handle = open(db_file_with_ec_tmp, 'w')
for each in open(db_file_raw):
    each_split = each.strip().split('\t')
    if len(each_split) > 1:
        if ',' in each_split[1]:
            enzyme_list = each_split[1].split(', ')
            for enzyme in enzyme_list:
                db_file_with_ec_handle.write('%s\t%s\n' % (enzyme, each_split[0]))
        else:
            db_file_with_ec_handle.write('%s\t%s\n' % (each_split[1], each_split[0]))
db_file_with_ec_handle.close()


os.system('cat %s | sort > %s' % (db_file_with_ec_tmp, db_file_with_ec))
os.system('rm %s' % db_file_with_ec_tmp)
