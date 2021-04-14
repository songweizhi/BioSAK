
ko_assign = '/Users/songweizhi/Desktop/user_ko.txt'

op_folder = '/Users/songweizhi/Desktop/ko_sep'

for each_gene in open(ko_assign):
    print(each_gene.strip())
    genome_id = '_'.join(each_gene.strip().split('\t')[0].split('_')[:-1])
    print(genome_id)
    pwd_ko_file = '%s/%s_ko.txt' % (op_folder, genome_id)

    pwd_ko_file_handle = open(pwd_ko_file, 'a')
    pwd_ko_file_handle.write(each_gene)
    pwd_ko_file_handle.close()

