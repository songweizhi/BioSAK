

SRA_Biosample_dict = dict()
SRA_Bioproject_dict = dict()
SRA_Biosample_Bioproject_txt = '/Users/songweizhi/Desktop/SRA_Biosample_Bioproject.txt'
for each_line in open(SRA_Biosample_Bioproject_txt):
    each_split = each_line.strip().split('\t')
    SRA_Biosample_dict[each_split[0]] = each_split[1]
    SRA_Bioproject_dict[each_split[0]] = each_split[2]

metadata_dict = dict()
metadata_txt = '/Users/songweizhi/Desktop/biosamples_uniq/metadata.txt'
for each_line in open(metadata_txt):
    each_split = each_line.strip().split('\t')
    metadata_dict[each_split[0]] = each_line.strip()




op_txt = '/Users/songweizhi/Desktop/op.txt'
op_txt_handle = open(op_txt, 'w')

for each in open("/Users/songweizhi/Desktop/150_MAGs.txt"):
    each_split = each.strip().split("\t")
    op_str = ''
    if ('SRR' in each_split[1]) or ('ERR' in each_split[1]):
        current_biosample = SRA_Biosample_dict.get(each_split[1], 'na')
        current_bioproject = SRA_Bioproject_dict.get(each_split[1], 'na')
        op_str = '%s\t%s\t%s\t%s\t%s' % (each_split[0], each_split[1], current_biosample, current_bioproject, metadata_dict.get(current_biosample, 'na'))
    else:
        op_str = '%s\t%s\t%s\t%s\t%s' % (each_split[0], 'na', each_split[1], 'na', metadata_dict.get(each_split[1], 'na'))
    op_txt_handle.write(op_str + '\n')
    print(op_str)
op_txt_handle.close()