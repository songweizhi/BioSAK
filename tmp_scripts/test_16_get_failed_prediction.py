
wd = '/Users/songweizhi/Desktop/new'
transfer_file = open('%s/donor2recip.txt' % wd)
core_genes = open('/Users/songweizhi/Desktop/get_core/qualified_core_genes.txt')


donor_list = []
for each in transfer_file:
    each_split = each.strip().split('\t')
    donor = each_split[0]
    if donor != 'donor_gene':
        donor_list.append(donor)


for each in core_genes:
    each_split = each.strip().split('\t')
    exist = 0
    for each1 in donor_list:
        if each1 in each_split:
            exist = 1
    if exist == 0:
        print(each.strip())