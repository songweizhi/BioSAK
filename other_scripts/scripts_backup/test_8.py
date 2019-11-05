
file = open('/Users/songweizhi/Desktop/googleVis_MetaBAT_MyCC_vs_CONCOCT_size_cutoff_0.5MB.csv')
handle = open('/Users/songweizhi/Desktop/googlevis_input.csv','w')
handle.write('A,B,C\n')
for each in file:
    print(each.strip())
    each_split = each.strip().split(',')
    meta_mycc = each_split[0]
    concoct = each_split[1]
    length = each_split[2]
    meta_mycc_split = meta_mycc.split('Cluster')
    metabat = meta_mycc_split[0]
    mycc = 'Cluster' + meta_mycc_split[1]

    metabat_binname = 'MetaBAT_' + metabat.split('.')[1]
    mycc_binname = 'MyCC_' + mycc.split('.')[-1]
    concoct_binname = 'CONCOCT_' + concoct.split('_')[1]
    print(metabat_binname)
    print(mycc_binname)
    print(concoct_binname)
    handle.write('%s,%s,%s\n' % (metabat_binname, mycc_binname, length))
    handle.write('%s,%s,%s\n' % (mycc_binname, concoct_binname, length))



    print('\n')

handle.close()
