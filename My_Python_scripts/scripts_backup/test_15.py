
mags_1030 = '/Users/songweizhi/Desktop/1030MAGs.txt'
gtdb_output = '/Users/songweizhi/Desktop/Kelp_NM_GTDB_r89.tsv'

mags_1030_list = []
for each in open(mags_1030):
    mags_1030_list.append(each.strip())




for each in open(gtdb_output):
    each_split = each.strip().split('\t')
    mag_id = each_split[0]
    taxon = each_split[1]

    if mag_id in mags_1030_list:

        if 'd__Archaea' in taxon:
            print(taxon)

# p__Thermoplasmatota, p__Halobacterota