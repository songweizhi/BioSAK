
import os

os.chdir('/Users/songweizhi/Desktop/test')

dataset = 'ns'

grouping_file = '%s_species_tree_grouping_modified.txt' % dataset
gtdbtk_op = '%s_gtdbtk.tsv' % dataset

out = '%s_out.txt' % dataset
out_sorted = '%s_out_sorted.txt' % dataset

bin_to_group_dict = {}
for each in open(grouping_file):
    each_split = each.strip().split(',')
    bin_to_group_dict[each_split[1]] = each_split[0]
print(bin_to_group_dict)


out_handle = open(out, 'w')
for each2 in open(gtdbtk_op):
    each2_split = each2.strip().split('\t')
    bin = each2_split[0]
    assign = each2_split[1]
    #print('%s\t%s\t%s' % (bin_to_group_dict[bin], bin, assign))
    out_handle.write('%s\t%s\t%s\n' % (bin_to_group_dict[bin], bin, assign))

out_handle.close()



os.system('cat %s | sort > %s' % (out, out_sorted))











