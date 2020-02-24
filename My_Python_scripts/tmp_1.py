
#
#


drep99_bin_list = []
for each in open('/Users/songweizhi/Desktop/bins.txt'):

    drep99_bin_list.append(each.strip())

print(drep99_bin_list)
print(len(drep99_bin_list))

n = 0
for each in open('/Users/songweizhi/Desktop/Tara_NM_CheckM_0.4_0.05.csv'):
    each_split = each.strip().split(',')



    if each_split[0] in drep99_bin_list:



        print(each_split)

        n += 1

print(n)