
list_p = []
list_c = []
list_o = []

for each in open('/Users/songweizhi/Desktop/taxonomy_r83_March2018.tsv'):
    each_split = each.strip().split(';')
    p = each_split[1]
    c = each_split[2]
    o = each_split[3]

    if p not in list_p:
        list_p.append(p)
    if c not in list_c:
        list_c.append(c)
    if o not in list_o:
        list_o.append(o)


print(len(list_p))
print(len(list_c))
print(len(list_o))

for each_p in list_p:
    print(each_p)
