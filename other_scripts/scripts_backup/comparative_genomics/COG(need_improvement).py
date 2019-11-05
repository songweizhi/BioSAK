__author__ = 'weizhisong'

cog_db = open('/Users/weizhisong/Desktop/whog')
matches = open('/Users/weizhisong/Desktop/matches.tab')
out = open ('/Users/weizhisong/Desktop/matches_cog.tab','w')

cog_catergory_dic = {}
cog_group_dic = {}

for each_line in cog_db:
    each_line = each_line.strip()
    if each_line.startswith('['):
        cog_value = each_line.split(' ')[1]
        cog_group = each_line.split(' ')[0][1]
        cog_group_dic[cog_value] = cog_group
    elif ':' in each_line:
        key = each_line.split(':')[1][2:]
        if ' ' in key:
            sub_keys = key.split(' ')
            for sub_key in sub_keys:
                cog_catergory_dic[sub_key] = cog_value
        else:
            cog_catergory_dic[key] = cog_value
    elif '_______' not in each_line and each_line != '':
        each_line_more = each_line.split(' ')
        for sub_key_more in each_line_more:
            cog_catergory_dic[sub_key_more] = cog_value

for match in matches:
    match = match.strip().split('\t')
    query = match[0]
    target = match[1]
    if target not in cog_catergory_dic.keys():
        pass
    else:
        out_print = query + '\t' + target + '\t' + cog_catergory_dic[target] + '\t' + cog_group_dic[cog_catergory_dic[target]] + '\n'
        out.write(out_print)
        print(out_print)
out.close()
  #      print query + '\t' + target + '\t' + cog_catergory_dic[target]

#print cog_catergory_dic['yi22_g1']