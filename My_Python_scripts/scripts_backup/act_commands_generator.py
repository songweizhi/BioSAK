matches = open("/Users/weizhisong/Desktop/28matches.txt")
out = open("/Users/weizhisong/Desktop/28matches_commands.txt", 'w')
for match in matches:
    match = match.strip()
    match_split = match.split('___')
    gene_1 = './' + match + '/' + match_split[0] + '.gbk '
    gene_2 = './' + match + '/' + match_split[1] + '.gbk '
    comparison_file = './' + match + '/' + match + '.txt '
    act = '/Users/weizhisong/Softwares/artemis/act '
    command = act + gene_1 + comparison_file + gene_2 + '\n'
    out.write(command)
    print(command)
out.close()
