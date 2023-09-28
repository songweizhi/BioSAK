

file_in = 'input_file.txt'

col_index = {}
for each_line in open(file_in):
    each_line_split = each_line.strip().split('\t')
    if each_line.startswith('\t'):
        col_index = {key: i for i, key in enumerate(each_line_split)}
    else:
        gene_1    = each_line_split[col_index['gene_1']]
        gene_2    = each_line_split[col_index['gene_2']]
        direction = each_line_split[col_index['direction']]

