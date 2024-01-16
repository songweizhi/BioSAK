
file_in = 'input_file.txt'

col_index = dict()
line_num_index = 0
for each_line in open(file_in):
    line_num_index += 1
    line_split = each_line.strip().split('\t')
    if line_num_index == 1:
        col_index = {key: i for i, key in enumerate(line_split)}
    else:
        gene_1    = line_split[col_index['gene_1']]
        gene_2    = line_split[col_index['gene_2']]
        direction = line_split[col_index['direction']]
