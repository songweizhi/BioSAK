
sam_file     = '/Users/songweizhi/Desktop/test.sam'
reads_id_txt = '/Users/songweizhi/Desktop/test.sam.reads.id.txt'

reads_id_txt_handle = open(reads_id_txt, 'w')
for each_line in open(sam_file):
    if not each_line.startswith('@'):
        each_line_split = each_line.strip().split('\t')
        read_id = each_line_split[0]
        reads_id_txt_handle.write(read_id + '\n')
reads_id_txt_handle.close()
