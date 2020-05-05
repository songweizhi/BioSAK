
file_in  = '/Users/songweizhi/Desktop/combined_R1.fasta'
file_out = '/Users/songweizhi/Desktop/combined_R1_formatted.fasta'
read_strand = '1'

file_out_handle = open(file_out, 'w')
for each in open(file_in):

    if each.startswith('>'):
        each_split = each.strip().split('_')
        file_out_handle.write('>read.%s/%s\n' % (each_split[0][2:], read_strand))
    else:
        file_out_handle.write(each)

file_out_handle.close()

