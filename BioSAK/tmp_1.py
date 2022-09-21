
sam_file = '/Users/songweizhi/Desktop/assembly.fasta.sam'
read_id_txt = '/Users/songweizhi/Desktop/reads.txt'
ref_id_txt = '/Users/songweizhi/Desktop/refs.txt'

sam_file = 'assembly.fasta.sam'
read_id_txt = 'reads.txt'
ref_id_txt = 'refs.txt'

interested_read_set = set()
for each_read in open(read_id_txt):
    interested_read_set.add(each_read.strip())


ref_set = set()
for each_aln in open(sam_file):
    each_aln_split = each_aln.strip().split('\t')
    read_id = each_aln_split[0]
    ref_id = each_aln_split[2]
    if read_id in interested_read_set:
        ref_set.add(ref_id)


ref_id_txt_handle = open(ref_id_txt, 'w')
for each_ref in ref_set:
    ref_id_txt_handle.write(each_ref + '\n')
ref_id_txt_handle.close()
