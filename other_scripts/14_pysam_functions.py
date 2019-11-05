import pysam


def get_reads_mapped_to_a_region(sam_file, seq_id, start_pos, end_pos):

    samfile = pysam.AlignmentFile(sam_file, "rb")
    iter = samfile.fetch(seq_id, start_pos, end_pos)

    mapped_read_id_list = []
    for mapped_read in iter:
        read_id = mapped_read.query_name
        mapped_read_id_list.append(read_id)

    return mapped_read_id_list


sam_file = '/Users/songweizhi/Desktop/pysam_test/ref_seq.bam'
mapped_read_list = get_reads_mapped_to_a_region(sam_file, 'ref_seq', 43, 44)
print(mapped_read_list)
