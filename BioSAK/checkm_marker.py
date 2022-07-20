from Bio import SeqIO


checkm_qa_6_txt = '/Users/songweizhi/Desktop/op_6.txt'
checkm_qa_9_fa  = '/Users/songweizhi/Desktop/op_9.fa'


mag_to_duplicated_marker_dict = {}
marker_name_to_identified_marker_dict = {}
for each_line in open(checkm_qa_6_txt):
    if not each_line.startswith('Bin Id	Marker Id	Gene Ids'):
        each_line_split = each_line.strip().split('\t')
        mag_id = each_line_split[0]
        marker_name = each_line_split[1]
        marker_name_with_mag_id = '%s___%s' % (mag_id, marker_name)
        identified_marker_list = each_line_split[2].split(',')

        if mag_id not in mag_to_duplicated_marker_dict:
            mag_to_duplicated_marker_dict[mag_id] = {marker_name}
        else:
            mag_to_duplicated_marker_dict[mag_id].add(marker_name)

        marker_name_to_identified_marker_dict[marker_name_with_mag_id] = identified_marker_list


print(mag_to_duplicated_marker_dict)
print(marker_name_to_identified_marker_dict)


for each_seq in SeqIO.parse(checkm_qa_9_fa, 'fasta'):
    seq_des = each_seq.description
    seq_des_split = seq_des.split(' ')
    print(seq_des_split)

    mag_id = seq_des_split[0]
    ctg_id = seq_des_split[1]






