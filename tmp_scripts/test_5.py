
NonRhodophyta_with_hitsTOnr = '/Users/songweizhi/Desktop/000/NonRhodophyta_with_hitsTOnr.fasta'
Rhodophyta_Specific_Final_transcriptome = '/Users/songweizhi/Desktop/000/Rhodophyta_Specific_Final_transcriptome.fasta'
id_list_file = '/Users/songweizhi/Desktop/000/annotated.txt'

seq_id_set = set()
for each_line in open(NonRhodophyta_with_hitsTOnr):
    if each_line.startswith('>'):
        seq_id = each_line[1:].split('|')[0]
        seq_id_set.add(seq_id)
print(len(seq_id_set))

for each_line in open(Rhodophyta_Specific_Final_transcriptome):
    if each_line.startswith('>'):
        seq_id = each_line[1:].split('|')[0]
        seq_id_set.add(seq_id)
print(len(seq_id_set))

id_list_file_handle = open(id_list_file, 'w')
for each_seq in seq_id_set:
    id_list_file_handle.write(each_seq + '\n')
id_list_file_handle.close()
