import os
from Bio import SeqIO

################################### CONFIGURATION ###################################

length_setting = 30000
contig_list_file = '/Users/weizhisong/Desktop/contig_list_2.txt'
assemble_file = '/Users/weizhisong/Desktop/combined_after.fas'
output_folder = '/Users/weizhisong/Desktop/circle3/'
act_executable_file = '/Users/weizhisong/Softwares/artemis/act'

#####################################################################################
os.system('mkdir ' + output_folder)
contigs = open(assemble_file)
contig_names = open(contig_list_file)

output_act = open(output_folder + 'act_commands.txt', 'w')

contig_list = []
for each_line in contig_names:
    each_line = each_line.strip()
    contig_list.append(each_line)

os.system('mkdir ' + output_folder + 'blast_results_outfmt_0_' + str(length_setting))
os.system('mkdir ' + output_folder + 'blast_results_outfmt_6_' + str(length_setting))
blast_folder_outfmt0 = 'blast_results_outfmt_0_' + str(length_setting)
blast_folder_outfmt6 = 'blast_results_outfmt_6_' + str(length_setting)
for seq_records in SeqIO.parse(contigs, 'fasta'):
    if seq_records.id in contig_list:
        sequence_length = len(seq_records.seq)
        if length_setting < sequence_length/2:
            length = length_setting
            starting_1 = 1
            ending_1 = length
            starting_2 = len(seq_records.seq) - length
            ending_2 = len(seq_records.seq)
            target_1 = seq_records.seq[starting_1 - 1:ending_1]
            target_2 = seq_records.seq[starting_2 - 1:ending_2]
            id = seq_records.id
            sequences_1 = str(target_1)
            sequences_2 = str(target_2)
            file_name_1 = str(seq_records.id) + '_' + str(starting_1) + '-' + str(ending_1) + '.fasta'
            file_name_2 = str(seq_records.id) + '_' + str(starting_2) + '-' + str(ending_2) + '.fasta'
            output_file_1 = open(output_folder + file_name_1, 'w')
            output_file_2 = open(output_folder + file_name_2, 'w')
            output_file_1.write('>' + str(seq_records.id) + '_' + str(starting_1) + '_' + str(ending_1) + '\n' + sequences_1 + '\n')
            output_file_2.write('>' + str(seq_records.id) + '_' + str(starting_2) + '_' + str(ending_2) + '\n' + sequences_2 + '\n')
            query_seq = output_folder + file_name_1
            subject_seq = output_folder + file_name_2
            output_blast_0 = output_folder + blast_folder_outfmt0 + '/' + str(seq_records.id) + '_blast_outfmt0.txt'
            output_blast_6 = output_folder + blast_folder_outfmt6 + '/' + str(seq_records.id) + '_blast_outfmt6.txt'
            parameters_outfmt0 = ' -evalue 1e-5'
            parameters_outfmt6 = ' -evalue 1e-5 -outfmt 6'
            blast_command_outfmt0 = 'blastn -query ' + query_seq + ' -subject ' + subject_seq + ' -out ' + output_blast_0 + parameters_outfmt0
            blast_command_outfmt6 = 'blastn -query ' + query_seq + ' -subject ' + subject_seq + ' -out ' + output_blast_6 + parameters_outfmt6
            os.system(blast_command_outfmt0)
            os.system(blast_command_outfmt6)
            output_file_1.close()
            output_file_2.close()
            length = 0
        else:
            length = sequence_length/2
            starting_1 = 1
            ending_1 = length
            starting_2 = len(seq_records.seq) - length
            ending_2 = len(seq_records.seq)
            target_1 = seq_records.seq[starting_1 - 1:ending_1]
            target_2 = seq_records.seq[starting_2:ending_2]
            id = seq_records.id
            sequences_1 = str(target_1)
            sequences_2 = str(target_2)
            file_name_1 = str(seq_records.id) + '_' + str(starting_1) + '-' + str(ending_1) + '.fasta'
            file_name_2 = str(seq_records.id) + '_' + str(starting_2) + '-' + str(ending_2) + '.fasta'
            output_file_1 = open(output_folder + file_name_1, 'w')
            output_file_2 = open(output_folder + file_name_2, 'w')
            output_file_1.write('>' + str(seq_records.id) + '_' + str(starting_1) + '_' + str(ending_1) + '\n' + sequences_1 + '\n')
            output_file_2.write('>' + str(seq_records.id) + '_' + str(starting_2) + '_' + str(ending_2) + '\n' + sequences_2 + '\n')
            query_seq = output_folder + file_name_1
            subject_seq = output_folder + file_name_2
            output_blast_0 = output_folder + blast_folder_outfmt0 + '/' + str(seq_records.id) + '_blast_outfmt0.txt'
            output_blast_6 = output_folder + blast_folder_outfmt6 + '/' + str(seq_records.id) + '_blast_outfmt6.txt'
            parameters_outfmt0 = ' -evalue 1e-5'
            parameters_outfmt6 = ' -evalue 1e-5 -outfmt 6'
            blast_command_outfmt0 = 'blastn -query ' + query_seq + ' -subject ' + subject_seq + ' -out ' + output_blast_0 + parameters_outfmt0
            blast_command_outfmt6 = 'blastn -query ' + query_seq + ' -subject ' + subject_seq + ' -out ' + output_blast_6 + parameters_outfmt6
            os.system(blast_command_outfmt0)
            os.system(blast_command_outfmt6)
            output_file_1.close()
            output_file_2.close()
            length = 0

        command_act = act_executable_file + ' ' + query_seq + ' ' + output_blast_6 + ' ' + subject_seq + '\n'
        output_act.write(command_act)
        print(command_act)
