import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


def get_recovered_iden100_transfers(d, match_summary):
    # read distribution into a dictionary
    transfer_to_recipient_dict = {}
    for each in open(d):
        each_split = each.strip().split(',')
        recipient_genome = each_split[0]
        transfers = each_split[1:]
        for transfer in transfers:
            transfer_to_recipient_dict[transfer] = recipient_genome

    recovered_num = 0
    recovered_transfers = []
    for match in open(match_summary):
        match_split = match.strip().split('\t')
        transfer_d = match_split[0]
        recipient_genome_d = match_split[1:]
        if transfer_to_recipient_dict[transfer_d] in recipient_genome_d:
            recovered_num += 1
            #print(transfer_d)
            recovered_transfers.append(transfer_d)

    return recovered_num, recovered_transfers


parser = argparse.ArgumentParser()

parser.add_argument('-a',
                    required=True,
                    help='Assemblies in multi-fasta format')

parser.add_argument('-t',
                    required=True,
                    help='Sequences of transferred genes in multi-fasta format')

parser.add_argument('-minf',
                    required=False,
                    default=1000,
                    type=int,
                    help='minimum length of flanking sequences')

parser.add_argument('-d',
                    required=True,
                    help='distribution of transfers to the recipient genomes')

parser.add_argument('-combined_genomes',
                    required=True,
                    help='distribution of transfers to the recipient genomes')

args = vars(parser.parse_args())
scaffold_file = args['a']
transfers_fasta = args['t']
minf = args['minf']
transfer_profile_file = args['d']
combined_genomes = args['combined_genomes']


# make blast database
wd = os.getcwd()
pwd_transfers_fasta = '%s/%s' % (wd, transfers_fasta)
pwd_scaffold_file = '%s/%s' % (wd, scaffold_file)
pwd_makeblastdb_exe = 'makeblastdb'
pwd_blastn_exe = 'blastn'
pwd_blast_result = '%s/blast_result.tab' % wd
os.system('%s -in %s -dbtype nucl -parse_seqids' % (pwd_makeblastdb_exe, pwd_scaffold_file))
outfmt = '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"'
print('\nRunning Blast, be patient...')
os.system('%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, pwd_transfers_fasta, pwd_scaffold_file, pwd_blast_result, outfmt))

matches = open('%s/matches.tab' % wd, 'w')


blast_results = open(pwd_blast_result)
recovered_transfers = []
#matches.write('query\tsubject\tleft_b\tleft_e\tright_b\tright_e\tleft_len\tright_len\n')
for match in blast_results:
    match_split = match.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    identity = float(match_split[2])
    align_len = int(match_split[3])
    sstart = int(match_split[8])
    send = int(match_split[9])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    coverage_q = float("{0:.2f}".format(float(align_len) * 100 / float(query_len)))
    coverage_s = float("{0:.2f}".format(float(align_len) * 100 / float(subject_len)))

    if send > sstart:
        if (identity >= 99) and (coverage_q >= 99) and ((sstart >= minf) or (subject_len - send >= minf)):
            #print(query)
            left_flanking_b = sstart - minf
            if left_flanking_b <= 0:
                left_flanking_b = 1
            left_flanking_e = sstart
            right_flanking_b = send
            right_flanking_e = send + minf
            if right_flanking_e > subject_len:
                right_flanking_e = subject_len
            left_flanking = [left_flanking_b, left_flanking_e]
            right_flanking = [right_flanking_b, right_flanking_e]

            left_flanking_length = left_flanking_e - left_flanking_b
            right_flanking_length = right_flanking_e - right_flanking_b
            matches.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (query, subject, left_flanking_b, left_flanking_e, right_flanking_b, right_flanking_e, left_flanking_length, right_flanking_length))
            if query not in recovered_transfers:
                recovered_transfers.append(query)
    if send < sstart:
        if (identity >= 99) and (coverage_q >= 99) and ((send >= minf) or (subject_len - sstart >= minf)):
            left_flanking_e = send
            left_flanking_b = send - minf
            if left_flanking_b <= 0:
                left_flanking_b = 1
            right_flanking_b = sstart
            right_flanking_e = sstart + minf
            if right_flanking_e > subject_len:
                right_flanking_e = subject_len
            left_flanking = [left_flanking_b, left_flanking_e]
            right_flanking = [right_flanking_b, right_flanking_e]
            left_flanking_length = left_flanking_e - left_flanking_b
            right_flanking_length = right_flanking_e - right_flanking_b
            matches.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (query, subject, left_flanking_b, left_flanking_e, right_flanking_b, right_flanking_e, left_flanking_length, right_flanking_length))
            if query not in recovered_transfers:
                recovered_transfers.append(query)
matches.close()


matches_r = open('%s/matches.tab' % wd)
matched_list = []
for each in matches_r:
    each_split = each.strip().split('\t')
    query_m = each_split[0]
    subject_m = each_split[1]
    left_b_m = int(each_split[2])
    left_e_m = int(each_split[3])
    right_b_m = int(each_split[4])
    right_e_m = int(each_split[5])
    left_len_m = int(each_split[6])
    right_len_m = int(each_split[7])
    matched_list.append(query_m)
    # read in assemblies
    for each_contig in SeqIO.parse(pwd_scaffold_file, 'fasta'):
        contig_id = each_contig.id
        if contig_id == subject_m:
            #print('\n%s\t%s' % (query_m, each_contig.id))
            flank_handle = open('%s/%s_flanking.fasta' % (wd, query_m), 'w')
            each_contig_seq = str(each_contig.seq)
            left_flanking_seq = each_contig_seq[left_b_m - 1:left_e_m - 1]
            right_flanking_seq = each_contig_seq[right_b_m:right_e_m]
            if len(left_flanking_seq) >= 50:
                left_flanking_seq_object = Seq(left_flanking_seq, generic_dna)
                left_flanking_seq_record = SeqRecord(left_flanking_seq_object)
                left_flanking_seq_record.id = '%slflk' % query_m
                left_flanking_seq_record.description = ''
                SeqIO.write(left_flanking_seq_record, flank_handle, 'fasta')
            if len(right_flanking_seq) >= 50:
                right_flanking_seq_object = Seq(right_flanking_seq, generic_dna)
                right_flanking_seq_record = SeqRecord(right_flanking_seq_object)
                right_flanking_seq_record.id = '%srflk' % query_m
                right_flanking_seq_record.description = ''
                SeqIO.write(right_flanking_seq_record, flank_handle, 'fasta')

            flank_handle.close()
            # print('left_flanking(%s):\n%s' % (len(left_flanking_seq), left_flanking_seq))
            # print('right_flanking(%s):\n%s' % (len(right_flanking_seq), right_flanking_seq))

# makeblastdb with combined genomes
pwd_combined_genomes = '%s/%s' % (wd, combined_genomes)
os.system('%s -in %s -dbtype nucl -parse_seqids' % (pwd_makeblastdb_exe, pwd_combined_genomes))

pwd_matched_summary = '%s/matched_summary.tab' % wd
matched_summary_handle = open(pwd_matched_summary, 'w')
for each_match in matched_list:
    #print(each_match)
    flank_file = '%s/%s_flanking.fasta' % (wd, each_match)
    blastdb = pwd_combined_genomes

    #print(flank_file)
    #print(combined_genomes)
    pwd_each_match_blast_result = '%s/%s_flanking.tab' % (wd, each_match)
    outfmt = '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"'
    os.system('%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, flank_file, blastdb, pwd_each_match_blast_result, outfmt))

    # read in blast results
    current_match_result = open(pwd_each_match_blast_result)
    current_match_list = []
    for each_current_match in current_match_result:
        #print(each_current_match)
        each_current_match_split = each_current_match.strip().split('\t')
        subject_genome = each_current_match_split[1]
        identity = float(each_current_match_split[2])
        if identity == 100:
            current_match_list.append(subject_genome)
    #print(current_match_list)
    current_match_list_uniq = []
    for each_match_genome in current_match_list:
        if each_match_genome not in current_match_list_uniq:
            current_match_list_uniq.append(each_match_genome)
    #print(current_match_list_uniq)
    #print('\t'.join(current_match_list_uniq))

    matched_summary_handle.write('%s\t%s\n' % (each_match, '\t'.join(current_match_list_uniq)))
    # if len(current_match_list_uniq) > 1:
    #     print(each_match)
matched_summary_handle.close()


recovered_iden100_transfer_num, recovered_iden100_transfer_ids = get_recovered_iden100_transfers(transfer_profile_file, pwd_matched_summary)

print('Recovered transfers(identity >= 99, coverage >= 99 and at least one flanking side longer than %sbp): %s' % (minf, recovered_iden100_transfer_num))

print(recovered_iden100_transfer_ids)
print(len(recovered_iden100_transfer_ids))



########## get length distribution ##########
list_lt600 = []
list_lt800 = []
list_lt1000 = []
list_lt1200 = []
list_lt1400 = []
list_lt1600 = []
list_lt1800 = []
list_lt2000 = []
for each in SeqIO.parse(pwd_transfers_fasta, 'fasta'):
    seq_id = each.id
    seq_length = len(each.seq)
    if 600 <=seq_length < 800:
        list_lt600.append(seq_id)
    elif 800 <=seq_length < 1000:
        list_lt800.append(seq_id)
    elif 1000 <=seq_length < 1200:
        list_lt1000.append(seq_id)
    elif 1200 <=seq_length < 1400:
        list_lt1200.append(seq_id)
    elif 1400 <=seq_length < 1600:
        list_lt1400.append(seq_id)
    elif 1600 <=seq_length < 1800:
        list_lt1600.append(seq_id)
    elif 1800 <=seq_length < 2000:
        list_lt1800.append(seq_id)
    elif seq_length >= 2000:
        list_lt2000.append(seq_id)

# get recovered genes for each group
list_lt600_recovered = []
list_lt800_recovered = []
list_lt1000_recovered = []
list_lt1200_recovered = []
list_lt1400_recovered = []
list_lt1600_recovered = []
list_lt1800_recovered = []
list_lt2000_recovered = []
print(recovered_iden100_transfer_ids)
print(len(recovered_iden100_transfer_ids))
for each in recovered_iden100_transfer_ids:
    if each in list_lt600:
        list_lt600_recovered.append(each)
    elif each in list_lt800:
        list_lt800_recovered.append(each)
    elif each in list_lt1000:
        list_lt1000_recovered.append(each)
    elif each in list_lt1200:
        list_lt1200_recovered.append(each)
    elif each in list_lt1400:
        list_lt1400_recovered.append(each)
    elif each in list_lt1600:
        list_lt1600_recovered.append(each)
    elif each in list_lt1800:
        list_lt1800_recovered.append(each)
    elif each in list_lt2000:
        list_lt2000_recovered.append(each)

print('600-800bp: %s' % (round(len(list_lt600_recovered)/len(list_lt600), 2)))
print('800-1000bp: %s' % (round(len(list_lt800_recovered)/len(list_lt800), 2)))
print('1000-1200bp: %s' % (round(len(list_lt1000_recovered)/len(list_lt1000), 2)))
print('1200-1400bp: %s' % (round(len(list_lt1200_recovered)/len(list_lt1200), 2)))
print('1400-1600bp: %s' % (round(len(list_lt1400_recovered)/len(list_lt1400), 2)))
print('1600-1800bp: %s' % (round(len(list_lt1600_recovered)/len(list_lt1600), 2)))
print('1800-2000bp: %s' % (round(len(list_lt1800_recovered)/len(list_lt1800), 2)))
print('>2000bp: %s' % (round(len(list_lt2000_recovered)/len(list_lt2000), 2)))

