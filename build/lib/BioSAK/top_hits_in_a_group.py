
def keep_top_n_blast_hit(file_in, file_out, top_hit_num):

    file_out_handle = open(file_out, 'w')

    current_query_id = ''
    current_query_wrote_hits = 0
    for blast_hit in open(file_in):

        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]

        if current_query_id == '':
            current_query_id = query_id
            file_out_handle.write(blast_hit)
            current_query_wrote_hits += 1

        elif current_query_id != '':
            if query_id == current_query_id:
                if current_query_wrote_hits < top_hit_num:
                    file_out_handle.write(blast_hit)
                    current_query_wrote_hits += 1
            else:
                current_query_id = query_id
                file_out_handle.write(blast_hit)
                current_query_wrote_hits = 1

    file_out_handle.close()


gtdb_ssu_file = '/Users/songweizhi/Desktop/ssu_all_r202.fna'
blast_op_best_hit = '/Users/songweizhi/Desktop/Pig_1-25_16S_uclust_0.999.QC_vs_GTDB_r202_blastn.tab'


ref_seq_to_taxon_dict = {}
for each_line in open(gtdb_ssu_file):
    if each_line.startswith('>'):

        each_line_split = each_line.strip()[1:].split(' ')
        ref_seq_id = each_line_split[0]

        # get ref_seq_des
        ref_seq_des = ' '.join(each_line_split[1:])
        if ' [' in ref_seq_des:
            ref_seq_des = ref_seq_des.split(' [')[0]
        if ' ' in ref_seq_des:
            ref_seq_des = '_'.join(ref_seq_des.split(' '))

        ref_seq_to_taxon_dict[ref_seq_id] = ref_seq_des

for each_hit in open(blast_op_best_hit):
    each_hit_split  = each_hit.strip().split('\t')
    query_id        = each_hit_split[0]
    ref_id          = each_hit_split[1]
    identity        = float(each_hit_split[2])
    aln_len         = int(each_hit_split[3])
    ref_taxon       = ref_seq_to_taxon_dict.get(ref_id, 'NA')
    if (query_id == 'Pig_subsample_25_5754') and ('Oligosphaerales' in ref_taxon):
    #if (query_id == 'Pig_subsample_25_5354'):
        print('%s\t%s\t%s\t%s\t%s' % (query_id, ref_id, identity, aln_len, ref_taxon))
