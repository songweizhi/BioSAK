

def best_hit(args):

    file_in = args['i']
    file_out = args['o']

    file_out_handle = open(file_out, 'w')
    best_hit_line = ''
    best_hit_query_id = ''
    best_hit_score = 0
    for blast_hit in open(file_in):
        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]
        bit_score = float(blast_hit_split[11])

        if best_hit_query_id == '':
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

        elif (query_id == best_hit_query_id) and (bit_score > best_hit_score):
            best_hit_score = bit_score
            best_hit_line = blast_hit

        elif query_id != best_hit_query_id:
            file_out_handle.write(best_hit_line)
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

    file_out_handle.write(best_hit_line)
    file_out_handle.close()


if __name__ == '__main__':

    file_in = '/Users/songweizhi/Desktop/ResistFlowBins/class_level/ResistFlow_c37_HGTs_PG.tab'
    file_out = '/Users/songweizhi/Desktop/ResistFlowBins/class_level/ResistFlow_c37_HGTs_PG_best_hit.tab'

    args = [file_in, file_out]

    best_hit(args)
