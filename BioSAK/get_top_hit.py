import argparse


def get_best_hit_all(file_in, file_out, taxonomy_dict):

    highest_bit_score_dict = {}
    for blast_hit in open(file_in):
        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]
        bit_score = float(blast_hit_split[12])

        if query_id not in highest_bit_score_dict:
            highest_bit_score_dict[query_id] = bit_score
        else:
            if bit_score > highest_bit_score_dict[query_id]:
                highest_bit_score_dict[query_id] = bit_score

    file_out_handle = open(file_out, 'w')
    file_out_handle.write('X.OTU.ID\tRef.ID\tIdentity\tAlignment_Length\tEvalue\tBit_Score\ttaxonomy\n')
    for blast_hit in open(file_in):
        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]
        ref_id = blast_hit_split[1]
        identity = float(blast_hit_split[2])
        Alignment_Length = int(blast_hit_split[3])
        evalue = float(blast_hit_split[10])
        bit_score = float(blast_hit_split[12])

        if bit_score == highest_bit_score_dict[query_id]:
            file_out_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (query_id, ref_id, identity, Alignment_Length, evalue, bit_score, taxonomy_dict[ref_id]))

    file_out_handle.close()


def get_best_hit_first(file_in, file_out, taxonomy_dict):

    file_out_handle = open(file_out, 'w')
    file_out_handle.write('X.OTU.ID\tRef.ID\tIdentity\tAlignment_Length\tEvalue\tBit_Score\ttaxonomy\n')
    best_hit_line = ''
    best_hit_query_id = ''
    best_hit_score = 0
    for blast_hit in open(file_in):
        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]
        bit_score = int(blast_hit_split[12])

        if best_hit_query_id == '':
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

        elif (query_id == best_hit_query_id) and (bit_score > best_hit_score):
            best_hit_score = bit_score
            best_hit_line = blast_hit

        elif query_id != best_hit_query_id:
            best_hit_line_split = best_hit_line.strip().split('\t')
            file_out_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (best_hit_line_split[0], best_hit_line_split[1], best_hit_line_split[2], best_hit_line_split[3], best_hit_line_split[10], float(best_hit_line_split[12]), taxonomy_dict[best_hit_line_split[1]]))
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

    best_hit_line_split = best_hit_line.strip().split('\t')
    file_out_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (best_hit_line_split[0], best_hit_line_split[1], best_hit_line_split[2], best_hit_line_split[3],best_hit_line_split[10], float(best_hit_line_split[12]), taxonomy_dict[best_hit_line_split[1]]))
    file_out_handle.close()


def get_best_hit(args):

    file_in =       args['i']
    file_out =      args['o']
    taxonomy_file = args['t']
    keep_all =      args['all']

    # store taxonomy info in dict
    taxonomy_dict = {}
    for ref in open(taxonomy_file):
        ref_split = ref.strip().split('\t')
        ref_id = ref_split[0]
        ref_taxon = ref_split[1]
        ref_taxon_split_no_rank = [i.split(':')[1] for i in ref_taxon.split(';')[:-1]]
        ref_taxon_no_rank = ';'.join(ref_taxon_split_no_rank[::-1])
        taxonomy_dict[ref_id] = ref_taxon_no_rank

    if keep_all is True:
        get_best_hit_all(file_in, file_out, taxonomy_dict)
    else:
        get_best_hit_first(file_in, file_out, taxonomy_dict)


if __name__ == '__main__':

    get_best_hit_parser = argparse.ArgumentParser()

    # arguments for COG_parser
    get_best_hit_parser.add_argument('-i',   required=True,                          help='blast results from BLCA')
    get_best_hit_parser.add_argument('-o',   required=True,                          help='output file')
    get_best_hit_parser.add_argument('-t',   required=True,                          help='taxonomy file')
    get_best_hit_parser.add_argument('-all', required=False, action='store_true',    help='keep all blast hits with highest bit score')

    args = vars(get_best_hit_parser.parse_args())

    get_best_hit(args)

