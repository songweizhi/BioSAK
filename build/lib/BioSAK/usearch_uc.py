import argparse


usearch_uc_usage = '''
=================== usearch_uc example commands ===================

BioSAK usearch_uc -uc uclust_0.6.uc -o uclust_0.6_min1.txt
BioSAK usearch_uc -uc uclust_0.6.uc -o uclust_0.6_min3.txt -n 3

===================================================================
'''


def usearch_uc(args):

    uc_file     = args['uc']
    min_seq_num = args['n']
    output_file = args['o']

    cluster_id_set = set()
    cluster_to_seq_member_dict = {}
    for each_line in open(uc_file):
        each_line_split = each_line.strip().split('\t')
        cluster_id = each_line_split[1]
        seq_id = each_line_split[8].split(' ')[0]
        cluster_id_set.add(int(cluster_id))

        if cluster_id not in cluster_to_seq_member_dict:
            cluster_to_seq_member_dict[cluster_id] = {seq_id}
        else:
            cluster_to_seq_member_dict[cluster_id].add(seq_id)

    # write out cluster sequence members
    qualified_cluster_num = 0
    output_file_handle = open(output_file, 'w')
    for each_cluster in sorted([i for i in cluster_id_set]):
        current_cluster_seqs = cluster_to_seq_member_dict[str(each_cluster)]
        current_cluster_seqs_sorted = sorted([i for i in current_cluster_seqs])

        if len(current_cluster_seqs) >= min_seq_num:
            output_file_handle.write('Cluster_%s\t%s\n' % (each_cluster, ','.join(current_cluster_seqs_sorted)))
            qualified_cluster_num += 1
    output_file_handle.close()

    if min_seq_num == 1:
        print('Found %s clusters with at least %s sequence.' % (qualified_cluster_num, min_seq_num))
    if min_seq_num > 1:
        print('Found %s clusters with at least %s sequences.' % (qualified_cluster_num, min_seq_num))


if __name__ == '__main__':

    usearch_uc_parser = argparse.ArgumentParser()
    usearch_uc_parser.add_argument('-uc', required=True,                        help='uc file from Usearch')
    usearch_uc_parser.add_argument('-n',  required=False, type=int, default=1,  help='minimum number of sequence in a cluster, default: 1')
    usearch_uc_parser.add_argument('-o',  required=True,                        help='output file')
    args = vars(usearch_uc_parser.parse_args())
    usearch_uc(args)

