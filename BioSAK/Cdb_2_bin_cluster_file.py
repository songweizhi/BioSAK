
def Cdb_2_bin_cluster_file(Cdb_file, bin_cluster_file):
    cluster_to_bin_dict = {}
    obtained_clusters = set()
    for each_bin in open(Cdb_file):
        if not each_bin.startswith('genome,'):
            each_bin_split = each_bin.strip().split(',')
            bin_id = each_bin_split[0]
            secondary_cluster = each_bin_split[-1]
            obtained_clusters.add(secondary_cluster)
            if secondary_cluster not in cluster_to_bin_dict:
                cluster_to_bin_dict[secondary_cluster] = [bin_id]
            else:
                cluster_to_bin_dict[secondary_cluster].append(bin_id)

    obtained_clusters_list = sorted([i for i in obtained_clusters])

    bin_cluster_file_handle = open(bin_cluster_file, 'w')
    for j in obtained_clusters_list:
        bin_cluster_file_handle.write('cluster_%s\t%s\n' % (j, '\t'.join(cluster_to_bin_dict[j])))
    bin_cluster_file_handle.close()


Cdb_2_bin_cluster_file('/Users/songweizhi/Desktop/Cdb_ANI95.csv', '/Users/songweizhi/Desktop/Cdb_ANI95.csv.txt')
Cdb_2_bin_cluster_file('/Users/songweizhi/Desktop/Cdb_ANI97.csv', '/Users/songweizhi/Desktop/Cdb_ANI97.csv.txt')
