
def subset_prokaryotes_csv(GTDB_genomes, GenBank_genomes, GenBank_genomes_subset):

    num_to_full_id_dict = {}
    GTDB_genome_id_list = set()
    for each in open(GTDB_genomes):
        if not each.startswith('"ID",'):
            each_split = each.strip().split(',')
            NCBI_ID_with_version = each_split[0][1:-1]
            NCBI_ID_without_version = '.'.join(NCBI_ID_with_version.split('.')[:-1])
            NCBI_ID_no_version_only_number = NCBI_ID_without_version[3:]
            GTDB_genome_id_list.add(NCBI_ID_no_version_only_number)
            num_to_full_id_dict[NCBI_ID_no_version_only_number] = NCBI_ID_with_version

    GenBank_genomes_subset_handle = open(GenBank_genomes_subset, 'w')
    found_in_GenBank_genome_id_list = set()
    for GenBank_genome in open(GenBank_genomes):
        if not GenBank_genome.startswith('#Organism Name'):
            GenBank_genome_split = GenBank_genome.strip().split(',')
            GenBank_genome_id_with_version = GenBank_genome_split[5][1:-1]
            GenBank_genome_id_without_version = '.'.join(GenBank_genome_id_with_version.split('.')[:-1])
            GenBank_genome_id_no_version_only_number = GenBank_genome_id_without_version[3:]
            if GenBank_genome_id_no_version_only_number in GTDB_genome_id_list:
                found_in_GenBank_genome_id_list.add(GenBank_genome_id_no_version_only_number)
                GenBank_genomes_subset_handle.write(GenBank_genome)
    GenBank_genomes_subset_handle.close()

    for each in GTDB_genome_id_list:
        if each not in found_in_GenBank_genome_id_list:
            print('not found: %s' % num_to_full_id_dict[each])


GTDB_genomes           = '/Users/songweizhi/Desktop/GTDB.csv'
GenBank_genomes        = '/Users/songweizhi/Desktop/prokaryotes.csv'
GenBank_genomes_subset = '/Users/songweizhi/Desktop/prokaryotes_subset.csv'

subset_prokaryotes_csv(GTDB_genomes, GenBank_genomes, GenBank_genomes_subset)

