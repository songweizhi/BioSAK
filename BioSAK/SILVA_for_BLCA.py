import argparse
from Bio import SeqIO
from BioSAK.global_functions import sep_path_basename_ext


SILVA_for_BLCA_usage = '''
================= SILVA_for_BLCA example commands =================

BioSAK SILVA_for_BLCA -i SILVA_138.2_SSURef_NR99_tax_silva.fasta

# output_file:
SILVA_138.2_SSURef_NR99_tax_silva.blca.fa
SILVA_138.2_SSURef_NR99_tax_silva.blca.tax

===================================================================
'''


def SILVA_for_BLCA(args):

    SILVA_db_file = args['i']

    db_path, db_base, db_ext = sep_path_basename_ext(SILVA_db_file)
    file_out_sequence = '%s/%s.blca.fa'     % (db_path, db_base)
    file_out_taxonomy = '%s/%s.blca.tax'    % (db_path, db_base)

    rank_list = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']

    file_out_sequence_handle = open(file_out_sequence, 'w')
    file_out_taxonomy_handle = open(file_out_taxonomy, 'w')
    for SILVA_seq in SeqIO.parse(SILVA_db_file, 'fasta'):

        SILVA_seq_taxon = ' '.join(SILVA_seq.description.split(' ')[1:])

        if not SILVA_seq_taxon.startswith('Eukaryota'):

            SILVA_seq_taxon_split = SILVA_seq_taxon.split(';')

            if len(SILVA_seq_taxon_split) < 7:

                if len(SILVA_seq_taxon_split) == 6:
                    SILVA_seq_taxon_split.append('unidentified')
                if len(SILVA_seq_taxon_split) == 5:
                    SILVA_seq_taxon_split.append('unidentified')
                    SILVA_seq_taxon_split.append('unidentified')
                if len(SILVA_seq_taxon_split) == 4:
                    SILVA_seq_taxon_split.append('unidentified')
                    SILVA_seq_taxon_split.append('unidentified')
                    SILVA_seq_taxon_split.append('unidentified')
                if len(SILVA_seq_taxon_split) == 3:
                    SILVA_seq_taxon_split.append('unidentified')
                    SILVA_seq_taxon_split.append('unidentified')
                    SILVA_seq_taxon_split.append('unidentified')
                    SILVA_seq_taxon_split.append('unidentified')
                if len(SILVA_seq_taxon_split) == 2:
                    SILVA_seq_taxon_split.append('unidentified')
                    SILVA_seq_taxon_split.append('unidentified')
                    SILVA_seq_taxon_split.append('unidentified')
                    SILVA_seq_taxon_split.append('unidentified')
                    SILVA_seq_taxon_split.append('unidentified')

            if len(SILVA_seq_taxon_split) > 7:
                SILVA_seq_taxon_split = [SILVA_seq_taxon_split[0], SILVA_seq_taxon_split[1], SILVA_seq_taxon_split[2], SILVA_seq_taxon_split[3], SILVA_seq_taxon_split[4], SILVA_seq_taxon_split[5], ' '.join(SILVA_seq_taxon_split[6:])]

            SILVA_seq_taxon_split_reverse = SILVA_seq_taxon_split[::-1]

            SILVA_seq_taxon_str = ''
            n = 0
            for taxon_rank in rank_list:
                SILVA_seq_taxon_str += '%s:%s;' % (taxon_rank, SILVA_seq_taxon_split_reverse[n])
                n += 1

            # write out to taxonomy file
            file_out_taxonomy_handle.write('%s\t%s\n' % (SILVA_seq.id, SILVA_seq_taxon_str))

            # write out to sequence file
            file_out_sequence_handle.write('>%s\n' % SILVA_seq.description)
            file_out_sequence_handle.write('%s\n' % SILVA_seq.seq)

    file_out_sequence_handle.close()
    file_out_taxonomy_handle.close()

    print('Done!')
    print('You may need to run:')
    print('makeblastdb -dbtype nucl -parse_seqids -in %s' % file_out_sequence)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(usage=SILVA_for_BLCA_usage)
    parser.add_argument('-i', required=True, help='SILVA SSU file, e.g. SILVA_138.2_SSURef_NR99_tax_silva.fasta')
    args = vars(parser.parse_args())
    SILVA_for_BLCA(args)
