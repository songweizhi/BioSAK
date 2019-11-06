import os
import glob
import shutil
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqFeature as SF
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
import multiprocessing as mp
from datetime import datetime


Prodigal_parser_usage = '''
==================== Prodigal_Runner example commands ====================

# for completed genome
BioSAK Prodigal_Runner -i genome_folder -x fa -p KelpGenome -t 6

# for metagenome-assembled genomes (MAGs) 
BioSAK Prodigal_Runner -i bin_folder -x fa -p KelpBins -t 6 -meta

Software dependencies:
module load python/3.5.2
module load prodigal/2.6.3 

==========================================================================
'''


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create)
    os.mkdir(folder_to_create)


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def export_aa_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.protein)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def prodigal_parser(seq_file, sco_file, prefix, output_folder_ffn, output_folder_faa, output_folder_gbk):

    bin_ffn_file =     '%s.ffn' % prefix
    bin_faa_file =     '%s.faa' % prefix
    bin_gbk_file =     '%s.gbk' % prefix
    pwd_bin_ffn_file = '%s/%s'  % (output_folder_ffn, bin_ffn_file)
    pwd_bin_faa_file = '%s/%s'  % (output_folder_faa, bin_faa_file)
    pwd_bin_gbk_file = '%s/%s'  % (output_folder_gbk, bin_gbk_file)

    # get sequence id list
    id_to_sequence_dict = {}
    sequence_id_list = []
    for each_seq in SeqIO.parse(seq_file, 'fasta'):
        id_to_sequence_dict[each_seq.id] = str(each_seq.seq)
        sequence_id_list.append(each_seq.id)


    # get sequence to cds dict and sequence to transl_table dict
    current_seq_id = ''
    current_transl_table = ''
    current_seq_csd_list = []
    seq_to_cds_dict = {}
    seq_to_transl_table_dict = {}
    for each_cds in open(sco_file):
        if each_cds.startswith('# Sequence Data'):

            # add to dict
            if current_seq_id != '':
                seq_to_cds_dict[current_seq_id] = current_seq_csd_list
                seq_to_transl_table_dict[current_seq_id] = current_transl_table

            # reset value
            current_seq_id = each_cds.strip().split('=')[-1][1:-1].split(' ')[0]
            current_transl_table = ''
            current_seq_csd_list = []

        elif each_cds.startswith('# Model Data'):
            current_transl_table = each_cds.strip().split(';')[-2].split('=')[-1]

        else:
            current_seq_csd_list.append('_'.join(each_cds.strip().split('_')[1:]))

    seq_to_cds_dict[current_seq_id] = current_seq_csd_list
    seq_to_transl_table_dict[current_seq_id] = current_transl_table


    bin_gbk_file_handle = open(pwd_bin_gbk_file, 'w')
    bin_ffn_file_handle = open(pwd_bin_ffn_file, 'w')
    bin_faa_file_handle = open(pwd_bin_faa_file, 'w')
    gene_index = 1
    for seq_id in sequence_id_list:

        # create SeqRecord
        current_sequence = Seq(id_to_sequence_dict[seq_id])
        current_SeqRecord = SeqRecord(current_sequence, id=seq_id)
        current_SeqRecord.seq.alphabet = generic_dna
        transl_table = seq_to_transl_table_dict[seq_id]

        current_SeqRecord_annotations = {}
        # current_SeqRecord_annotations['data_file_division'] =   'UNA'
        current_SeqRecord_annotations['date'] =                 (datetime.now().strftime('%d-%b-%Y')).upper()
        # current_SeqRecord_annotations['definition'] =           '%s_%s' % (prefix, seq_id)
        current_SeqRecord_annotations['accession'] =            ''
        current_SeqRecord_annotations['version'] =              ''
        current_SeqRecord_annotations['keywords'] =             ['.']
        current_SeqRecord_annotations['source'] =               prefix
        current_SeqRecord_annotations['organism'] =             prefix
        current_SeqRecord_annotations['taxonomy'] =             ['Unclassified']
        current_SeqRecord_annotations['comment'] =              '.'
        current_SeqRecord.annotations = current_SeqRecord_annotations


        # add SeqFeature to SeqRecord
        for cds in seq_to_cds_dict[seq_id]:

            # define locus_tag id
            locus_tag_id = '%s_%s' % (prefix, "{:0>5}".format(gene_index))

            # define FeatureLocation
            cds_split = cds.split('_')
            cds_start = SF.ExactPosition(int(cds_split[0]))
            cds_end = SF.ExactPosition(int(cds_split[1]))
            cds_strand = cds_split[2]
            current_strand = None
            if cds_strand == '+':
                current_strand = 1
            if cds_strand == '-':
                current_strand = -1
            current_feature_location = FeatureLocation(cds_start, cds_end, strand=current_strand)

            # get nc sequence
            sequence_nc = ''
            if cds_strand == '+':
                sequence_nc = id_to_sequence_dict[seq_id][cds_start-1:cds_end]
            if cds_strand == '-':
                sequence_nc = str(Seq(id_to_sequence_dict[seq_id][cds_start-1:cds_end], generic_dna).reverse_complement())

            # translate to aa sequence
            sequence_aa = str(SeqRecord(Seq(sequence_nc)).seq.translate(table=transl_table))

            # remove * at the end
            sequence_aa = sequence_aa[:-1]

            # export nc and aa sequences
            export_dna_record(sequence_nc, locus_tag_id, '', bin_ffn_file_handle)
            export_aa_record(sequence_aa, locus_tag_id, '', bin_faa_file_handle)

            # Define feature type
            current_feature_type = 'CDS'

            # Define feature qualifiers
            current_qualifiers_dict = {}
            current_qualifiers_dict['locus_tag'] = locus_tag_id
            current_qualifiers_dict['transl_table'] = transl_table
            current_qualifiers_dict['translation'] = sequence_aa

            # Create a SeqFeature
            current_feature = SeqFeature(current_feature_location, type=current_feature_type, qualifiers=current_qualifiers_dict)

            # Append Feature to SeqRecord
            current_SeqRecord.features.append(current_feature)
            gene_index += 1

        # export to gbk file
        SeqIO.write(current_SeqRecord, bin_gbk_file_handle, 'genbank')

    bin_ffn_file_handle.close()
    bin_faa_file_handle.close()
    bin_gbk_file_handle.close()


def prodigal_worker(argument_list):

    input_genome = argument_list[0]
    input_genome_folder = argument_list[1]
    meta_mode = argument_list[2]
    pwd_prodigal_output_folder = argument_list[3]
    pwd_prodigal_output_ffn_folder = argument_list[4]
    pwd_prodigal_output_faa_folder = argument_list[5]
    pwd_prodigal_output_gbk_folder = argument_list[6]

    # prepare command (according to Prokka)
    input_genome_basename, input_genome_ext = os.path.splitext(input_genome)
    pwd_input_genome = '%s/%s' % (input_genome_folder, input_genome)
    pwd_output_sco = '%s/%s.sco' % (pwd_prodigal_output_folder, input_genome_basename)

    prodigal_cmd_meta = 'prodigal -f sco -q -c -m -g 11 -p meta -i %s -o %s' % (pwd_input_genome, pwd_output_sco)
    prodigal_cmd_nonmeta = 'prodigal -f sco -q -c -m -g 11 -i %s -o %s' % (pwd_input_genome, pwd_output_sco)

    if meta_mode is True:
        prodigal_cmd = prodigal_cmd_meta
    else:
        prodigal_cmd = prodigal_cmd_nonmeta

    os.system(prodigal_cmd)

    # prepare ffn, faa and gbk files from prodigal output
    prodigal_parser(pwd_input_genome, pwd_output_sco, input_genome_basename, pwd_prodigal_output_ffn_folder, pwd_prodigal_output_faa_folder, pwd_prodigal_output_gbk_folder)


def Annotation_Prodigal(args):

    input_genome_folder = args['i']
    file_extension = args['x']
    output_prefix = args['p']
    meta_mode = args['meta']
    num_threads = args['t']

    # create output folder
    output_folder_sco = '%s_prodigal_sco' % output_prefix
    output_folder_ffn = '%s_prodigal_ffn' % output_prefix
    output_folder_faa = '%s_prodigal_faa' % output_prefix
    output_folder_gbk = '%s_prodigal_gbk' % output_prefix

    force_create_folder(output_folder_sco)
    force_create_folder(output_folder_ffn)
    force_create_folder(output_folder_faa)
    force_create_folder(output_folder_gbk)

    # get input genome list
    input_genome_re = '%s/*.%s' % (input_genome_folder, file_extension)
    input_genome_file_list = [os.path.basename(file_name) for file_name in glob.glob(input_genome_re)]

    # prepare command list
    list_for_multiple_arguments_Prodigal = []
    for input_genome in input_genome_file_list:
        list_for_multiple_arguments_Prodigal.append([input_genome, input_genome_folder, meta_mode, output_folder_sco, output_folder_ffn, output_folder_faa, output_folder_gbk])

    # run prodigal with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(prodigal_worker, list_for_multiple_arguments_Prodigal)
    pool.close()
    pool.join()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')

    # Annotation modules
    Prodigal_parser = subparsers.add_parser('Prodigal_Runner', description='Wrapper for running Prodigal',usage=Prodigal_parser_usage)

    Prodigal_parser.add_argument('-i',          required=True,  help='input genome folder')
    Prodigal_parser.add_argument('-x',          required=False, default='fasta', help='file extension')
    Prodigal_parser.add_argument('-p',          required=True,  help='output prefix')
    Prodigal_parser.add_argument('-meta',       required=False, action="store_true", help='annotation mode for metagenome assembled genomes (MAGs)')
    Prodigal_parser.add_argument('-t',          required=False, type=int, default=1, help='number of threads')

    args = vars(parser.parse_args())

    Annotation_Prodigal(args)
