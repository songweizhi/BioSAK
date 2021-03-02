import os
import argparse
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


sequence_manipulator_usage = '''
================================ Sequence manipulator example commands ===============================

BioSAK gbk2fa -gbk bin_1.gbk
BioSAK gbk2ffn -gbk bin_1.gbk
BioSAK gbk2faa -gbk bin_1.gbk
BioSAK ffn2faa -ffn bin_1.ffn

# get reverse complement sequence(s) for single and multiple sequence(s)
BioSAK get_rc -seq AAAAATTTTTGGGGGCCCCC
BioSAK get_rc -seq bin_1.ffn

======================================================================================================
'''


convert_align_format_usage = '''
================================ convert_align_format example commands ===============================

# convert alignment in fasta format to phylip-relaxed format
BioSAK convert_align_format -in NorthSea.aln -inf fasta -out NorthSea.phylip -outf phylip-relaxed

# Alignment format:
  clustal, emboss, fasta, fasta-m10, ig, maf, mauve, nexus, 
  phylip, phylip-sequential, phylip-relaxed, stockholm

# More details about alignment format is here: 
  https://biopython.org/wiki/AlignIO

======================================================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def export_seq_record(gene_seq, gene_id, gene_description, seq_type, output_handle):

    seq_object = Seq(gene_seq)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def gbk2fa(args):

    gbk_in = args['gbk']
    gbk_path, gbk_basename, gbk_extension = sep_path_basename_ext(gbk_in)
    ffn_out = '%s/%s.fa' % (gbk_path, gbk_basename)
    SeqIO.convert(gbk_in, 'genbank', ffn_out, 'fasta')


def gbk2ffn(args):

    gbk_in = args['gbk']
    gbk_path, gbk_basename, gbk_extension = sep_path_basename_ext(gbk_in)
    ffn_out = '%s/%s.ffn' % (gbk_path, gbk_basename)

    ffn_out_handle = open(ffn_out, 'w')
    for seq_record in SeqIO.parse(gbk_in, 'genbank'):
        for feature in seq_record.features:
            feature_locus_tag = feature.qualifiers['locus_tag'][0]

            if feature.location.strand == 1:
                feature_seq = str(seq_record.seq)[feature.location.start-1:feature.location.end]
            else:
                feature_seq_rc = str(seq_record.seq)[feature.location.start-1:feature.location.end]
                feature_seq = str(Seq(feature_seq_rc).reverse_complement())

            export_seq_record(feature_seq, feature_locus_tag, '', 'N', ffn_out_handle)

    ffn_out_handle.close()


def gbk2faa(args):

    gbk_in = args['gbk']
    gbk_path, gbk_basename, gbk_extension = sep_path_basename_ext(gbk_in)
    faa_out = '%s/%s.faa' % (gbk_path, gbk_basename)

    faa_out_handle = open(faa_out, 'w')
    for seq_record in SeqIO.parse(gbk_in, 'genbank'):
        for feature in seq_record.features:
            feature_locus_tag = feature.qualifiers['locus_tag'][0]
            feature_translation = feature.qualifiers['translation'][0]
            export_seq_record(feature_translation, feature_locus_tag, '', 'P', faa_out_handle)

    faa_out_handle.close()


def ffn2faa(args):

    ffn_in = args['ffn']
    ffn_path, ffn_basename, ffn_extension = sep_path_basename_ext(ffn_in)
    faa_out = '%s/%s.faa' % (ffn_path, ffn_basename)

    faa_out_handle = open(faa_out, 'w')
    for seq_record in SeqIO.parse(ffn_in, 'fasta'):
        seq_record_seq = seq_record.seq
        seq_record_seq_aa = seq_record_seq.translate()
        export_seq_record(str(seq_record_seq_aa), str(seq_record.id), '', 'P', faa_out_handle)
    faa_out_handle.close()


def get_rc(args):

    seq_in = args['seq']

    # check whether seq_in is a file
    if os.path.isfile(seq_in) is True:

        seq_in_path, seq_in_basename, seq_in_extension = sep_path_basename_ext(seq_in)
        seq_out = '%s/%s_rc%s' % (seq_in_path, seq_in_basename, seq_in_extension)

        seq_out_handle = open(seq_out, 'w')
        for seq_record in SeqIO.parse(seq_in, 'fasta'):
            seq_record_rc_id = '%s_rc' % seq_record.id
            seq_record_seq = seq_record.seq
            seq_record_seq_rc = seq_record_seq.reverse_complement()
            export_seq_record(str(seq_record_seq_rc), seq_record_rc_id, '', 'N', seq_out_handle)
        seq_out_handle.close()

    else:
        # check whether input is DNA sequence
        nc_bases = ['A', 'a', 'T', 't', 'G', 'g', 'C', 'c']
        nc_bases_count = 0
        for nc_base in nc_bases:
            nc_bases_count += seq_in.count(nc_base)

        if nc_bases_count == len(seq_in):
            seq_in_rc = str(Seq(seq_in).reverse_complement())
            print('>seq_in_rc\n%s\n' % seq_in_rc)

        else:
            print('Sequence file not exist or non DNA sequence detected, program exited!')
            exit()


def fq2fa(args):

    fq_in = args['fq']
    fq_path, fq_basename, fq_extension = sep_path_basename_ext(fq_in)
    fa_out = '%s/%s.fa' % (fq_path, fq_basename)

    SeqIO.convert(fq_in, 'fastq', fa_out, 'fastq')


def convert_align_format(args):

    aln_in =            args['in']
    aln_in_format =     args['inf']
    aln_out =           args['out']
    aln_out_format =    args['outf']

    AlignIO.convert(aln_in, aln_in_format, aln_out, aln_out_format)

