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
==================== Barrnap_Runner example commands ====================

# module needed
module load perl/5.28.0
module load hmmer/3.2.1
module load bedtools/2.27.1
module load barrnap/0.9

# for completed genome
BioSAK Barrnap_Runner -i genome_folder -x fna -t 6


barrnap MBAD.fna > barrnap_MBAD.gff3
bedtools getfasta -fi MBAD.fna -bed barrnap_MBAD.gff3 -fo rRNA_seq_MBAD.fasta -name -s

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


def Annotation_Barrnap(args):
    pass


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')

    # Annotation modules
    Prodigal_parser = subparsers.add_parser('Prodigal_Runner', description='Wrapper for running Prodigal',usage=Prodigal_parser_usage)

    Prodigal_parser.add_argument('-i',          required=True,  help='input genome folder')
    Prodigal_parser.add_argument('-x',          required=False, default='fna', help='file extension')
    Prodigal_parser.add_argument('-t',          required=False, type=int, default=1, help='number of threads')

    args = vars(parser.parse_args())

    Annotation_Barrnap(args)
