from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna
from Bio import SeqFeature as SF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


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


def get_genome_length(genome_file):

    total_length_bp = 0
    for each_seq in SeqIO.parse(genome_file, 'fasta'):
        total_length_bp += len(each_seq.seq)

    total_length_Mbp = float("{0:.2f}".format(total_length_bp / (1024 * 1024)))

    return total_length_Mbp


def sco2gbk():
    pass

def sco2ffn():
    pass

def sco2faa():
    pass

def gbk2ffn():
    pass


def gbk2faa():
    pass


def ffn2faa():
    pass

