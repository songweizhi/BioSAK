from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


gene_seq = 'AAAAAAAAAAAAAAAAAAA'
gene_id = 'Hsp'
gene_description = ''
pwd_output_file = '/Users/songweizhi/Desktop/test_plot_flankings/AAM_00209___BNM_00032/test.fasta'


output_handle = open(pwd_output_file, 'w')
export_dna_record(gene_seq, gene_id, gene_description, output_handle)
output_handle.close()




from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
    fasta_out.write_header()
    fasta_out.write_record(seq_record)
    fasta_out.write_footer()