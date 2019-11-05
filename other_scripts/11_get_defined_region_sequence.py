import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


# usage example
# python3 11_get_defined_region_sequence.py -f input.fasta -i contig_1 -s 5000 -e 10000 -o output.fasta
# python3 /Users/songweizhi/PycharmProjects/python_scripts/11_get_defined_region_sequence.py -f /Users/songweizhi/Desktop/FC/combined_ref.fasta -i 2.10_chromosome -s 866221 -e 941272 -o /Users/songweizhi/Desktop/2.10_chromosome_866221-941272.fasta


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


parser = argparse.ArgumentParser(description='', add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
optional.add_argument('-h', action='help', help='Show this help message and exit')
required.add_argument('-f', dest='FILE', nargs='?', required=True,  type=str, help='Sequence file')
required.add_argument('-i', dest='SEQ_ID', nargs='?', required=True, type=str, help='Sequence ID')
required.add_argument('-s', dest='START_POS', nargs='?', required=True,  type=int, help='Start position')
required.add_argument('-e', dest='END_POS', nargs='?', required=True,  type=int, help='End position')
required.add_argument('-o', dest='OUTPUT', nargs='?', required=True,  type=str, help='End position')

args = vars(parser.parse_args())
ref_genome = args['FILE']
ref_seq = args['SEQ_ID']
start_pos = args['START_POS']
end_pos = args['END_POS']
output_file = args['OUTPUT']


output_file_handle = open(output_file, 'w')
for each_seq in SeqIO.parse(ref_genome, 'fasta'):
    if each_seq.id == ref_seq:
        seq_to_export = str(each_seq.seq[start_pos-1:end_pos])
        output_seq_id = '%s_%s_%s' % (ref_seq, start_pos, end_pos)
        export_dna_record(seq_to_export, output_seq_id, '', output_file_handle)
output_file_handle.close()

