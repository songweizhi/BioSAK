import os
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


os.chdir('/Users/songweizhi/Desktop/ccc')

blast_result = 'BS107_ref.txt'

for each in open(blast_result):
    #print(each)
    each_split = each.strip().split('\t')
    aln_length = int(each_split[3])
    if aln_length >= 5000:
        print(each.strip())






# out_handle = open('tig1_aligned.fasta', 'w')
#
# for each in SeqIO.parse('BS107_canu.fasta', 'fasta'):
#     if each.id == 'tig00000001':
#         each_aligned = each.seq[32114:(621808 - 580000)]
#         print(len(each_aligned))
#         export_dna_record(str(each_aligned), each.id, '', out_handle)
#
#
# for each2 in SeqIO.parse('BS107_ref.fasta', 'fasta'):
#     if each2.id == 'BS107_641382347':
#         each2_aligned = each2.seq[0:(589711 - 580000)]
#         print(len(each2_aligned))
#         export_dna_record(str(each2_aligned), each2.id, '', out_handle)
#
# for each3 in SeqIO.parse('DSM17395_ref.fasta', 'fasta'):
#     if each3.id == 'NC_018290.1':
#         each3_aligned = each3.seq[1253594:(1843303 - 580000)]
#         print(len(each3_aligned))
#         export_dna_record(str(each3_aligned), each3.id, '', out_handle)
#
# out_handle.close()
#
# # tig00000001		32115	    621809
# # BS107_641382347	1		    589711
# # NC_018290		    1189375		1843304


