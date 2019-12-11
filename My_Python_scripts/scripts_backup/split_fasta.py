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


#output_handle = open(pwd_output_file, 'w')

# get assignment dict
needed_ctgs = []
for each_assignment in open('/Users/songweizhi/Desktop/DSM17395.haplotigs/DSM17395.haplotigs.txt'):
    print(each_assignment)
    ctg_id = each_assignment.strip().split('\t')[0]
    ctg_assign = each_assignment.strip().split('\t')[1]
    if ctg_assign == '2.10':
        needed_ctgs.append(ctg_id)

print(needed_ctgs)


for each in SeqIO.parse('/Users/songweizhi/Desktop/DSM17395.haplotigs/DSM17395.haplotigs.fas', 'fasta'):
    if each.id in needed_ctgs:
        print(each.id)
        output_handle = open('/Users/songweizhi/Desktop/DSM17395.haplotigs/%s.fasta' % each.id, 'w')
        export_dna_record(str(each.seq), each.id, '', output_handle)

        output_handle.close()
