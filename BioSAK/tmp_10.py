from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

seq = 'ATGC'

seq = Seq(seq).reverse_complement()
print(seq)