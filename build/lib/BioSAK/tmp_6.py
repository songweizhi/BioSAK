from Bio import SeqIO

def fq2fa(fq_in, fa_out):
    fa_out_handle = open(fa_out, 'w')
    for each_long_read in SeqIO.parse(fq_in, 'fastq'):
        fa_out_handle.write('>%s\n%s\n' % (each_long_read.id, each_long_read.seq))
    fa_out_handle.close()


fq_in = '/Users/songweizhi/Desktop/abd_demo/AphrocallistesBeatrix_subset.fastq'
fa_out = '/Users/songweizhi/Desktop/abd_demo/AphrocallistesBeatrix_subset.fasta'
fq2fa(fq_in, fa_out)
