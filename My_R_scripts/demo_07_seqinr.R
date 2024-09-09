library('seqinr')


fasta_file_in = '/Users/songweizhi/Desktop/seq_demo.fasta'
fasta_file_out = '/Users/songweizhi/Desktop/seq_demo_out.fasta'
ctgs_to_keep = list('seq_1', 'seq_3', 'seq_2')


# initialize fasta file
close(file(fasta_file_out, open = "w"))

for (each_seq in read.fasta(file = fasta_file_in, seqtype = 'DNA', as.string = TRUE, forceDNAtolower = 0)){
    each_seq_id = getName(each_seq)
    
    if (any(ctgs_to_keep == each_seq_id)){
      write.fasta(each_seq, each_seq_id, open = 'a', file.out = fasta_file_out)
    }
}  
  
