# Extract sequences with provided id
BioSAK select_seq -seq ctg.fasta -id seq_id.txt -out output.1.fa -option 1 
BioSAK select_seq -seq reads.fastq -id read_id.txt -out output.1.fq -option 1 -fq

# Extract sequences except those in seq_id.txt
BioSAK select_seq -seq ctg.fasta -id seq_id.txt -out output.0.fa -option 0
BioSAK select_seq -seq ctg.fasta -id seq_id.txt -out output.0.fa -option 0 -oneline

# seq_id.txt file format: one id per line, great than symbol excluded.
