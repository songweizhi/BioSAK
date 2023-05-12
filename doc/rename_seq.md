
+ Rename "NODE_941_length_17600_cov_52.7123" to "NODE_941"

      BioSAK rename_seq -in Contigs.fa -sep_in "_" -n 2

+ Rename "Seawater|NODE|941|length|17600|cov|52.7123" to "Seawater_NODE_941"

      BioSAK rename_seq -in Contigs.fa -sep_in "|" -sep_out "_" -n 3

+ Add prefix to all sequences in a fasta file

      BioSAK rename_seq -in Contigs.fa -prefix seawater

+ Rename "NODE_941_length_17600_cov_52.7123" to "Seawater_NODE_941"

      BioSAK rename_seq -in Contigs.fa -sep_in "_" -n 2 -prefix Seawater

+ Rename sequences in multiple files

      BioSAK rename_seq -in seq_folder -x fa -sep_in "_" -n 2 -t 12
      BioSAK rename_seq -in seq_folder -x fa -prefix prefix_file.txt

+ Prefix file format: sequence file name (with extension) followed by prefixed to be added, tab separated, one file per line
  
      genome_1.fa prefix_1
      genome_1.fa prefix_2
      genome_1.fa prefix_3
