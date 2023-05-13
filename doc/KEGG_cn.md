
中文版
---


### Dependencies: 


blast+ or diamond

### Annotation with NCBI blastp (default, for small dataset)

    BioSAK KEGG -db_dir path/to/your/KEGG_db_dir -t 6 -seq_in input.faa -depth input.depth

### Annotation with Diamond blastp (for big dataset)

    BioSAK KEGG -db_dir path/to/your/KEGG_db_dir -t 12 -seq_in faa_folder -x faa -depth depth_files -diamond

### Get summary for BlastKOALA/GhostKOALA produced results

    BioSAK KEGG -db_dir path/to/your/KEGG_db_dir -t 9 -ko_in user_ko.txt
    BioSAK KEGG -db_dir path/to/your/KEGG_db_dir -t 9 -ko_in user_ko_folder -x txt

### Prepare DB files, you need to have the following three files in your KEGG_db_dir:

1. Sequence file, only needed for "-seq_in" mode, DECOMPRESS and RENAME it to kegg_db_seq.fasta
   e.g. prokaryotes.pep.gz (https://www.kegg.jp/kegg/download/Readme/README.fasta)
2. seq2ko file, only needed for "-seq_in" mode, DECOMPRESS and RENAME it to kegg_db_seq2ko.txt
   e.g. prokaryotes.dat.gz (https://www.kegg.jp/kegg/download/Readme/README.fasta)
3. ko00001.keg
   https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir=

### seq2ko file format(tab separated)
aaa:Acav_4596	K05375
zpr:ZPR_0691	K21572

### How it works:
1. KEGG module uses Blast+/Diamond to get the best hits of query genes in the database with user defined e-value cutoff (default 0.001).
2. The TotalDepth of a KO is calculated by summing up the depth of all genes assigned to it.
3. The percentage of GeneNumber/TotalDepth of genes assigned to a KO is calculated by dividing them 
   by the total number/depth of genes with KO assignment (default) or by all genes in a genome ("-pct_by_all"). 

### Note!!!
1. If you run KEGG annotation for multiple files in a batch manner and want to have their depth info incorporated into the results, 
   you need to provide a folder containing individual depth files for each of your input sequence file.
   Name of the depth file needs to be exactly the same as its corresponding sequence file, except the extension which is ".depth".
2. Diamond requires quite a lot of memory for sequence comparison, especially for huge db file (e.g. KEGG db).
   Remember to request sufficient memory (e.g. 90 or 120gb) in your job script and specify a small number (e.g. -t 6) 
   of jobs executing in parallel. Otherwise, you may see some of your query genomes with no gene been annotated.

### Depth file format (one gene per line, tab separated)
gene_1	30
gene_2	10.58

+ To do:
1. level C stats: separate stats for Pathway, Brite and the rests
