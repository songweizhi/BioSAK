### Dependencies
  + blast+
  + diamond

### Annotate protein sequences

      BioSAK COG2020 -m P -t 6 -db_dir path/to/your/COG_db_dir -i genes.faa 
      BioSAK COG2020 -m P -t 6 -db_dir path/to/your/COG_db_dir -i faa_files -x faa -depth depth_files

### Annotate DNA sequences (ORFs)

      BioSAK COG2020 -m N -t 6 -db_dir path/to/your/COG_db_dir -i genes.ffn -depth gene.depth
      BioSAK COG2020 -m N -t 6 -db_dir path/to/your/COG_db_dir -i ffn_files -x ffn

### Prepare DB files (version 2020):

      cd path/to/your/COG_db_dir
      wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.fa.gz
      wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.cog.csv
      wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab
      wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab
      wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/Readme.2020-11-25.txt
      gunzip cog-20.fa.gz
      module load blast
      makeblastdb -in cog-20.fa -dbtype prot -parse_seqids -logfile cog-20.fa.log
      module load diamond
      diamond makedb --in cog-20.fa --db cog-20.fa.dmnd --quiet

### How it works:
1. COG2020 module uses Blast+/Diamond to get the best hits of query genes in the database 
   with users defined e-value cutoff (default 0.001).

1. The TotalDepth of a COG id/category is obtained by summing up the depth of all genes assigned to it.

1. The percentage of GeneNumber/TotalDepth of genes assigned to a COG is calculated by dividing them 
   by the total number/depth of genes with COG assignment (default) or all query genes in a file (if "-pct_by_all" specified). 

### Note!!!

If you run COG2020 for multiple files in a batch manner and want to have their depth info incorporated into the results, 
you need to provide a folder containing individual depth file for each of your input sequence file.
Name of the depth file needs to be exactly the same as its corresponding sequence file, except the extension is ".depth".

#### Depth file format (one gene per line, tab separated)

    gene_1	30
    gene_2	10.58
