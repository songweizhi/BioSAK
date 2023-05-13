    module load blast+/2.9.0
    BioSAK top_16S_hits -p 16S_vs_GTDB -q 16S_sequences.fa -r GTDB_ssu_all_r95.fna -t 6
    BioSAK top_16S_hits -p 16S_vs_SILVA -q 16S_sequences.fa -r SILVA_138_SSURef_NR99_tax_silva.fasta -t 6 -top 5

# prepare GTDB SSU database file
    wget https://data.gtdb.ecogenomic.org/releases/release95/95.0/genomic_files_all/ssu_all_r95.tar.gz
    tar xvzf ssu_all_r95.tar.gz
    makeblastdb -in GTDB_ssu_all_r95.fna -dbtype nucl -parse_seqids -logfile /dev/null

# prepare SILVA SSU database file
    wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva.fasta.gz
    gunzip SILVA_138_SSURef_NR99_tax_silva.fasta.gz
    makeblastdb -in SILVA_138_SSURef_NR99_tax_silva.fasta -dbtype nucl -parse_seqids -logfile /dev/null
