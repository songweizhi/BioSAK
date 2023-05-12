
# Usage:
Before you start, you need to download https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt
and provide it to get_genome_NCBI with '-db'.

# Download genomes in the gnms_to_download.txt
    BioSAK get_genome_NCBI -db prokaryotes.txt -out gnm_dir -id gnms_to_download.txt -t 8 -fna
    BioSAK get_genome_NCBI -db prokaryotes.txt -out gnm_dir -id gnms_to_download.txt -t 8 -fna -faa -name

# Format of gnms_to_download.txt (genome ID can be found in the "Assembly Accession" column of the prokaryotes.txt)
GCA_009840555.1
GCA_009840575.1

# You can get the id of genomes from a specific taxon using this link, IDs are in the "Assembly" column.
https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/refseq_category:reference
