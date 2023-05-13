
# Get gene depth by contig depth, together with gbk or gff (prefered) file
BioSAK get_gene_depth -gff test.gff -ctg_depth contig_depth.txt
BioSAK get_gene_depth -gff test.gff -ctg_depth contig_depth.txt -skip_header
BioSAK get_gene_depth -gbk test.gbk -ctg_depth contig_depth.txt -skip_header

# Contig depth file format (one contig per line, tab separated)
contig_1    30.16
contig_2    26
contig_3    0.726