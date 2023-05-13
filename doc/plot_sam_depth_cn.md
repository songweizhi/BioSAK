
BioSAK plot_sam_depth -r ref.fa -d 12D9.bam.depth
BioSAK plot_sam_depth -r ref.fa -d 12D9.bam.depth -i contig_01 -s 100 -e 900 -k 300 -l 500,550,600

# get depth file 
module load samtools/1.15
samtools depth 12D9.bam > 12D9.bam.depth
