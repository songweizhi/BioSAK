
    module load samtools/1.15
    BioSAK split_sam -p single_ctg -i input.bam -r contig_1 -t 12
    BioSAK split_sam -p multi_ctgs -i input.bam -r ctgs.txt -t 12

# Output files:
    prefix.sam, prefix_sorted.bam and prefix_sorted.bam.bai

# ctgs.txt file format: one id per line, ">" excluded.
