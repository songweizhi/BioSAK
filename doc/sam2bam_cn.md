
    module load samtools
    BioSAK sam2bam -sam input.sam -t 12

# This is a wrapper for the following three steps:
    samtools view -bS input.sam -o input.bam
    samtools sort input.bam -o input_sorted.bam
    samtools index input_sorted.bam
