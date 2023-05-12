module load samtools/1.15

    BioSAK bam2reads -b sorted.bam -o bam2reads_wd_ctg1 -r ctg_1
    BioSAK bam2reads -b sorted.bam -o bam2reads_wd_ctg2 -r ctg_2:200-5000
    BioSAK bam2reads -b sorted.bam -o bam2reads_wd_ctg3 -r ctg_3:200-
    BioSAK bam2reads -b sorted.bam -o bam2reads_wd_ctg4 -r ctg_4:-500 -s reads.fastq
