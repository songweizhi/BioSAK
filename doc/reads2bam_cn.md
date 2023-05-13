
+ Dependencies: bowtie2 and samtools

      BioSAK reads2bam -ref ref.fa -r1 R1.fq -r2 R1.fq -fq -index -t 12
      BioSAK reads2bam -ref ref.fa -r1 R1.fa -r2 R1.fa -up unpaired.fa -index -t 12
      BioSAK reads2bam -ref ref.fa -up unpaired.fa -index -local -t 12 -tmp
      BioSAK reads2bam -ref ref.fa -up unpaired_R1.fa,unpaired_R2.fa -index -t 12 -tmp
