
### Description 



### Example commands

+ Extract reads for a single MAG

      BioSAK get_MAG_reads_long -r LongReads.fastq -s assembly.sam -m MAG_1.fa -o MAG_1_reads
    
+ Extract reads for multiple MAGs

      BioSAK get_MAG_reads_long -r LongReads.fastq -s assembly.sam -m MAG_dir -x fa -o extracted_reads -l 0

+ Extract reads for multiple MAGs, ignore reads shorter than 5000 bp

      BioSAK get_MAG_reads_long -r LongReads.fastq -s assembly.sam -m MAG_dir -x fa -o extracted_reads -l 5000


### Notes

+ The input sam file is obtained by mapping quality-filtered reads to all (not a subset of) contigs derived 
   from these reads (e.g., flye produced assembly.fasta).


+ The input MAGs and Sam file need to be generated from the same set of contigs (without renaming).


+ When a read is assigned to contigs from multiple MAGs, it will be added to the extracted reads for each of those MAGs.


+ If you use minimap2 for reads mapping, you might want to specify "--secondary=no" to minimise low-quality secondary alignments:

      minimap2 -d assembly.fasta.mmi assembly.fasta
      minimap2 --secondary=no -t 6 -a assembly.fasta.mmi long_reads.fastq > assembly.sam
