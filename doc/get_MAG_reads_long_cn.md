BioSAK get_MAG_reads_long -r LongReads.fastq -s assembly.sam -m MAG_1.fa -o MAG_1_reads
BioSAK get_MAG_reads_long -r LongReads.fastq -s assembly.sam -m MAG_dir -x fa -o extracted_reads -l 0
BioSAK get_MAG_reads_long -r LongReads.fastq -s assembly.sam -m MAG_dir -x fa -o extracted_reads -l 5000

# Note
1. The input sam file is obtained by mapping quality-filtered reads to all (not a subset of) contigs derived 
   from these reads (e.g., flye produced assembly.fasta).
2. The input MAGs and Sam file need to be generated from the same set of contigs (without renaming).
3. If one read was mapped to multiple contigs from different MAGs, this read will be exported to the extracted 
   reads of all the above MAGs.
4. If you use minimap2 for reads mapping, you might want to specify "--secondary=no" to minimise low-quality 
   secondary alignments:
   minimap2 -d assembly.fasta.mmi assembly.fasta
   minimap2 --secondary=no -t 6 -a assembly.fasta.mmi long_reads.fastq > assembly.sam
