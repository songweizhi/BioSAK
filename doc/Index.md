
Help information for BioSAK modules
---


### Functional annotation

+ __[KEGG](KEGG.md)__                   			- KEGG annotation
+ __[COG2020](COG2020.md)__                			- COG annotation (v2020, by blastp/diamond)
+ __[arCOG](arCOG.md)__                  			- COG annotation for archaea (version ar18)
+ __[dbCAN](dbCAN.md)__                  			- CAZy annotation with dbCAN 
+ __[CheckM](CheckM.md)__                 			- Parse CheckM outputs

### Metagenomics

+ __[magabund](magabund.md)__               		- Calculate MAG abundance
+ __[mean_MAG_cov](mean_MAG_cov.md)__           	- Get mean MAG depth (based on MetaBAT produced depth)
+ __[RunGraphMB](RunGraphMB.md)__             		- Prepare input files for GraphMB
+ __[get_MAG_reads_long](get_MAG_reads_long.md)__   - Extract MAG-specific long reads for reassembling
+ __[get_gnm_size](get_gnm_size.md)__           	- Get the total length of genome(s)
+ __[get_gene_depth](get_gene_depth.md)__         	- Get gene depth by contig depth
+ __[MeanMappingDepth](MeanMappingDepth.md)__       - Get mean mapping depth 
+ __[Plot_MAG](Plot_MAG.md)__               		- plot MAGs, (GC vs depth)

### Genome databases

+ __[get_GTDB_taxon_gnm](get_GTDB_taxon_gnm.md)__ 	- Get id of genomes from specified GTDB taxons
+ __[get_genome_GTDB](get_genome_GTDB.md)__        	- Batch download GTDB genomes
+ __[get_genome_NCBI](get_genome_NCBI.md)__        	- Batch download GenBank genomes
+ __[sampling_GTDB_gnms](sampling_GTDB_gnms.md)__   - Select GTDB genomes from a specified taxon at specified sampling rank
+ __[subset_GTDB_meta](subset_GTDB_meta.md)__       - Subset metadata of GTDB reference genomes

### 16S rRNA related

+ __[top_16S_hits](top_16S_hits.md)__           	- Classify 16S by top-blast-hits approach
+ __[SILVA_for_BLCA](SILVA_for_BLCA.md)__        	- Prepare BLCA-compatible SILVA SSU database
+ __[GTDB_for_BLCA](GTDB_for_BLCA.md)__          	- Prepare BLCA-compatible GTDB SSU database
+ __[BLCA_op_parser](BLCA_op_parser.md)__         	- Make the BLCA outputs bit easier to read
+ __[Tax4Fun2IndOTU](Tax4Fun2IndOTU.md)__         	- Get functional profile of individual OTUs (to be added)

### Sequence/Alignment manipulator

+ __[gbk2fa/gbk2faa/gbk2ffn](.md)__ 				- Format convertors
+ __[ffn2faa/gfa2fa/get_rc](.md)__  				- Format convertors
+ __[fq2fa](fq2fa.md)__                  			- Convert fastq to fasta
+ __[fa2id](fa2id.md)__                  			- Export sequence id
+ __[slice_seq](slice_seq.md)__              		- Get specified region of a sequence
+ __[rename_seq](rename_seq.md)__            	 	- Rename sequences in a file
+ __[select_seq](select_seq.md)__             		- Select sequences by their ids
+ __[split_fasta](split_fasta.md)__            		- Split one fasta file into multiple files
+ __[merge_seq](merge_seq.md)__              		- Merge sequence files, remove duplicated ones if any
+ __[cat_fa](cat_fa.md)__                 			- Combine sequence files, prefix sequence id with file name

### Sam/Bam

+ __[plot_sam_depth](plot_sam_depth.md)__         	- Plot SAM depth
+ __[reads2bam](reads2bam.md)__              		- Mapping and sorting
+ __[sam2bam](sam2bam.md)__                			- Sam to BAM with samtools
+ __[split_sam](split_sam.md)__              		- Split SAM/BAM file by reference
+ __[bam2reads](bam2reads.md)__              		- Extract reads (id) from sam file

### Others

+ __[iTOL](iTOL.md)__                   			- Prepare iTOL-compatible files for tree visualization
+ __[js_cmds](js_cmds.md)__                			- Commands to job scripts
+ __[exe_cmds](exe_cmds.md)__               		- Execute commands with multiprocessing (one cmd per line)
+ __[BestHit](BestHit.md)__                			- Keep best blast hits (outfmt 6)
+ __[VisGeneFlk](VisGeneFlk.md)__             		- Visualize gene flanking regions
+ __[split_folder](split_folder.md)__           	- Split folder
+ __[usearch_uc](usearch_uc.md)__             		- Parse Usearch uc file
+ __[get_Pfam_hmms](get_Pfam_hmms.md)__          	- Get Pfam profiles by id
+ __[Reads_simulator](Reads_simulator.md)__        	- Simulate NGS reads
+ __[SubsampleLongReads](SubsampleLongReads.md)__   - Subsample Long Reads
+ __[rename_reads_Reago](rename_reads_Reago.md)__   - Rename paired reads for Reago
+ __[cross_link_seqs](cross_link_seqs.md)__        	- Cross link matched regions between two sequences (not working well!)

### Note

Phylogenetic tree-related modules have been moved to TreeSAK, please go to its github page for details.
