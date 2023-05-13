
帮助
---


### Functional annotation

+ __[KEGG](KEGG_cn.md)__                   			    - 使用 KEGG 数据库进行功能注释
+ __[COG2020](COG2020_cn.md)__                			- 使用 COG 数据库进行功能注释 (v2020)
+ __[arCOG](arCOG_cn.md)__                  			- 使用 arCOG 数据库进行功能注释 (ar18)
+ __[dbCAN](dbCAN_cn.md)__                  			- 使用 CAZy 数据库进行功能注释 (需要调用dbCAN)

### Metagenomics

+ __[CheckM](CheckM_cn_cn.md)__                     	- 解析 CheckM 输出结果
+ __[magabund](magabund_cn.md)__               		    - Calculate MAG abundance
+ __[mean_MAG_cov](mean_MAG_cov_cn.md)__           	    - Get mean MAG depth (based on MetaBAT produced depth)
+ __[RunGraphMB](RunGraphMB_cn.md)__             		- Prepare input files for GraphMB
+ __[get_MAG_reads_long](get_MAG_reads_long_cn.md)__    - Extract MAG-specific long reads for reassembling
+ __[get_gnm_size](get_gnm_size_cn.md)__           	    - Get the total length of genome(s)
+ __[get_gene_depth](get_gene_depth_cn.md)__         	- Get gene depth by contig depth
+ __[MeanMappingDepth](MeanMappingDepth_cn.md)__        - Get mean mapping depth 
+ __[Plot_MAG](Plot_MAG_cn.md)__               		    - plot MAGs, (GC vs depth)

### Genome databases

+ __[get_GTDB_taxon_gnm](get_GTDB_taxon_gnm_cn.md)__ 	- Get id of genomes from specified GTDB taxons
+ __[get_genome_GTDB](get_genome_GTDB_cn.md)__        	- Batch download GTDB genomes
+ __[get_genome_NCBI](get_genome_NCBI_cn.md)__        	- Batch download GenBank genomes
+ __[sampling_GTDB_gnms](sampling_GTDB_gnms_cn.md)__    - Select GTDB genomes from a specified taxon at specified sampling rank
+ __[subset_GTDB_meta](subset_GTDB_meta_cn.md)__        - Subset metadata of GTDB reference genomes

### 16S rRNA related

+ __[top_16S_hits](top_16S_hits_cn.md)__           	    - Classify 16S by top-blast-hits approach
+ __[SILVA_for_BLCA](SILVA_for_BLCA_cn.md)__        	- Prepare BLCA-compatible SILVA SSU database
+ __[GTDB_for_BLCA](GTDB_for_BLCA_cn.md)__          	- Prepare BLCA-compatible GTDB SSU database
+ __[BLCA_op_parser](BLCA_op_parser_cn.md)__         	- Make the BLCA outputs bit easier to read
+ __[Tax4Fun2IndOTU](Tax4Fun2IndOTU_cn.md)__         	- Get functional profile of individual OTUs (to be added)

### Sequence/Alignment manipulator

+ __[gbk2fa/gbk2faa/gbk2ffn](_cn.md)__ 				    - Format convertors
+ __[ffn2faa/gfa2fa/get_rc](_cn.md)__  				    - Format convertors
+ __[fq2fa](fq2fa_cn.md)__                  			- Convert fastq to fasta
+ __[fa2id](fa2id_cn.md)__                  			- Export sequence id
+ __[slice_seq](slice_seq_cn.md)__              		- Get specified region of a sequence
+ __[rename_seq](rename_seq_cn.md)__            	 	- Rename sequences in a file
+ __[select_seq](select_seq_cn.md)__             		- Select sequences by their ids
+ __[split_fasta](split_fasta_cn.md)__            		- Split one fasta file into multiple files
+ __[merge_seq](merge_seq_cn.md)__              		- Merge sequence files, remove duplicated ones if any
+ __[cat_fa](cat_fa_cn.md)__                 			- Combine sequence files, prefix sequence id with file name

### Sam/Bam

+ __[plot_sam_depth](plot_sam_depth_cn.md)__         	- Plot SAM depth
+ __[reads2bam](reads2bam_cn.md)__              		- Mapping and sorting
+ __[sam2bam](sam2bam_cn.md)__                			- Sam to BAM with samtools
+ __[split_sam](split_sam_cn.md)__              		- Split SAM/BAM file by reference
+ __[bam2reads](bam2reads_cn.md)__              		- Extract reads (id) from sam file

### Others

+ __[iTOL](iTOL_cn.md)__                   			    - Prepare iTOL-compatible files for tree visualization
+ __[js_cmds](js_cmds_cn.md)__                			- Commands to job scripts
+ __[exe_cmds](exe_cmds_cn.md)__               		    - Execute commands with multiprocessing (one cmd per line)
+ __[BestHit](BestHit_cn.md)__                			- Keep best blast hits (outfmt 6)
+ __[VisGeneFlk](VisGeneFlk_cn.md)__             		- Visualize gene flanking regions
+ __[split_folder](split_folder_cn.md)__           	    - Split folder
+ __[usearch_uc](usearch_uc_cn.md)__             		- Parse Usearch uc file
+ __[get_Pfam_hmms](get_Pfam_hmms_cn.md)__          	- Get Pfam profiles by id
+ __[Reads_simulator](Reads_simulator_cn.md)__        	- Simulate NGS reads
+ __[SubsampleLongReads](SubsampleLongReads_cn.md)__    - Subsample Long Reads
+ __[rename_reads_Reago](rename_reads_Reago_cn.md)__    - Rename paired reads for Reago
+ __[cross_link_seqs](cross_link_seqs_cn.md)__        	- Cross link matched regions between two sequences (not working well!)

### Note

Phylogenetic tree-related modules have been moved to TreeSAK, please go to its github page for details.
