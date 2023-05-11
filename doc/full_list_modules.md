

    % BioSAK -h

                        ...::: BioSAK :::...

    Functional annotation
       KEGG                    ->  KEGG annotation
       COG2020                 ->  COG annotation (v2020, by blastp/diamond)
       arCOG                   ->  COG annotation for archaea (version ar18)
       dbCAN                   ->  CAZy annotation with dbCAN
       CheckM                  ->  Parse CheckM outputs

    Metagenomics
       magabund                ->  Calculate MAG abundance
       mean_MAG_cov            ->  Get mean MAG depth (based on MetaBAT produced depth)
       RunGraphMB              ->  Prepare input files for GraphMB
       get_MAG_reads_long      ->  Extract MAG-specific long reads for reassembling
       get_gnm_size            ->  Get the total length of genome(s)
       get_gene_depth          ->  Get gene depth by contig depth
       MeanMappingDepth        ->  Get mean mapping depth 
       Plot_MAG                ->  plot MAGs, (GC vs depth)

    Genome databases
       get_GTDB_taxon_gnm      ->  Get id of genomes from specified GTDB taxons
       get_genome_GTDB         ->  Batch download GTDB genomes
       get_genome_NCBI         ->  Batch download GenBank genomes
       sampling_GTDB_gnms      ->  Select GTDB genomes from a specified taxon at specified sampling rank
       subset_GTDB_meta        ->  Subset metadata of GTDB reference genomes

    16S rRNA related
       top_16S_hits            ->  Classify 16S by top-blast-hits approach
       SILVA_for_BLCA          ->  Prepare BLCA-compatible SILVA SSU database
       GTDB_for_BLCA           ->  Prepare BLCA-compatible GTDB SSU database
       BLCA_op_parser          ->  Make the BLCA outputs bit easier to read
       Tax4Fun2IndOTU          ->  Get functional profile of individual OTUs (to be added)
    
    Sequence/Alignment manipulator
       gbk2fna/gbk2faa/gbk2ffn ->  Format convertors
       ffn2faa/gfa2fa/get_rc   ->  Format convertors
       fq2fa                   ->  Convert fastq to fasta
       fa2id                   ->  Export sequence id
       slice_seq               ->  Get specified region of a sequence
       rename_seq              ->  Rename sequences in a file
       select_seq              ->  Select sequences by their ids
       split_fasta             ->  Split one fasta file into multiple files
       merge_seq               ->  Merge sequence files, remove duplicated ones if any
       cat_fa                  ->  Combine sequence files, prefix sequence id with file name
    
    Sam/Bam
       plot_sam_depth          ->  Plot SAM depth
       reads2bam               ->  Mapping and sorting
       sam2bam                 ->  Sam to BAM with samtools
       split_sam               ->  Split SAM/BAM file by reference
       bam2reads               ->  Extract reads (id) from sam file
    
    Others
       iTOL                    ->  Prepare iTOL-compatible files for tree visualization
       js_cmds                 ->  Commands to job scripts
       exe_cmds                ->  Execute commands with multiprocessing (one cmd per line)
       BestHit                 ->  Keep best blast hits (outfmt 6)
       VisGeneFlk              ->  Visualize gene flanking regions
       split_folder            ->  Split folder
       usearch_uc              ->  Parse Usearch uc file
       get_Pfam_hmms           ->  Get Pfam profiles by id
       Reads_simulator         ->  Simulate NGS reads
       SubsampleLongReads      ->  Subsample Long Reads
       rename_reads_Reago      ->  Rename paired reads for Reago
       cross_link_seqs         ->  Cross link matched regions between two sequences (not working well!)
       
    Phylogenetic tree-related modules have been moved to TreeSAK, please go to its github page for details.
