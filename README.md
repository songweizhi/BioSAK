
## BioSAK (A Swiss-Army-Knife for Bioinformaticians)

[![pypi licence ](https://img.shields.io/pypi/l/BioSAK.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version ](https://img.shields.io/pypi/v/BioSAK.svg)](https://pypi.python.org/pypi/BioSAK) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4070001.svg)](https://doi.org/10.5281/zenodo.4070001)


Contact
---

Shan Zhang<sup>1,2</sup> and Weizhi Song<sup>1,2</sup>

<sup>1</sup> Center for Marine Science & Innovation, University of New South Wales, Sydney, Australia

<sup>2</sup> Department of Ocean Science, Hong Kong University of Science and Technology, Hong Kong

E-mail: shanzhang@ust.hk and songwz03@gmail.com


Installation
---

+ BioSAK has been tested on Linux/Mac, but **NOT** yet on Windows.


+ BioSAK is implemented in python3, you can

  + install it with `pip3 install BioSAK`
  
  + upgrade it with `pip3 install --upgrade BioSAK`


+ If you are an UNSW [Katana](https://research.unsw.edu.au/katana) user, [this](doc/katana.md) might be helpful.


 Getting help
---

+ You can get example commands for most of the modules by typing, for example, `BioSAK iTOL -h`.


+ Please refer to the [documentation](doc/Index.md) page for help (in preparation, very messy).


[//]: # (+ **中文版**帮助文件[在此]&#40;doc/Index_cn.md&#41;.)


+ A license statement is [here](LICENSE).


BioSAK modules
---

+ Type `BioSAK -h` to see a full list of modules


                ...::: BioSAK v1.105.0 :::...

      Genome databases
         get_GTDB_taxon_gnm      ->  Get id of genomes from specified GTDB taxons
         get_genome_GTDB         ->  Batch download GTDB genomes
         get_genome_NCBI         ->  Batch download GenBank genomes
         sampling_GTDB_gnms      ->  Select GTDB genomes
         subset_GTDB_meta        ->  Subset metadata of GTDB reference genomes
         metaAssembly            ->  Get metadata of NCBI assembly records
         metaBiosample           ->  Get metadata of NCBI biosample records
         statsTaxa               ->  stats GTDB taxa
  
      Metagenomics
         metabat2concoct         ->  convert MetaBAT depth to CONCOCT depth
         metabat2maxbin          ->  convert MetaBAT depth to MaxBin depth
         CheckM                  ->  Parse CheckM outputs
         Plot_MAG                ->  plot MAGs, (GC vs depth)
         magabund                ->  Calculate MAG abundance
         mean_MAG_cov            ->  Get mean MAG depth (by MetaBAT depth)
         RunGraphMB              ->  Prepare input files for GraphMB
         gc                      ->  Get GC content
         get_gnm_size            ->  Get the total length of genome(s)
         get_gene_depth          ->  Get gene depth by contig depth
         MeanMappingDepth        ->  Get mean mapping depth 
         get_MAG_reads_long      ->  Extract MAG-specific long reads for reassembling
         mmseqs                  ->  Classify metagenomic contigs with mmseqs
         parse_mmseqs_tsv        ->  Parse mmseqs tsv
         fastaai                 ->  A wrapper for FastAAI
         abd                     ->  get MAG abundance across metagenomes (Wenxiu Wang et al. 2024)

      Functional annotation
         KEGG                    ->  KEGG annotation
         koala                   ->  Separate the combined BlastKOALA or GhostKOALA output
         COG2020                 ->  COG annotation (v2020, by blastp/diamond)
         arCOG                   ->  COG annotation for archaea (version ar18)
         dbCAN                   ->  CAZy annotation with dbCAN
         Combine_KEGG_arCOG      ->  Combine KEGG and arCOG annotation results
         enrich                  ->  Functional enrichment analysis
         gapseq                  ->  Data matrix GapSeq predicted pathways
         stats_ko                ->  get stats for a list of provided KO
         stats_cog2020           ->  get stats for a list of provided COG (to be added)
         stats_arcog             ->  get stats for a list of provided arCOG 
  
      16S rRNA related
         top_16S_hits            ->  Classify 16S by top-blast-hits approach
         SILVA_for_BLCA          ->  Prepare BLCA-compatible SILVA SSU database
         GTDB_for_BLCA           ->  Prepare BLCA-compatible GTDB SSU database
         UNITE_for_BLCA          ->  Prepare BLCA-compatible UNITE SSU database
         BLCA_op_parser          ->  Make the BLCA outputs bit easier to read
         Tax4Fun2IndOTU          ->  Get functional profile for individual OTUs (to be added)
         UsearchNovogene         ->  Usearch for Novogene 16S amplicon sequencing results
          
      Sequence manipulator
         gbk2fna/gbk2faa/gbk2ffn ->  Format convertors
         ffn2faa/gfa2fa/get_rc   ->  Format convertors
         fq2fa                   ->  Convert fastq to fasta
         fa2id                   ->  Export sequence id
         slice_seq               ->  Get specified region of a sequence
         rename_seq              ->  Rename sequences in a file
         prefix_seq_by_file_name ->  prefix sequences by file name
         select_seq              ->  Select sequences by their ids
         split_fasta             ->  Split one fasta file into multiple files
         merge_seq               ->  Merge sequence files, remove duplicated ones if any
         cat_fa                  ->  Combine fasta files, prefix sequence id with file name
             
      Sam and Bam
         reads2bam               ->  Mapping and sorting
         sam2bam                 ->  Sam to BAM with samtools
         split_sam               ->  Split SAM/BAM file by reference
         bam2reads               ->  Extract reads (id) from sam file
         plot_sam_depth          ->  Plot SAM depth
      
      Dataframe and Statistics
         subset_df               ->  Subset dataframe
         merge_df                ->  Merge dataframes
         add_desc                ->  Add function description to input of the iTOL module
         transpose               ->  Transpose dataframe
         wilcox                  ->  Wilcoxon signed-rank test (non-parametric paired T-test)
         mannwhitneyu            ->  Mann-Whitney U rank test on two independent samples
      
      Others
         js_cmds                 ->  Put commands in job scripts
         js_hpc3                 ->  Put commands in job scripts (HKUST hpc3)
         hpc3                    ->  Submit jobs on HKUST hpc3
         srun                    ->  srun one-line commands on HKUST hpc3
         exe_cmds                ->  Execute commands with multiprocessing
         split_folder            ->  Split folder
         prefix_file             ->  Prefix file
         BestHit                 ->  Keep best blast hits (outfmt 6)
         VisGeneFlk              ->  Visualize gene flanking regions
         usearch_uc              ->  Parse Usearch uc file
         get_Pfam_hmms           ->  Get Pfam profiles by id
         Reads_simulator         ->  Simulate NGS reads
         SubsampleLongReads      ->  Subsample Long Reads
         rename_reads_Reago      ->  Rename paired reads for Reago
         cross_link_seqs         ->  Cross link matched regions between two sequences
         submitHPC               ->  A wrapper for submitHPC.sh
         KeepRemovingTmp         ->  Keep removing old files in a folder
         ribbon                  ->  Make a ribbon diagram
         compare_sets            ->  compare_sets
         sankey                  ->  get sankey plot
         sra                     ->  Download reads with sratoolkit
         vis_color_scheme        ->  Visualize color scheme
         trim                    ->  a wrapper for trimmomatic
         FasterqDump             ->  a wrapper for fasterq-dump
  

+ A changelog is [here](BioSAK/VERSION).

