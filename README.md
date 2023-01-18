
## BioSAK (A Swiss-Army-Knife for Bioinformaticians)

[![pypi licence ](https://img.shields.io/pypi/l/BioSAK.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version ](https://img.shields.io/pypi/v/BioSAK.svg)](https://pypi.python.org/pypi/BioSAK) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4070001.svg)](https://doi.org/10.5281/zenodo.4070001)


Contact
---

+ **Weizhi Song**, Postdoctoral Researcher
+ Center for Marine Science & Innovation, University of New South Wales, Sydney, Australia
+ E-mail: songwz03@gmail.com

    
Installation
---

1. BioSAK has been tested on Linux/Mac, but **NOT** on Windows.

1. BioSAK is implemented in python3, you can install it with pip3:

       # for the first time installation
       pip3 install BioSAK
      
       # for updating
       pip3 install --upgrade BioSAK
      
1. If you are an UNSW Katana user, install with:

       ######### Install BioSAK with Python's virtual environment ########

       module load python/3.7.3
       mkdir ~/mypython3env_BioSAK
       python3 -m venv --system-site-packages ~/mypython3env_BioSAK
       source ~/mypython3env_BioSAK/bin/activate
       pip3 install BioSAK

       ####################### For later running #########################

       module load python/3.7.3
       source ~/mypython3env_BioSAK/bin/activate
       BioSAK -h
              
       # Type "deactivate" in your terminal to leave Python's virtual environment.


BioSAK modules
---
   
                 ...::: BioSAK v1.60.2 :::...

    Annotation modules
       CheckM                 ->  Parse CheckM outputs
       KEGG                   ->  KEGG annotation
       COG2020                ->  COG annotation (v2020, by blastp/diamond)
       dbCAN                  ->  CAZy annotation with dbCAN

    Metagenomics
       magabund               ->  Calculate MAG abundance
       mean_MAG_cov           ->  Get mean MAG depth (based on MetaBAT produced depth)
       RunGraphMB             ->  Prepare input files for GraphMB
       get_MAG_reads_long     ->  Extract MAG-specific long reads for reassembling

    16S rRNA related
       top_16S_hits           ->  Classify 16S by top-blast-hits approach
       SILVA_for_BLCA         ->  Prepare BLCA-compatible SILVA SSU database
       GTDB_for_BLCA          ->  Prepare BLCA-compatible GTDB SSU database
       BLCA_op_parser         ->  Make the BLCA outputs bit easier to read
    
    Sequence manipulator
       gbk2fa/gbk2faa         ->  Format convertors
       gbk2ffn/ffn2faa        ->  Format convertors
       get_rc                 ->  Get reverse complement sequence
       slice_seq              ->  Get specified region of a sequence
       rename_seq             ->  Rename sequences in a file
       select_seq             ->  Select sequences by their ids
       split_fasta            ->  Split one fasta file into multiple files
       merge_seq              ->  Merge sequence files, remove duplicated ones if any
       cat_fa                 ->  Combine sequence files, prefix sequence id with file name
       get_gene_depth         ->  Get gene depth by contig depth
       get_gnm_size           ->  Get the total length of genome(s)
       fa2id                  ->  Export sequence id

    Multiple Sequence Alignment (MSA)
       convert_align_format   ->  Convert alignment format
       OneLineAln             ->  One-line fasta format alignments
       SubsetAlnCols          ->  Subset MSA by column    
    
    Sam and Bam
       plot_sam_depth         ->  Plot SAM depth
       reads2bam              ->  Mapping and sorting
       sam2bam                ->  Sam to BAM with samtools
       split_sam              ->  Split SAM/BAM file by reference
       bam2reads              ->  Extract reads (id) from sam file
       MeanMappingDepth       ->  Get mean mapping depth 

    Tree-related
       GTDB_tree_r207         ->  Infer GTDB (r207) archaeal/bacterial tree
       subset_tree            ->  Subset tree
       compare_trees          ->  Compare trees with Mantel test
       rename_leaves          ->  Rename tree leaves
       FLN                    ->  Format leaf names (e.g. remove spaces in names)
       iTOL                   ->  Prepare iTOL-compatible files for tree visualization
    
    Genome databases
       get_GTDB_taxon_gnm     ->  Get id of genomes from specified GTDB taxons
       get_genome_GTDB        ->  Batch download GTDB genomes
       get_genome_NCBI        ->  Batch download GenBank genomes
       sampling_GTDB_gnms     ->  Select GTDB genomes from a specified taxon at specified sampling rank
       subset_GTDB_meta       ->  Subset metadata of GTDB reference genomes
        
    Visualization
       Plot_MAG               ->  plot MAGs, (GC vs depth)
       VisGeneFlk             ->  Visualize gene flanking regions
       cross_link_seqs        ->  Cross link matched regions between two sequences
          
    Other modules
       BestHit                ->  Keep best blast hits (outfmt 6)
       js_cmds                ->  Commands to job scripts
       exe_cmds               ->  Execute commands in a file with multiprocessing (one cmd per line)
       split_folder           ->  Split folder
       usearch_uc             ->  Parse Usearch uc file
       get_Pfam_hmms          ->  Get Pfam profiles by id
       Reads_simulator        ->  Simulate NGS reads
       rename_reads_Reago     ->  Rename paired reads for Reago

Get help
---

1. Example usage for most of BioSAK modules can be found from their own help information, e.g.
        
       BioSAK COG2020 -h
       BioSAK iTOL -h 

2. Short tutorials for a few modules:

    https://github.com/songweizhi/BioSAK/tree/master/BioSAK_tutorial
