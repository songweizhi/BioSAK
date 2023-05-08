
## BioSAK (A Swiss-Army-Knife for Bioinformaticians)

[![pypi licence ](https://img.shields.io/pypi/l/BioSAK.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version ](https://img.shields.io/pypi/v/BioSAK.svg)](https://pypi.python.org/pypi/BioSAK) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4070001.svg)](https://doi.org/10.5281/zenodo.4070001)


Contact
---

Shan Zhang<sup>1</sup> and Weizhi Song<sup>1,2</sup>

<sup>1</sup> Center for Marine Science & Innovation, University of New South Wales, Sydney, Australia

<sup>2</sup> School of Life Sciences, The Chinese University of Hong Kong, Hong Kong, China

E-mail: zzfanyi@gmail.com and songwz03@gmail.com

    
Installation
---

1. BioSAK has been tested on Linux/Mac, but **NOT** on Windows.

1. BioSAK is implemented in python3, you can install it with pip3:

       # for the first time installation
       pip3 install BioSAK
      
       # for updating
       pip3 install --upgrade BioSAK

1. If you are an UNSW Katana user, please install as described [here](doc/README_installation_UNSW.md):


BioSAK modules
---
   
                 ...::: BioSAK v1.63.5 :::...

    Functional annotation
       KEGG                   ->  KEGG annotation
       COG2020                ->  COG annotation (v2020, by blastp/diamond)
       arCOG                  ->  COG annotation for archaea (version ar18)
       dbCAN                  ->  CAZy annotation with dbCAN

    Metagenomics
       CheckM                 ->  Parse CheckM outputs
       magabund               ->  Calculate MAG abundance
       mean_MAG_cov           ->  Get mean MAG depth (based on MetaBAT produced depth)
       RunGraphMB             ->  Prepare input files for GraphMB
       get_MAG_reads_long     ->  Extract MAG-specific long reads for reassembling
       get_gnm_size           ->  Get the total length of genome(s)
       get_gene_depth         ->  Get gene depth by contig depth
       MeanMappingDepth       ->  Get mean mapping depth 
       Plot_MAG               ->  plot MAGs, (GC vs depth)

    Genome databases
       get_GTDB_taxon_gnm     ->  Get id of genomes from specified GTDB taxons
       get_genome_GTDB        ->  Batch download GTDB genomes
       get_genome_NCBI        ->  Batch download GenBank genomes
       sampling_GTDB_gnms     ->  Select GTDB genomes from a specified taxon at specified sampling rank
       subset_GTDB_meta       ->  Subset metadata of GTDB reference genomes

    16S rRNA related
       top_16S_hits           ->  Classify 16S by top-blast-hits approach
       SILVA_for_BLCA         ->  Prepare BLCA-compatible SILVA SSU database
       GTDB_for_BLCA          ->  Prepare BLCA-compatible GTDB SSU database
       BLCA_op_parser         ->  Make the BLCA outputs bit easier to read
       Tax4Fun2IndOTU         ->  Get functional profile of individual OTUs (to be added)
    
    Sequence manipulator
       gbk2fa/gbk2faa/gbk2ffn ->  Format convertors
       ffn2faa/gfa2fa/get_rc  ->  Format convertors
       fq2fa                  ->  Convert fastq to fasta
       fa2id                  ->  Export sequence id
       slice_seq              ->  Get specified region of a sequence
       rename_seq             ->  Rename sequences in a file
       select_seq             ->  Select sequences by their ids
       split_fasta            ->  Split one fasta file into multiple files
       merge_seq              ->  Merge sequence files, remove duplicated ones if any
       cat_fa                 ->  Combine sequence files, prefix sequence id with file name
    
    Sam/Bam
       plot_sam_depth         ->  Plot SAM depth
       reads2bam              ->  Mapping and sorting
       sam2bam                ->  Sam to BAM with samtools
       split_sam              ->  Split SAM/BAM file by reference
       bam2reads              ->  Extract reads (id) from sam file
    
    Others
       iTOL                   ->  Prepare iTOL-compatible files for tree visualization
       js_cmds                ->  Commands to job scripts
       exe_cmds               ->  Execute commands with multiprocessing (one cmd per line)
       BestHit                ->  Keep best blast hits (outfmt 6)
       VisGeneFlk             ->  Visualize gene flanking regions
       split_folder           ->  Split folder
       usearch_uc             ->  Parse Usearch uc file
       get_Pfam_hmms          ->  Get Pfam profiles by id
       Reads_simulator        ->  Simulate NGS reads
       SubsampleLongReads     ->  Subsample Long Reads
       rename_reads_Reago     ->  Rename paired reads for Reago
       cross_link_seqs        ->  Cross link matched regions between two sequences (not working well!)
       
    Phylogenetic tree-related modules have been moved to TreeSAK, please go to its github page for details.


Get help
---
1. Example usage for most of BioSAK modules can be found from their own help information, e.g.
        
       BioSAK COG2020 -h
       BioSAK iTOL -h 

2. Short tutorials for a few modules:

    https://github.com/songweizhi/BioSAK/tree/master/BioSAK_tutorial
