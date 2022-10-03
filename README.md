
## BioSAK (A Swiss-Army-Knife for Biologists)

[![pypi licence ](https://img.shields.io/pypi/l/BioSAK.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version ](https://img.shields.io/pypi/v/BioSAK.svg)](https://pypi.python.org/pypi/BioSAK) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4070001.svg)](https://doi.org/10.5281/zenodo.4070001)


Contact
---

+ **Weizhi Song**, Postdoctoral Researcher
+ Center for Marine Science & Innovation, University of New South Wales, Sydney, Australia
+ E-mail: songwz03@gmail.com


Dependencies
---

Module dependent, see module-specific help info for details. E.g.

    BioSAK COG2020 -h
    BioSAK dbCAN -h
    
    
Installation
---

1. BioSAK has been tested on Linux/Mac, but **NOT** on Windows.

1. BioSAK is implemented in python3, you can install it with pip3:

       # for the first time installation
       pip3 install BioSAK
      
       # for updating
       pip3 install --upgrade BioSAK
      
1. For UNSW Katana users

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
   
                 ...::: BioSAK v1.48.0 :::...

    Annotation modules
       CheckM                ->  Parse CheckM outputs
       KEGG                  ->  KEGG annotation
       COG2020               ->  COG annotation (v2020, by blastp/diamond)
       dbCAN                 ->  CAZy annotation with dbCAN

    16S rRNA related
       top_16S_hits          ->  Classify 16S by top-blast-hits approach
       SILVA_for_BLCA        ->  Prepare BLCA-compatible SILVA SSU database
       GTDB_for_BLCA         ->  Prepare BLCA-compatible GTDB SSU database
       BLCA_op_parser        ->  Make the BLCA outputs bit easier to read
    
    Sequence manipulator
       gbk2fa                ->  gbk to fasta
       gfa2fa                ->  gfa to fasta (to be added)
       gbk2ffn               ->  gbk to ffn
       gbk2faa               ->  gbk to faa
       ffn2faa               ->  ffn to faa
       get_rc                ->  Get reverse complement sequence
       slice_seq             ->  Get specified region of a sequence
       rename_seq            ->  Rename sequences in a file
       select_seq            ->  Select sequences by their ids
       split_fasta           ->  Split a fasta file into multiple files
       fa2id                 ->  export sequence id
       get_gene_depth        ->  Get gene depth by contig depth
       convert_align_format  ->  Convert alignment format
       OneLineAln            ->  One-line fasta format alignments
       SubsetAlnCols         ->  Subset MSA by column
       rename_reads_Reago    ->  Rename paired reads for Reago
       MeanMappingDepth      ->  Get mean mapping depth 
       mean_MAG_cov          ->  get mean MAG depth, weighted by contig length 
       get_gnm_size          ->  Get the total length of genome(s)
       get_bin_abundance     ->  Get bin abundance
       merge_seq             ->  Merge sequence files, remove duplicated ones if any

    Sam and Bam
       plot_sam_depth        ->  Plot SAM depth
       reads2bam             ->  Mapping and sorting
       sam2bam               ->  Sam to BAM with samtools
       split_sam             ->  Split SAM/BAM file by reference
       bam2reads             ->  Extract reads (id) from sam file

    Tree related
       GTDB_tree_r207        ->  Infer GTDB (r207) archaeal/bacterial tree
       subset_tree           ->  Subset tree
       compare_trees         ->  Compare trees with Mantel test
       rename_leaves         ->  Rename tree leaves
       FLN                   ->  format leaf names (e.g. remove spaces in names)
       iTOL                  ->  Prepare iTOL-compatible files for tree visualization
    
    Genome databases
       get_GTDB_taxon_gnm    ->  Get id of genomes from specified GTDB taxons
       get_genome_GTDB       ->  Batch download GTDB genomes
       get_genome_NCBI       ->  Batch download GenBank genomes
       sampling_GTDB_gnms    ->  Select GTDB genomes from a specified taxon at specified sampling rank
       subset_GTDB_meta      ->  Subset metadata of GTDB reference genomes
        
    Visualization
       Plot_MAG              ->  plot MAGs, (GC vs depth)
       VisGeneFlk            ->  Visualize gene flanking regions
       SankeyTaxon           ->  Plot taxonomic classification with Sankey plot
       cross_link_seqs       ->  Cross link matched regions between two sequences
          
    Other modules
       split_folder          ->  Split folder
       js_cmds               ->  Commands to job scripts
       BestHit               ->  Keep best blast hits (outfmt 6)
       get_Pfam_hmms         ->  Get Pfam profiles by id
       Reads_simulator       ->  Simulate NGS reads
       usearch_uc            ->  Parse Usearch uc file
       exe_cmds              ->  Execute commands in a file (one cmd per line)
       get_MAG_reads_long    ->  Extract MAG-specific long reads for reassembling

Get help
---

1. Example usage for most of BioSAK modules can be found from their own help information, e.g.
        
       BioSAK COG2020 -h
       BioSAK iTOL -h 

1. Short tutorials for some BioSAK modules can be found here:

    https://github.com/songweizhi/BioSAK/tree/master/BioSAK_tutorial
