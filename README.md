
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
       Prodigal                ->   Wrapper for running Prodigal
       Prokka                  ->   Wrapper for running Prokka
       CheckM                  ->   Wrapper for running CheckM
       CheckM_op_parser        ->   Parse (combined) CheckM outputs
       COG2020                 ->   COG annotation (v2020, by blastp/diamond)
       arCOG                   ->   to be added
       KEGG                    ->   KEGG annotation
       dbCAN                   ->   CAZy annotation with dbCAN

    16S rRNA related
       top_16S_hits            ->   Classify 16S by top-blast-hits approach
       SILVA_for_BLCA          ->   Prepare BLCA-compatible SILVA SSU database
       GTDB_for_BLCA           ->   Prepare BLCA-compatible GTDB SSU database
       BLCA_op_parser          ->   Make the BLCA outputs bit easier to read
    
    Sequence manipulator
       gbk2fa                  ->   gbk to fasta
       gbk2ffn                 ->   gbk to ffn
       gbk2faa                 ->   gbk to faa
       ffn2faa                 ->   ffn to faa
       get_rc                  ->   Get reverse complement sequence
       slice_seq               ->   Get specified region of a sequence
       rename_seq              ->   Rename sequences in a file
       select_seq              ->   Select sequences by their ids
       split_fasta             ->   Split a fasta file into multiple files
       get_gene_depth          ->   Get gene depth by contig depth
       convert_align_format    ->   Convert alignment format
       OneLineAln              ->   One-line fasta format alignments
       SubsetAlnCols           ->   Subset MSA by column
       rename_reads_for_Reago  ->   Rename paired reads for Reago
       MeanMappingDepth        ->   Get mean mapping depth 
    
    SAM and BAM
       plot_sam_depth          ->   Plot SAM depth
       reads2bam               ->   Mapping and sorting
       sam2bam                 ->   Sam to BAM with samtools
       split_sam               ->   Split SAM/BAM file by reference

    Tree manipulator
       GTDB_tree_r207          ->   Infer GTDB (r207) archaeal/bacterial tree
       label_tree              ->   Add labels to tree leaves (does not work)
       subset_tree             ->   Subset tree
       iTOL                    ->   Prepare iTOL-compatible files for tree visualization
       compare_trees           ->   Compare trees with mantel test
    
    Genome databases
       get_GTDB_taxon_gnm      ->   Get id of all genomes from specified GTDB taxons
       get_genome_GTDB         ->   Batch download GTDB representative genomes
       get_genome_NCBI         ->   Batch download GenBank genomes
   
    Visualization
       SankeyTaxon             ->   Plot taxonomic classification with Sankey plot
       VisGeneFlk              ->   Visualize gene flanking regions
       Plot_MAG                ->   plot MAGs, (GC vs depth)
                 
    Other modules
       split_folder            ->   Split folder
       js_cmds                 ->   Commands to job scripts
       BestHit                 ->   Keep Best Hits only (blast outfmt 6)
       get_bin_abundance       ->   Get bin abundance
       mean_MAG_cov            ->   MAG depth based on weighted mean contig depth 
       get_Pfam_hmms           ->   Get Pfam profiles by id
       Reads_simulator         ->   Simulate NGS reads
       usearch_uc              ->   Usearch uc file parser
       get_gnm_size            ->   Get the total length of genome(s)
       exe_cmds                ->   Execute commands in a file (one per line)

Get help
---

1. Example usage for most of BioSAK modules can be found from their own help information, e.g.
        
       BioSAK COG2020 -h
       BioSAK iTOL -h 

1. Short tutorials for some BioSAK modules can be found here:

    https://github.com/songweizhi/BioSAK/tree/master/BioSAK_tutorial
