
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

    Annotation modules
       Prodigal               ->   Wrapper for running Prodigal
       CheckM                 ->   Wrapper for running CheckM
       CheckM_op_parser       ->   Parse (combined) CheckM outputs
       COG2020                ->   COG annotation (v2020, by blastp/diamond)
       KEGG                   ->   KEGG annotation
       dbCAN                  ->   CAZy annotation with dbCAN
       NetEnzymes             ->   Get network of enzymes (based on MetaCyc, under development)   
       Enrichment             ->   Gene set enrichment analysis (to be added)
    
    16S related modules
       GTDB_16S               ->   Classify 16S against GTDB via best-hit approach
       SILVA_for_BLCA         ->   Prepare BLCA-compatible SILVA SSU database
       GTDB_for_BLCA          ->   Prepare BLCA-compatible GTDB SSU database
       
    Sequence manipulator
       gbk2fa                 ->   gbk to fasta
       gbk2ffn                ->   gbk to ffn
       gbk2faa                ->   gbk to faa
       ffn2faa                ->   ffn to faa
       get_rc                 ->   Get reverse complement sequence
       rename_seq             ->   Rename sequences in a file
       select_seq             ->   Select sequences by their id
       get_gene_depth         ->   Get gene depth by contig depth
       convert_align_format   ->   Convert alignment format
       OneLineAln             ->   One-line fasta format alignments
       SubsetAlnCols          ->   Subset MSA by column
       rename_reads_for_Reago ->   Rename paired reads for Reago
       MeanMappingDepth       ->   Get mean mapping depth 

    Tree manipulator
       get_SCG_tree           ->   Construct SCG tree for query genomes
       label_tree             ->   Add labels to tree leaves
       subset_tree            ->   Subset tree
       iTOL                   ->   Prepare iTOL-compatible files for tree visualization
         
    Other modules
       split_folder           ->   Split folder
       SankeyTaxon            ->   Plot taxonomic classification with Sankey plot
       BestHit                ->   Keep Best Hits only (blast outfmt 6)
       get_bin_abundance      ->   Get bin abundance
       sra_reads_downloader   ->   Download SRA read files
       dwnld_GenBank_genome   ->   Batch download GenBank genomes
       get_Pfam_hmms          ->   Get Pfam profiles by id
       reads_simulator        ->   Simulate NGS reads
       plot_sam_depth         ->   Plot sam depth
       reads2bam              ->   mapping and sorting
       sam2bam                ->   sam to bam with samtools
       VisGeneFlk             ->   visualize flanking regions of specified gene
       usearch_uc             ->   Usearch uc file parser

Get help
---

1. Example usage for most of BioSAK modules can be found from their specific help information, e.g.
        
       BioSAK COG2020 -h
       BioSAK iTOL -h 
       BioSAK VisGeneFlk -h 

1. Several short tutorials for specific modules of BioSAK can also be found here:

    https://github.com/songweizhi/BioSAK/tree/master/BioSAK_tutorial
