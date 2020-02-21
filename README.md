
## BioSAK (A Swiss Army Knife for Biologists)

[![pypi licence ](https://img.shields.io/pypi/l/BioSAK.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version ](https://img.shields.io/pypi/v/BioSAK.svg)](https://pypi.python.org/pypi/BioSAK) 

Contact
---

+ **Weizhi Song**, Postdoctoral Researcher
+ Center for Marine Science & Innovation, University of New South Wales, Sydney, Australia
+ E-mail: songwz03@gmail.com

Dependencies
---

Module dependent, see module-specific help info for details. E.g.

    BioSAK KEGG -h
    BioSAK COG2014 -h

Installation
---

1. BioSAK has been tested on Linux/Mac, but **NOT** supported on Windows.

1. BioSAK is implemented in python3, you can install it with pip3:

       # for the first time installation
       pip3 install BioSAK
      
       # for later updating
       pip3 install --upgrade BioSAK
      
1. For UNSW Katana users

       ############# install BioSAK with virtual environment ############

       module load python/3.7.3
       mkdir ~/mypython3env_BioSAK
       python3 -m venv --system-site-packages ~/mypython3env_BioSAK
       source ~/mypython3env_BioSAK/bin/activate
       pip3 install BioSAK

       ####################### For later running ########################

       module load python/3.7.3
       source ~/mypython3env_BioSAK/bin/activate
       BioSAK -h

BioSAK modules
---

    Annotation modules
       Prodigal               ->   Wrapper for running Prodigal
       CheckM                 ->   Wrapper for running CheckM
       CheckM_op_parser       ->   Parse (combined) CheckM outputs
       COG2003                ->   COG annotation (v2003, by rpsblast)
       COG2014                ->   COG annotation (v2014, by blastp/diamond)
       KEGG                   ->   KEGG annotation
       dbCAN                  ->   CAZy annotation with dbCAN
       NetEnzymes             ->   Get network of enzymes (based on MetaCyc)   
    
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

    Tree manipulator
       get_SCG_tree           ->   Construct SCG tree for query genomes
       label_tree             ->   Add labels to tree leaves
       subset_tree            ->   Subset tree
       iTOL                   ->   Plot tree with iTOL
                      
    Other modules
       split_folder           ->   Split folder
       SankeyTaxon            ->   Plot taxonomic classification with Sankey plot
       BestHit                ->   Keep Best Hits only (blast outfmt 6)
       get_bin_abundance      ->   Get bin abundance
       sra_reads_downloader   ->   Download SRA read files
       dwnld_GenBank_genome   ->   Batch download GenBank genomes
       get_Pfam_hmms          ->   Get Pfam profiles by id
       reads_simulator        ->   Simulate NGS reads

Help information
---

    BioSAK -h
    BioSAK COG2014 -h
    BioSAK select_seq -h
    BioSAK dwnld_GenBank_genome -h
    etc.
