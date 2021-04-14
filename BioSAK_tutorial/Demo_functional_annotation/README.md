
## Short note for functional annotation against the COG, KEGG and CAZy databases using BioSAK (on Katana)


### Install BioSAK with Python's virtual environment 

BioSAK is a place where I put together some of my frequently used python scripts. 
You can install it on Katana with Python's virtual environment. 

1. Log into katana

       ssh z1234567@katana.restech.unsw.edu.au -o "ServerAliveInterval 10"
        
1. Start an interactive job (running programs on Katana head node is not allowed)    
        
       qsub -I -l nodes=1:ppn=6,mem=60gb,walltime=05:59:00

1. Create a Python3 virtual environment and install BioSAK

       module load python/3.7.3
       mkdir ~/mypython3env_BioSAK
       python3 -m venv --system-site-packages ~/mypython3env_BioSAK
       source ~/mypython3env_BioSAK/bin/activate
       pip3 install BioSAK
  
       # for future updating
       pip3 install --upgrade BioSAK
       
             
1. BioSAK is ready for running now. If you want to run BioSAK in the future, just run the following commands to activate python's virtual environment.

       module load python/3.7.3
       source ~/mypython3env_BioSAK/bin/activate
       BioSAK -h


### Prepare database files

1. COG (version 2020)

       # create a folder and download db files to this folder
       cd db_COG2020
       wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.fa.gz
       wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.cog.csv
       wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab
       wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab
       wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/Readme.2020-11-25.txt
       gunzip cog-20.fa.gz
        
       # make blast db
       module load blast+/2.11.0
       makeblastdb -in cog-20.fa -dbtype prot -parse_seqids -logfile cog-20.fa.log
        
       # make diamond db
       module load diamond/0.9.31
       diamond makedb --in cog-20.fa --db cog-20.fa.dmnd --quiet

1. dbCAN (version 9, released August 2020)

       cd db_dbCAN_V9
       wget http://bcb.unl.edu/dbCAN2/download/Databases/V9/hmmscan-parser.sh
       wget http://bcb.unl.edu/dbCAN2/download/Databases/V9/CAZyDB.07302020.fam-activities.txt
       wget http://bcb.unl.edu/dbCAN2/download/Databases/V9/dbCAN-HMMdb-V9.txt
       mv CAZyDB.07302020.fam-activities.txt CAZyDB.fam-activities.txt
       mv dbCAN-HMMdb-V9.txt dbCAN-fam-HMMs.txt
       module load hmmer/3.3
       hmmpress dbCAN-fam-HMMs.txt
       
1. KEGG

   + Please note that the KEGG database files (e.g. [genus_prokaryotes.pep.gz](https://www.kegg.jp/kegg/download/Readme/README.fasta)) are not for free, if you don't have permission to use these files, 
      you can use [BlastKOALA](https://www.kegg.jp/blastkoala/) or [GhostKOALA](https://www.kegg.jp/ghostkoala/) for KEGG annotation, which are free web-based online tools developed by the KEGG team.
      BlastKOALA and GhostKOALA assign K numbers to the user's sequence data by BLAST and GHOSTX searches, respectively, against a nonredundant set of KEGG GENES.
     
   + BlastKOALA: https://www.kegg.jp/blastkoala/ 
     (Up to **ten thousand** sequences may be uploaded)
   
   + GhostKOALA: https://www.kegg.jp/ghostkoala/ 
     (The file size of up to **300 MB** (one million sequences with average length of 300 or three million sequences with average length of 100) may be uploaded)
   
   + BlastKOALA Step-by-step Instructions: https://www.kegg.jp/blastkoala/help_blastkoala.html


### Functional annotation

1. COG

       module load python/3.7.3
       source ~/mypython3env/bin/activate
       module load diamond/0.9.31
       cd /srv/scratch/z5039045/Kelp_coassembly
       BioSAK COG2020 -i faa_files -x faa -db_dir path/to/db_COG2020 -m P -t 6 -diamond

1. CAZy (using dbCAN)

       module load python/3.7.3
       source ~/mypython3env/bin/activate
       module load hmmer/3.3
       cd /srv/scratch/z5039045/Kelp_coassembly
       BioSAK dbCAN -i faa_files -x faa -db_dir path/to/db_dbCAN_V9 -m P -t 6

1. KEGG

   If you have done KEGG annotation with BlastKOALA, you can still use BioSAK to summarize your annotation results (e.g. the number/percentage of query genes annotated to each KO at different levels).  
  
       module load python/3.7.3
       source ~/mypython3env/bin/activate
       module load diamond/0.9.24
       cd /srv/scratch/z5039045/Kelp_coassembly
       BioSAK KEGG -seq_in faa_files -x faa -db_dir /srv/scratch/z5039045/DB/KEGG_2016-09-26 -t 6 -diamond
       BioSAK KEGG -ko_in BlastKOALA_output_files -x txt -db_dir /srv/scratch/z5039045/DB/KEGG_2016-09-26 -t 6


### Output files (take COG as an example)

1. Annotation result for each query gene.

    | Query | ID | Category | Description |
    |:---:|:---:|:---:|---|
    | gene_1 | COG0801 | H | 7,8-dihydro-6-hydroxymethylpterin pyrophosphokinase |
    | gene_2 ||||
    | gene_3 | COG0514 | L | Superfamily II DNA helicase RecQ |
    | gene_4 | COG1309 | K | DNA-binding protein, AcrR family, includes nucleoid occlusion protein SlmA |
    | gene_5 | COG0514 | L | Superfamily II DNA helicase RecQ |

1. The number/percentage of genes annotated to each COG ID/category.

    | ID | GeneNumber | Description |
    |:---:|:---:|---|
    | COG0514 | 2 | Superfamily II DNA helicase RecQ |
    | COG0801 | 1 | 7,8-dihydro-6-hydroxymethylpterin pyrophosphokinase |
    | COG1309 | 1 | DNA-binding protein, AcrR family, includes nucleoid occlusion protein SlmA |

    | Category | GeneNumber | Description |
    |:---:|:---:|---|
    | L | 2 | Replication, recombination and repair |
    | H | 1 | Coenzyme transport and metabolism |
    | K | 1 | Transcription |

1. A data matrix of the number/percentage of genes annotated to each COG ID/category for each of your query file, if you annotated multiple sequence files together.

    | | COG0514 | COG0801 | COG1309 | ... |
    |:---:|:---:|:---:|:---:|:---:|
    | file_1 | 8 | 0 | 1 | ... |
    | file_2 | 3 | 0 | 2 | ... |
    | file_3 | 5 | 1 | 4 | ... |
    | file_4 | 3 | 9 | 1 | ... |
    | file_5 | 4 | 0 | 3 | ... |
    |  ...   |   |   |   |     |

       
### References and online resources:

+ dbCAN2: http://bcb.unl.edu/dbCAN2/blast.php


