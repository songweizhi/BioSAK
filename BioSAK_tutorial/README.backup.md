
## Tutorial for running BioSAK on Katana

[![pypi   licence        ](https://img.shields.io/pypi/l/BioSAK.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi   version        ](https://img.shields.io/pypi/v/BioSAK.svg)](https://pypi.python.org/pypi/BioSAK) 


### Setup Python virtual environment


1. Log into katana

        ssh zID@katana.restech.unsw.edu.au -o "ServerAliveInterval 10"
        
1. Start a interactive job (running programs on Katana head node is not allowed)    
        
        qsub -I -l nodes=1:ppn=6,mem=60gb,walltime=02:59:00

1. Create a Python3 virtual environment and install BioSAK

        module load python/3.7.3
        mkdir ~/mypython3env_BioSAK
        python3 -m venv --system-site-packages ~/mypython3env_BioSAK
        source ~/mypython3env_BioSAK/bin/activate
        
        # for the first time installation
        pip3 install BioSAK
  
        # for later updating
        pip3 install --upgrade BioSAK

1. get help information of BioSAK

        BioSAK -h
        BioSAK COG2014 -h
        BioSAK KEGG -h
        BioSAK dbCAN -h


### Prepare database files for COG, KEGG and CAZy (dbCAN) annotation

1. Create a DB folder in your scratch, all database files needed for COG, KEGG and CAZy annotation 
will be stored in separate folders within DB.
        
        cd /srv/scratch/zID
        mkdir BioSAK_db
        mkdir BioSAK_db/COG2014
        mkdir BioSAK_db/KEGG
        mkdir BioSAK_db/dbCAN
 
1. Prepare COG2014 database files

        cd /srv/scratch/zID/BioSAK_db/COG2014
        wget ftp://ftp.ncbi.nlm.nih.gov//pub/COG/COG2014/data/prot2003-2014.fa.gz
        wget ftp://ftp.ncbi.nlm.nih.gov//pub/COG/COG2014/data/prot2003-2014.tab
        wget ftp://ftp.ncbi.nlm.nih.gov//pub/COG/COG2014/data/cog2003-2014.csv
        wget ftp://ftp.ncbi.nlm.nih.gov//pub/COG/COG2014/data/cognames2003-2014.tab
        wget ftp://ftp.ncbi.nlm.nih.gov//pub/COG/COG2014/data/fun2003-2014.tab        
        gunzip prot2003-2014.fa.gz
        
1. Prepare KEGG database files

        cd /srv/scratch/zID/BioSAK_db/KEGG
        
        # prokaryotes.pep.gz and prokaryotes.dat.gz (not free)
        https://www.kegg.jp/kegg/download/Readme/README.fasta
        
        # ko00001.keg
        wget https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir=
 
1. Prepare dbCAN database files

        cd /srv/scratch/zID/BioSAK_db/dbCAN
        wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/hmmscan-parser.sh
        wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/dbCAN-fam-HMMs.txt
        wget http://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.07312019.fam-activities.txt
        mv CAZyDB.07312019.fam-activities.txt CAZyDB.fam-activities.txt
        
        
### Get gene depth

1. Mapping reads back to assemblies and get their depth with MetaBAT's jgi_summarize_bam_contig_depths

        module load metabat/0.32.4
        jgi_summarize_bam_contig_depths --outputDepth Ctg_depth.txt sample_1.bam
        
2. Predict genes from assemblies with Prokka

        module load perl/5.28.0
        module load infernal/1.1.2 
        module load blast+/2.2.31 
        module load hmmer/3.2.1
        module load prodigal/2.6.3
        module load tbl2asn/25.7 
        module load parallel/20190522 
        module load aragorn/1.2.38 
        module load prokka/1.13.3
        prokka --force --metagenome --prefix Sample_1 --locustag Sample_1 --outdir Sample_1 Sample1_ctgs.fa

3. Get gene depth according to contig depth (together with gbk file)

        BioSAK get_gene_depth -gbk Sample_1.gbk -ctg_depth Ctg_depth.txt -skip_header
        

### Run BioSAK

1. Download demo dataset

    https://www.dropbox.com/sh/cllqvy6uuw7e3oj/AAA-WqNhAcNwsALKXdVnby8qa?dl=0

1. COG annotation

        module load diamond/0.9.24
        BioSAK COG2014 -db_dir /srv/scratch/zID/BioSAK_db/COG2014 -m P -t 4 -i faa_files -x faa -diamond -depth FiveGenomes_depth
        
1. KEGG annotation

        module load diamond/0.9.24
        BioSAK KEGG -db_dir /srv/scratch/zID/BioSAK_db/KEGG -t 4 -seq_in faa_files -x faa -diamond -depth FiveGenomes_depth
        
        # or, if you already have you proteins annotated with BlastKOALA/GhostKOALA
        BioSAK KEGG -db_dir /srv/scratch/zID/BioSAK_db/KEGG -t 4 -ko_in ko_files -x faa -diamond -depth FiveGenomes_depth

1. CAZy annotation

        module load hmmer/3.2.1
        BioSAK dbCAN -db_dir /srv/scratch/zID/BioSAK_db/dbCAN -m P -t 4 -i faa_files -x faa -depth FiveGenomes_depth
        
    
### References and online resources:

+ WebMGA: http://weizhong-lab.ucsd.edu/webMGA/server/cog/
+ BlastKOALA: https://www.kegg.jp/blastkoala/
+ GhostKOALA: https://www.kegg.jp/ghostkoala/
+ dnCAN2: http://bcb.unl.edu/dbCAN2/blast.php
+ jgi_summarize_bam_contig_depths: https://bitbucket.org/berkeleylab/metabat/issues/36/how-the-depth-of-contig-was-calculated-why

