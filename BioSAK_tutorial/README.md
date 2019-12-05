
## Tutorial for running BioSAK on Katana

[![pypi licence ](https://img.shields.io/pypi/l/BioSAK.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version ](https://img.shields.io/pypi/v/BioSAK.svg)](https://pypi.python.org/pypi/BioSAK) 


### Install BioSAK with Python virtual environment 

1. Log into katana

       ssh z1234567@katana.restech.unsw.edu.au -o "ServerAliveInterval 10"
        
1. Start a interactive job (running programs on Katana head node is not allowed)    
        
       qsub -I -l nodes=1:ppn=4,mem=60gb,walltime=02:59:00

1. Create a Python3 virtual environment and install BioSAK

       module load python/3.7.3
       mkdir ~/mypython3env_BioSAK
       python3 -m venv --system-site-packages ~/mypython3env_BioSAK
       source ~/mypython3env_BioSAK/bin/activate
        
       # for the first time installation
       pip3 install BioSAK
  
       # for later updating
       pip3 install --upgrade BioSAK
       
       # to leave Python's virtual environment
       deactivate 
             
1. BioSAK is ready for running now. If you want to run BioSAK later, just run the following commands to activate the virtual environment.

       module load python/3.7.3
       source ~/mypython3env_BioSAK/bin/activate
       BioSAK -h


### Prepare database files

1. Download [qsub_prepare_DB.sh](https://github.com/songweizhi/BioSAK/blob/master/BioSAK_tutorial/qsub_prepare_DB.sh) and change the email address in line 5 to your own.

1. Upload it to your Katana Scratch and submit with qsub.

       cd /srv/scratch/z1234567
       qsub qsub_prepare_DB.sh


### Run BioSAK

1. Download [demo data](https://www.dropbox.com/s/ur9c0vsbndl5lop/BioSAK_demo.tar.gz?dl=0)

       cd /srv/scratch/z1234567
       wget https://www.dropbox.com/s/ur9c0vsbndl5lop/BioSAK_demo.tar.gz
       tar -xzvf BioSAK_demo.tar.gz
       cd BioSAK_demo
              
1. Predict genes from assemblies with Prokka

       module load perl/5.28.0
       module load infernal/1.1.2 
       module load blast+/2.2.31 
       module load hmmer/3.2.1
       module load prodigal/2.6.3
       module load tbl2asn/25.7 
       module load parallel/20190522 
       module load aragorn/1.2.38 
       module load prokka/1.13.3
       mkdir CtgSeq_Prokka    
       prokka --metagenome --prefix Kelp --locustag Kelp --outdir CtgSeq_Prokka/Kelp CtgSeq/Kelp_ctg.fa
       prokka --metagenome --prefix Sponge --locustag Sponge --outdir CtgSeq_Prokka/Sponge CtgSeq/Sponge_ctg.fa
       prokka --metagenome --prefix Seawater --locustag Seawater --outdir CtgSeq_Prokka/Seawater CtgSeq/Seawater_ctg.fa
       prokka --metagenome --prefix Sediment --locustag Sediment --outdir CtgSeq_Prokka/Sediment CtgSeq/Sediment_ctg.fa
    
    Prokka help info:    
    + --metagenome:       Improve gene predictions for highly fragmented genomes
    + --locustag:         Locus tag prefix (prefix of gene id)

1. copy faa and gff files into separate folders

       mkdir CtgSeq_faa
       cp CtgSeq_Prokka/*/*.faa CtgSeq_faa/
       
       mkdir CtgSeq_gff
       cp CtgSeq_Prokka/*/*.gff CtgSeq_gff/      

1. get gene depth according to the depth of the contig they sit in 
    
       module load python/3.7.3
       source ~/mypython3env_BioSAK/bin/activate
       BioSAK get_gene_depth -gff CtgSeq_gff/Kelp.gff -ctg_depth CtgDepth/Kelp_ctg.depth -skip_header
       BioSAK get_gene_depth -gff CtgSeq_gff/Sponge.gff -ctg_depth CtgDepth/Sponge_ctg.depth -skip_header
       BioSAK get_gene_depth -gff CtgSeq_gff/Seawater.gff -ctg_depth CtgDepth/Seawater_ctg.depth -skip_header
       BioSAK get_gene_depth -gff CtgSeq_gff/Sediment.gff -ctg_depth CtgDepth/Sediment_ctg.depth -skip_header

       # move generated protein depth files into a separate folder
       mkdir CtgSeq_faa_depth
       mv CtgSeq_gff/*.depth CtgSeq_faa_depth/

1. run COG, KEGG and CAZy annotation

       module load diamond/0.9.24
       module load hmmer/3.2.1
       	   
       zID="z1234567"
       BioSAK COG2014 -db_dir /srv/scratch/$zID/BioSAK_db/COG2014 -m P -t 4 -i CtgSeq_faa -x faa -diamond -depth CtgSeq_faa_depth
       BioSAK KEGG -db_dir /srv/scratch/$zID/BioSAK_db/KEGG -t 4 -seq_in CtgSeq_faa -x faa -diamond -depth CtgSeq_faa_depth
       BioSAK dbCAN -db_dir /srv/scratch/$zID/BioSAK_db/dbCAN -m P -t 4 -i CtgSeq_faa -x faa -depth CtgSeq_faa_depth


### Commands for getting contig depth with MetaBAT's jgi_summarize_bam_contig_depths

    # Mapping filtered reads back to their assemblies
    module load bowtie/2.3.4.2
    bowtie2-build -f Seawater.fa Seawater
    bowtie2 -x Seawater -1 Seawater_R1_Q25_P.fastq -2 Seawater_R2_Q25_P.fastq -S Seawater.sam -p 6 -q

    # turn SAM file into BAM file
    module load samtools/1.9
    samtools view -bS Seawater.sam -o Seawater.bam
    rm Seawater.sam
    samtools sort Seawater.bam -o Seawater_sorted.bam
    rm Seawater.bam
    samtools index Seawater_sorted.bam

    # get contig depth from BAM file
    module load metabat/2.12.1
    jgi_summarize_bam_contig_depths --outputDepth Seawater_ctg.depth Seawater_sorted.bam


### Some other modules of BioSAK

1. Download GenBank genomes

    + Go to https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/refseq_category:reference
    + Search genomes you want to download (e.g. Alteromonas, Archaea)
    + Click "Download" on the top right corner of the table
    + Download genome with dwnld_GenBank_genome
    
           cd /srv/scratch/$zID/BioSAK_demo/OtherFiles
           BioSAK dwnld_GenBank_genome -csv Alteromonas.csv -fna -faa -gbff -name

1. For more:

       BioSAK -h


### References and online resources:

+ WebMGA: http://weizhong-lab.ucsd.edu/webMGA/server/cog/
+ BlastKOALA: https://www.kegg.jp/blastkoala/
+ GhostKOALA: https://www.kegg.jp/ghostkoala/
+ dnCAN2: http://bcb.unl.edu/dbCAN2/blast.php
+ jgi_summarize_bam_contig_depths: https://bitbucket.org/berkeleylab/metabat/issues/36/how-the-depth-of-contig-was-calculated-why

