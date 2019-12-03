
## Tutorial for running BioSAK on Katana

[![pypi licence ](https://img.shields.io/pypi/l/BioSAK.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version ](https://img.shields.io/pypi/v/BioSAK.svg)](https://pypi.python.org/pypi/BioSAK) 


### Setup Python virtual environment


1. Log into katana

       ssh zID@katana.restech.unsw.edu.au -o "ServerAliveInterval 10"
        
1. Start a interactive job (running programs on Katana head node is not allowed)    
        
       qsub -I -l nodes=1:ppn=6,mem=60gb,walltime=02:59:00

1. Create a Python3 virtual environment on Katana and install BioSAK

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
       
       
1. get help information of BioSAK

       module load python/3.7.3
       source ~/mypython3env_BioSAK/bin/activate
       BioSAK -h


### Prepare database files

1. Download [qsub_prepare_DB.sh](https://github.com/songweizhi/BioSAK/blob/master/BioSAK_tutorial/qsub_prepare_DB.sh)

1. Change the email address in line 5 to your own

1. Upload it to your Katana Scratch and submit with qsub

       cd /srv/scratch/$zID
       qsub qsub_prepare_DB.sh


### Run BioSAK

1. Download demo data

       zID="z5039045"
       cd /srv/scratch/$zID
       wget https://www.dropbox.com/s/ah99bsqi1cb7592/BioSAK_demo.tar.gz
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
       mkdir Metagenomic_assemblies_Prokka
       prokka --force --metagenome --prefix Kelp --locustag Kelp --outdir Metagenomic_assemblies_Prokka/Kelp Metagenomic_assemblies/Kelp_ctg.fa
       prokka --force --metagenome --prefix Sponge --locustag Sponge --outdir Metagenomic_assemblies_Prokka/Sponge Metagenomic_assemblies/Sponge_ctg.fa
       prokka --force --metagenome --prefix Seawater --locustag Seawater --outdir Metagenomic_assemblies_Prokka/Seawater Metagenomic_assemblies/Seawater_ctg.fa
       prokka --force --metagenome --prefix Sediment --locustag Sediment --outdir Metagenomic_assemblies_Prokka/Sediment Metagenomic_assemblies/Sediment_ctg.fa

    note:
    --metagenome       Improve gene predictions for highly fragmented genomes
    --locustag         Locus tag prefix (prefix of gene id)
    --prefix [X]       Filename output prefix

1. copy faa and gff files into separate folders

       mkdir Metagenomic_assemblies_faa
       cp Metagenomic_assemblies_Prokka/*/*.faa Metagenomic_assemblies_faa
       
       mkdir Metagenomic_assemblies_gff
       cp Metagenomic_assemblies_Prokka/*/*.gff Metagenomic_assemblies_gff      

1. get gene depth according to the depth of the contig they sit in 
    
       module load python/3.7.3
       source ~/mypython3env_BioSAK/bin/activate
       BioSAK get_gene_depth -gff Metagenomic_assemblies_gff/Kelp.gff -ctg_depth Metagenomic_assemblies_depth/Kelp_ctg.depth -skip_header
       BioSAK get_gene_depth -gff Metagenomic_assemblies_gff/Sponge.gff -ctg_depth Metagenomic_assemblies_depth/Sponge_ctg.depth -skip_header
       BioSAK get_gene_depth -gff Metagenomic_assemblies_gff/Seawater.gff -ctg_depth Metagenomic_assemblies_depth/Seawater_ctg.depth -skip_header
       BioSAK get_gene_depth -gff Metagenomic_assemblies_gff/Sediment.gff -ctg_depth Metagenomic_assemblies_depth/Sediment_ctg.depth -skip_header

       # move generated protein depth files into a separate folder
       mv Metagenomic_assemblies_gff/*.depth Metagenomic_assemblies_faa_depth/

1. run COG, KEGG and CAZy annotation

       module load diamond/0.9.24
       module load hmmer/3.2.1
       	   
       BioSAK COG2014 -db_dir /srv/scratch/$zID/BioSAK_db/COG2014 -m P -t 4 -i Metagenomic_assemblies_faa -x faa -diamond -depth Metagenomic_assemblies_faa_depth
       
       BioSAK KEGG -db_dir /srv/scratch/$zID/BioSAK_db/KEGG -t 4 -seq_in Metagenomic_assemblies_faa -x faa -diamond -depth Metagenomic_assemblies_faa_depth
       
       BioSAK dbCAN -db_dir /srv/scratch/$zID/BioSAK_db/dbCAN -m P -t 4 -i Metagenomic_assemblies_faa -x faa -depth Metagenomic_assemblies_faa_depth


### Steps for getting contig depth with MetaBAT's jgi_summarize_bam_contig_depths

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


### References and online resources:

+ WebMGA: http://weizhong-lab.ucsd.edu/webMGA/server/cog/
+ BlastKOALA: https://www.kegg.jp/blastkoala/
+ GhostKOALA: https://www.kegg.jp/ghostkoala/
+ dnCAN2: http://bcb.unl.edu/dbCAN2/blast.php
+ jgi_summarize_bam_contig_depths: https://bitbucket.org/berkeleylab/metabat/issues/36/how-the-depth-of-contig-was-calculated-why
