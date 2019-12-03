#!/bin/bash
#PBS -l nodes=1:ppn=3
#PBS -l mem=30gb
#PBS -l walltime=11:59:00
#PBS -M weizhi.song@unsw.edu.au
#PBS -j oe
#PBS -m ae

module load diamond/0.9.24


# Change into directory that you submitted the job from
cd $PBS_O_WORKDIR


# create folders
mkdir DB_BioSAK
mkdir DB_BioSAK/COG2014
mkdir DB_BioSAK/KEGG
mkdir DB_BioSAK/dbCAN


################################## prepare COG db files ##################################

cd $PBS_O_WORKDIR/DB_BioSAK/COG2014
wget ftp://ftp.ncbi.nlm.nih.gov//pub/COG/COG2014/data/prot2003-2014.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov//pub/COG/COG2014/data/prot2003-2014.tab
wget ftp://ftp.ncbi.nlm.nih.gov//pub/COG/COG2014/data/cog2003-2014.csv
wget ftp://ftp.ncbi.nlm.nih.gov//pub/COG/COG2014/data/cognames2003-2014.tab
wget ftp://ftp.ncbi.nlm.nih.gov//pub/COG/COG2014/data/fun2003-2014.tab        
gunzip prot2003-2014.fa.gz

# make diamond db
diamond makedb --in prot2003-2014.fa --db prot2003-2014.fa.dmnd --quiet


################################## prepare KEGG db files #################################

cd $PBS_O_WORKDIR/DB_BioSAK/KEGG

# download ko00001.keg
wget https://www.dropbox.com/s/8dmhocivgpe0pok/ko00001.keg

# download, decompress and rename prokaryotes.dat.gz
wget https://www.dropbox.com/s/4mu8v0xbp4e8pdo/prokaryotes.dat.gz
gunzip prokaryotes.dat.gz
mv prokaryotes.dat kegg_db_seq2ko.txt

# download, decompress and rename prokaryotes.pep.gz
wget https://www.dropbox.com/s/04uqf2ju6mei457/prokaryotes.pep.gz
gunzip prokaryotes.pep.gz
mv prokaryotes.pep kegg_db_seq.fasta

# make diamond db
diamond makedb --in kegg_db_seq.fasta --db kegg_db_seq.fasta.dmnd --quiet


################################# prepare CAZy db files ##################################

cd $PBS_O_WORKDIR/DB_BioSAK/dbCAN
wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/hmmscan-parser.sh
wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/dbCAN-fam-HMMs.txt
wget http://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.07312019.fam-activities.txt
mv CAZyDB.07312019.fam-activities.txt CAZyDB.fam-activities.txt

