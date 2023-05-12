+ Dependencies
  + hmmer

+ annotate protein sequences

      BioSAK dbCAN -m P -t 6 -db_dir path/to/your/dbCAN_db_dir -i recipient.faa 
      BioSAK dbCAN -m P -t 6 -db_dir path/to/your/dbCAN_db_dir -i recipient.faa -depth recipient.depth
      BioSAK dbCAN -m P -t 6 -db_dir path/to/your/dbCAN_db_dir -i faa_files -x faa
      BioSAK dbCAN -m P -t 6 -db_dir path/to/your/dbCAN_db_dir -i faa_files -x faa -depth depth_files

+ annotate DNA sequences

      BioSAK dbCAN -m N -t 6 -db_dir /srv/scratch/z5039045/DB/dbCAN_V9 -i recipient.ffn
      BioSAK dbCAN -m N -t 6 -db_dir /srv/scratch/z5039045/DB/dbCAN_V9 -i ffn_files -x ffn

+ Prepare DB files (versions V9):

      cd path/to/your/dbCAN_db_dir
      wget http://bcb.unl.edu/dbCAN2/download/Databases/V9/hmmscan-parser.sh
      wget http://bcb.unl.edu/dbCAN2/download/Databases/V9/CAZyDB.07302020.fam-activities.txt
      wget http://bcb.unl.edu/dbCAN2/download/Databases/V9/dbCAN-HMMdb-V9.txt
      mv CAZyDB.07302020.fam-activities.txt CAZyDB.fam-activities.txt
      mv dbCAN-HMMdb-V9.txt dbCAN-fam-HMMs.txt
      hmmpress dbCAN-fam-HMMs.txt

+ How it works:
1. http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/readme.txt
2. The TotalDepth of a CAZy family is calculated by summing up the depth of all genes assigned to it.
3. The percentage of GeneNumber/TotalDepth of genes assigned to a CAZy family is calculated by dividing it 
   by the summation of GeneNumber/TotalDepth of all identified CAZy families. 

+ Note!!!
If you run dbCAN for multiple files in a batch manner and want to have their depth info incorporated into the results, 
you need to provide a folder containing individual depth files for each of your input sequence file.
Name of the depth file needs to be exactly the same as its corresponding sequence file, except the extension which is ".depth".

+ Depth file format (one gene per line, tab separated)
gene_1	30
gene_2	10.58