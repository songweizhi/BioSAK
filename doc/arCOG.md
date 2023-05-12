### Dependencies: 
  + blast+, 
  + diamond


### Prepare DB files:

    cd path/to/your/arCOG_db_dir
    wget https://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/tmp.ar18/ar18.fa.gz
    wget https://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/tmp.ar18/ar18.ar14.02.csv
    wget https://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/tmp.ar18/arCOG_names_220807.txt
    wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab
    gunzip ar18.fa.gz
    makeblastdb -in ar18.fa -dbtype prot -parse_seqids -logfile ar18.fa.log
    diamond makedb --in ar18.fa --db ar18.fa.dmnd --quiet

### Annotate protein sequences

    BioSAK arCOG -m P -t 6 -db_dir /Users/songweizhi/DB/arCOG18 -i genes.faa 
    BioSAK arCOG -m P -t 6 -db_dir /Users/songweizhi/DB/arCOG18 -i faa_files -x faa -depth depth_files

### Annotate DNA sequences (ORFs)

    BioSAK arCOG -m N -t 6 -db_dir /Users/songweizhi/DB/arCOG18 -i genes.ffn -depth gene.depth
    BioSAK arCOG -m N -t 6 -db_dir /Users/songweizhi/DB/arCOG18 -i ffn_files -x ffn

+ For how it works, please refers to the help information of the COG2020 module (BioSAK COG2020 -h)
