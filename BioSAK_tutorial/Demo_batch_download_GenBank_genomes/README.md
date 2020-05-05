Download GenBank genomes

   + Go to https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/refseq_category:reference
   + Search genomes you want to download (e.g. Alteromonas, Archaea)
   + Click "Download" on the top right corner of the table
   + Download genome with dwnld_GenBank_genome
    
         cd /srv/scratch/$zID/BioSAK_demo/OtherFiles
         BioSAK dwnld_GenBank_genome -csv Alteromonas.csv -fna -faa -gbff -name
