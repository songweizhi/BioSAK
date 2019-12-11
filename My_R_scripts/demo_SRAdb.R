# SRAdb
source (" http://bioconductor.org/biocLite.R")
biocLite("SRAdb")
library(SRAdb)

sqlfile = getSRAdbFile()
sra_con = dbConnect(SQLite(), sqlfile)


#rs <- getSRA (search_terms ='"antarctic", "metagenome"', out_types=c('run','study'), sra_con=sra_con, acc_only=0)
rs <- getSRA (search_terms ='"soil", "metagenome"', out_types=c('run','study'), sra_con=sra_con, acc_only=0)

write.csv(rs, "soil_metagenome_output.txt")


## Download sra files from NCBI SRA using ftp protocol:
getSRAfile( in_acc = c("SRR000648","SRR000657"), sra_con = sra_con, destDir = getwd(), fileType = 'sra' )



?getSRA
