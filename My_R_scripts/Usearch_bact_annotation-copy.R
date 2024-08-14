# 16S script
# author: "JL"
# date: "22/6/2023"
# 
# This script includes the use of the new database GTDB and BLCA method.
# ---

#install the seqinr package first then load libraries and usearch software
#install.packages("seqinr")
#install.packages("dplyr")
  
library(seqinr)
library(dplyr)

system("ln -s 00_DataNeeded/usearch/usearch11.0.667_i86osx64 usearch11.0.667 ; chmod +x usearch11.0.667") # MacOS

system('./usearch11.0.667') # Usearch info




# usearch v11.0.667_i86osx64, 34.4Gb RAM, 10 cores
# (C) Copyright 2013-18 Robert C. Edgar, all rights reserved.
# https://drive5.com/usearch
# 
# License: t.thomas@unsw.edu.au, non-profit use, max 1 process(es)



library(tidyverse)
# install.packages("reshape2")
library(reshape2)

##############################################################################################################
#################################################### GTDB BLCA ###############################################
##############################################################################################################

# add clustalo and muscle to system path (required by BLCA)
system("cp ./00_DataNeeded/clustalo/clustal-omega-1.2.3-macosx ./00_DataNeeded/clustalo/clustalo ; chmod +x ./00_DataNeeded/clustalo/clustalo") # MacOS
system("cp ./00_DataNeeded/muscle/muscle3.8.31_i86darwin64 ./00_DataNeeded/muscle/muscle ; chmod +x ./00_DataNeeded/muscle/muscle") # MacOS

Sys.setenv(PATH = paste(Sys.getenv("PATH"), "./00_DataNeeded/clustalo", sep = ":"))
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "./00_DataNeeded/muscle", sep = ":"))


# run BLCA agsinst GTDB SSU database

# If you saw a "BioPython is not detected!" error, decomment the following line to install biopython
# system(paste(system("which python3", intern = TRUE), '-m pip install biopython --user'))


# r214.1


# system('python3 ./00_DataNeeded/BLCA/2.blca_main.py  -r ./00_DataNeeded/db_GTDB214.1/ssu_all_r214_BLCAparsed.taxonomy -q ./00_DataNeeded/db_GTDB214.1/ssu_all_r214_BLCAparsed.fasta -i 9_UchimeReference/AllSamples_unoise_nc.fasta -o 11_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.1.txt')#-c 0.95 --iset 0.95
system('python3 ./00_DataNeeded/BLCA/2.blca_main.py  -r ./00_DataNeeded/db_GTDB214.1/ssu_all_r214_BLCAparsed.taxonomy -q ./00_DataNeeded/db_GTDB214.1/ssu_all_r214_BLCAparsed.fasta -i b_rep_set.fasta -o 11_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.1.txt')#-c 0.95 --iset 0.95

# 
# clustalo is located in your PATH!
#   >  > Fasta file read in!
#   >  > Reading in taxonomy information! ....
# blastn is located in your PATH!
#   > > Running blast!!
#   > > Blastn Finished!!
#   >  > read in blast file...
# >  > blastn file opened
# >  > blast output read in
# >  > Start aligning reads...
# >  > Taxonomy file generated!!
#   Time elapsed: 1832 minutes

# remove taxonomic ranks from BLCA's classification and put number in parenthesis (make it easy to read)
m = 1
for (each_line in readLines(file('11_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.1.txt',open="r")) ){
  each_line_split = strsplit(each_line, '\t')
  OTU_ID = each_line_split[[1]][1]
  taxonomy = each_line_split[[1]][2]
  taxonomy_split = strsplit(taxonomy, ';')
  taxonomy_no_rank = ''
  n = 1
  for (taxon in taxonomy_split[[1]]){
    if (n%%2 == 1){
      taxon_split = strsplit(taxon, ':')
      if (length(taxon_split[[1]]) ==2)
      {taxon_no_rank = taxon_split[[1]][2]} 
      else 
      {taxon_no_rank = taxon_split[[1]][1]}
      taxonomy_no_rank = paste(taxonomy_no_rank, taxon_no_rank, sep = ");")} 
    else 
    {taxonomy_no_rank = paste(taxonomy_no_rank, taxon, sep = "(")}
    n = n + 1
  }
  taxonomy_no_rank = paste(taxonomy_no_rank, ')', sep = "")
  taxonomy_no_rank = substr(taxonomy_no_rank, 3, nchar(taxonomy_no_rank))
  if (taxonomy_no_rank == "Unclassified)"){taxonomy_no_rank = "Unclassified"}
  
  taxonomy_no_rank_with_OTU = paste(OTU_ID, taxonomy_no_rank, sep = "\t")
  
  if (m == 1)
  {cat(taxonomy_no_rank_with_OTU,file='11_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.2.txt',sep="\n", append=FALSE)} 
  else 
  {cat(taxonomy_no_rank_with_OTU,file='11_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.2.txt',sep="\n", append=TRUE)}
  m = m +1
}

 
# # Merge Table and Taxonomy
dir.create("12_FinalOtuTableGTDB_BLCA")
OTU = read.delim("b_allsample_notaxa.txt", header = T)
#I manually removed the numbers in parentesis after the taxonomy e.g. (100.0000). Then I saved the file as "AllSamples_unoise_BLCA_out.2filtered2_split.txt". here the taxonomy is split.
TAX = read.delim("11_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.2.txt", header = F)
names(TAX) = c("X.OTU.ID", "Taxonomy")
#
OTU_TAX = merge(OTU, TAX, by = "X.OTU.ID")
write.table(OTU_TAX, "12_FinalOtuTableGTDB_BLCA/AllSamples_unoise_otu_table_BLCA_filtered2.txt", sep = "\t", row.names = F, col.names = T, quote = F)
#


