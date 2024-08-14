# TEST

# Please decomment
#system("cp 00_DataNeeded/usearch/usearch11.0.667_i86linux64 usearch11.0.667 ; chmod +x usearch11.0.667") # Linux
#system("cp 00_DataNeeded/usearch/usearch11.0.667_i86osx64 usearch11.0.667 ; chmod +x usearch11.0.667") # MacOS
#shell("copy 00_DataNeeded\\usearch\\usearch11.0.667_win64.exe usearch11.0.667") # Windows

system('./usearch11.0.667') # Test

# Quality filtering
dir.create("02_PreProcessedData")
dir.create("03_MergedData")
dir.create("04_FilteredData")
dir.create("05_TrimmedData")
dir.create("06_DereplicatedData")

files = dir(path = '01_RawData/', pattern = '_R1_')
for (file1 in files)
{
  file2 = gsub(x = file1, pattern = '_R1_', replacement = '_R2_')
  filename = gsub(x = file1, pattern = "_R1_001.fastq.*", replacement = '', perl = T)
  
  #	Quality trimming
  command = paste('java -jar 00_DataNeeded/Trimmomatic-0.38/trimmomatic-0.38.jar  PE 01_RawData/', file1, " 01_RawData/", file2, " 02_PreProcessedData/", filename, "_pF.fastq 02_PreProcessedData/", filename, "_upF.fastq 02_PreProcessedData/",filename, "_pR.fastq 02_PreProcessedData/",filename, "_upR.fastq SLIDINGWINDOW:4:12 MINLEN:100", sep = "")
  system(command)
  
  #	Merge paired reads (range needs to be adjusted)
  command = paste("./usearch11.0.667 -fastq_mergepairs 02_PreProcessedData/", filename, "_pF.fastq -reverse 02_PreProcessedData/", filename, "_pR.fastq -fastqout 03_MergedData/",filename, ".fastq -relabel @ -fastq_maxdiffs 5 -fastq_pctid 80 -fastq_minmergelen 300 -fastq_maxmergelen 500 -sample ", filename, sep = "")
  system(command)
  
  #	Quality filtering
  command = paste("./usearch11.0.667  -fastq_filter 03_MergedData/", filename, ".fastq -fastaout 04_FilteredData/", filename, ".fasta -fastq_maxns 1 -fastq_maxee 1", sep = "")
  system(command)

  #	Remove primers
  command = paste("./usearch11.0.667 -fastx_truncate 04_FilteredData/",filename, ".fasta -stripleft 17 -stripright 17 -fastaout 05_TrimmedData/",filename, ".fasta", sep = "")
  system(command)
  
  #	Dereplication
  command = paste("./usearch11.0.667  -fastx_uniques 05_TrimmedData/", filename, ".fasta -fastaout 06_DereplicatedData/", filename, ".fasta -sizeout", sep = "")
  system(command)
}

rm(list=ls())

#Check merging range
library(seqinr)
seq_lengths = NULL
for(f in dir(path = '05_TrimmedData/', full.names = T)[1:3])
{
  print(f)
  fasta_file = read.fasta(f)
  seq_lengths = c(seq_lengths, as.numeric(lengths(fasta_file)))
}
hist(seq_lengths) # USe in line 27/28 to adjust range
quantile(seq_lengths)


#	Join
dir.create("07_JoinedData")

# Please decomment
#system("cat 06_DereplicatedData/*.fasta > 07_JoinedData/AllSamples.fasta") # Linux and MacOS
#shell("type 06_DereplicatedData\\*.fasta > 07_JoinedData\\AllSamples.fasta") # Windows

#	Dereplication

dir.create("08_UniqueSequences")

system("./usearch11.0.667 -fastx_uniques 07_JoinedData/AllSamples.fasta -fastaout 08_UniqueSequences/AllSamples_uniques.fasta -sizein -sizeout -strand both")

#	Generating unique sequences using UNOISE
dir.create("09_DenoisedSequences")

system("./usearch11.0.667 -unoise3 08_UniqueSequences/AllSamples_uniques.fasta -zotus 09_DenoisedSequences/AllSamples_denoised.fasta")

#	Chimera Removal
dir.create("10_UchimeReference")

system("./usearch11.0.667 -uchime2_ref 09_DenoisedSequences/AllSamples_denoised.fasta -db 00_DataNeeded/SILVA/silva132_97.fna -strand plus -mode high_confidence -notmatched 10_UchimeReference/AllSamples_unoise_nc.fasta")

#	OTU table generation
# An OTU table is made by the otutab command
dir.create("11_OtuTable")

system("./usearch11.0.667 -otutab 07_JoinedData/AllSamples.fasta -db 10_UchimeReference/AllSamples_unoise_nc.fasta -id 0.97 -otutabout 11_OtuTable/AllSamples_unoise_otu_table1.txt")
#system("./usearch11.0.66 -otutab 07_JoinedData/AllSamples.fasta -db 10_UchimeReference/AllSamples_unoise_nc.fasta -id 0.97 -otutabout 11_OtuTable/AllSamples_unoise_otu_table2.txt -maxaccepts 0 -maxrejects 0") #very slow but more accurate


##############################################################################################################
####################	Mapping of OTUs on Reference Database with best-blast-match approach ###################
##############################################################################################################

# dir.create("12_TaxAssignment")
# 
# system("makeblastdb -in 00_DataNeeded/SILVA/silva132_97.fna -dbtype nucl")
# system("blastn -query 10_UchimeReference/AllSamples_unoise_nc.fasta -outfmt 6 -out 12_TaxAssignment/AllSamples_unoise_nc.txt -db 00_DataNeeded/SILVA/silva132_97.fna -evalue 1e-20 -num_threads 2")
# 
# blast_tab = read.delim(file = "12_TaxAssignment/AllSamples_unoise_nc.txt", header = F)
# best_hit = NULL
# for(id in unique(blast_tab$V1))
# {
#   best_hit = c(best_hit, which(blast_tab$V1 == id)[1])
# }
# blast_tab_red = blast_tab[best_hit,]
# write.table(x = blast_tab_red, file = "12_TaxAssignment/AllSamples_unoise_nc2.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = F)
# 
# rm(list=ls(all=TRUE))

##############################################################################################################
########################### Mapping of OTUs on Reference Database with BLCA  (start) #########################
##############################################################################################################

# Follow instructions in section "prepare GTDB SSU database" in the link below to prepare db files for BLCA 
# https://github.com/songweizhi/Katana_cmds/blob/master/BLCA_GTDB.sh
# copy the generated ssu_all_r95_BLCAparsed.taxonomy and ssu_all_r95_BLCAparsed.fasta to 00_DataNeeded/BLCA/db_GTDB_SSU
# 00_DataNeeded/BLCA/db_GTDB_SSU/ssu_all_r95_BLCAparsed.taxonomy
# 00_DataNeeded/BLCA/db_GTDB_SSU/ssu_all_r95_BLCAparsed.fasta

library(tidyverse)
library(reshape2)

# add clustalo and muscle to system path (required by BLCA)
system("cp 00_DataNeeded/clustalo/clustal-omega-1.2.3-macosx 00_DataNeeded/clustalo/clustalo ; chmod +x 00_DataNeeded/clustalo/clustalo") # MacOS
system("cp 00_DataNeeded/muscle/muscle3.8.31_i86darwin64 00_DataNeeded/muscle/muscle ; chmod +x 00_DataNeeded/muscle/muscle") # MacOS
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "00_DataNeeded/clustalo", sep = ":"))
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "00_DataNeeded/muscle", sep = ":"))

dir.create("12_TaxAssignmentGTDB_BLCA")

# run BLCA agsinst GTDB SSU database

# If you saw a "BioPython is not detected!" error, decomment the following line to install biopython
system(paste(system("which python3", intern = TRUE), '-m pip install biopython --user'))

# r95
system('python3 00_DataNeeded/BLCA/2.blca_main.py -r 00_DataNeeded/BLCA/db_GTDB_SSU/ssu_all_r95_BLCAparsed.taxonomy -q 00_DataNeeded/BLCA/db_GTDB_SSU/ssu_all_r95_BLCAparsed.fasta -i 10_UchimeReferenceGTDB/zOTUs_nc.fasta -o 12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.1.txt')

# remove taxonomic ranks from BLCA's classification and put number in parenthesis (make it easy to read)
m = 1
for (each_line in readLines(file('12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out1.txt',open="r")) ){
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
  {cat(taxonomy_no_rank_with_OTU,file='12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.2.txt',sep="\n", append=FALSE)} 
  else 
  {cat(taxonomy_no_rank_with_OTU,file='12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.2.txt',sep="\n", append=TRUE)}
  m = m +1
}

##############################################################################################################
######################### Mapping of OTUs on Reference Database with BLCA (end) ##############################
##############################################################################################################


#Merge Table and Taxonomy
dir.create("13_FinalOtuTable")

SLV = read.delim("00_DataNeeded/SILVA/raw_taxonomy.txt", header = F)
names(SLV) = c("SilvaId", "taxonomy")

OTU = read.delim("11_OtuTable/AllSamples_unoise_otu_table1.txt", header = T)
TAX = read.delim("12_TaxAssignment/AllSamples_unoise_nc2.txt", header = F)

TAX = TAX[,c(1:4,11)]
names(TAX) = c("X.OTU.ID", "SilvaId", "Identity", "Alignment_Length", "Evalue")

TAX_SLV = merge(TAX, SLV, by = "SilvaId")
OTU_TAX = merge(OTU, TAX_SLV, by = "X.OTU.ID")

write.table(OTU_TAX, "13_FinalOtuTable/AllSamples_unoise_otu_table.txt", sep = "\t", row.names = F, col.names = T, quote = F)

rm(list=ls(all=TRUE))
