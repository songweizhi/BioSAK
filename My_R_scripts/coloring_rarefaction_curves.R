Author:"Marwan Majzoub"
Description:"This file contains all the commands for further analysis of the OTU table"
Date:"30/05/2021"

library("ggplot2")
library("vegan")
library("reshape")
library("gridExtra")
library("tidyverse")
library("dplyr")
library(RColorBrewer)
library(picante)
library(seqinr)
library(drc)
library(lattice)
library(permute)


# Read otu table
AllSamples <- read.delim("coloring_rarefaction_curves_otu_table.txt")
#another way to read the OTU table
#otu_table = t(read.delim('otu_table.txt', row.names = 1))

#RArefaction curve
#OTU Table only
otu_table <- dplyr::select(AllSamples, -OTU_ID)
otu_table <- as.matrix(t(otu_table))
colnames(otu_table) <- AllSamples[,1]

#sample with lowest no of seqs
raremax <- min(rowSums(otu_table))

####################################################################################

out <- rarecurve(otu_table, step=50, label=T,ylab = "OTU Count")

rare <- lapply(out, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

names(rare) <- rownames(otu_table)

rare <- map_dfr(rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

ggplot(data = rare)+
  geom_line(aes(x = raw.read, y = OTU, color = sample)) +
  geom_text(data = rare %>% #here we need coordinates of the labels
              group_by(sample) %>% #first group by samples
              summarise(max_OTU = max(OTU), #find max OTU
                        max_raw = max(raw.read)), #find max raw read
            aes(x = max_raw, y = max_OTU, label = sample), check_overlap = T, hjust = 0)+
  scale_x_continuous(labels =  scales::scientific_format()) +
  labs(colour = "black", x = "Sample Size", y = "OTU Count")

####################################################################################


