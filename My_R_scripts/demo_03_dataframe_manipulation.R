library(dplyr)


#################### How To Create A Simple Data Frame in R? ####################

# make some vectors first
Died.At <- c(22,40,72,41)
Writer.At <- c(16, 18, 36, 36)
First.Name <- c("John", "Edgar", "Walt", "Jane")
Second.Name <- c("Doe", "Poe", "Whitman", "Austen")
Sex <- c("MALE", "MALE", "MALE", "FEMALE")
Date.Of.Death <- c("2015-05-10", "1849-10-07", "1892-03-26","1817-07-18")

# combine the vectors with the data.frame()
writers_df <- data.frame(Died.At, Writer.At, First.Name, Second.Name, Sex, Date.Of.Death)

# Use str() to get more info about the generated data frame
str(writers_df)
writers_df


################################### add row name ###################################

setwd('/Users/songweizhi/R_scripts')
csv_in ='Demo_03_ESBL.csv'

# read in csv
AR_matrix = read.csv(csv_in)
AR_matrix

row.names(AR_matrix) = paste('Sample', seq.int(nrow(AR_matrix)), sep = "_")
AR_matrix


############################ subset df to get metadata ############################

AR_matrix_metadata = AR_matrix[,1:4]
AR_matrix_metadata


########################## subset df to get value matrix ##########################

AR_matrix_value = AR_matrix[,5:ncol(AR_matrix)]
AR_matrix_value


###################################################################################

rm(list=ls())

