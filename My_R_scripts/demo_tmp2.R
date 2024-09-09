library(assertr) # to Concatenate all columns of each row in data frame


bin_subfolder = 'home/Desktop'
bin_subfolder
bin_subfolder_split = strsplit(bin_subfolder, '/')
bin_subfolder_split

needed = bin_subfolder_split[[1]][-1]
needed
typeof(bin_subfolder_split)



m0 <- matrix(NA, 4, 0)
rownames(m0)
m0

m2 <- cbind(1, 1:4)
colnames(m2, do.NULL = FALSE)
colnames(m2) <- c("x","Y")
rownames(m2) <- rownames(m2, do.NULL = FALSE, prefix = "Obs.")
m2


my_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(my_df)  = c("MetaBAT", "MyCC", "CONCOCT")
my_df['lima', 'MyCC'] = 'aa'
my_df['song', 'MetaBAT'] = 'bb'
my_df['lima', 'MetaBAT'] = 'cc'
my_df['lima', 'CONCOCT'] = 3
my_df['bins', 'MyCC'] = 'bb'
my_df['song', 'CONCOCT'] = 1
my_df['bins', 'CONCOCT'] = 10
my_df

my_df_reordered = my_df[order(my_df$CONCOCT),]
my_df_reordered






my_df$na_count <- apply(my_df, 1, function(x) sum(is.na(x)))
my_df

my_df$concate = ifelse(apply(my_df, 1, function(x) sum(is.na(x))) == 0, col_concat(my_df, sep = '___'), NA)
my_df


my_df$concate = ifelse(is.na(my_df) == TRUE, col_concat(my_df, sep = '___'), NA)
my_df


my_df$concate = ifelse(any(is.na(my_df) == TRUE) == TRUE, col_concat(my_df, sep = '___'), NA)
my_df


my_df$concate = col_concat(my_df, sep = '___')
my_df


if (any(is.na(my_df) == TRUE) == TRUE){
  my_df$concate= col_concat(my_df, sep = '___')
} else {
  my_df$concate = NA
}

my_df$concate[any(is.na(my_df) == TRUE) == TRUE] = col_concat(my_df, sep = '___')
my_df





my_df$song

any(my_df$song == <NA>)

?any


my_str = c(1,2,3,4)

any(my_df['lima',] == 1)


any(is.na(c(1,2,NA)))

    
    
any(is.na(my_df['song',]) == TRUE) == TRUE

df <- data.frame(B=c("A","B","C","C"), C=c("A","C","B","B"), D=c("B","A","C","A") )
df
df$A<-ifelse(df$B==df$C,paste(df$D),paste(df$C))

return(df)




