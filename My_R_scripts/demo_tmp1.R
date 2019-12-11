library('seqinr')
library('rlist')


setwd('/Users/songweizhi/Desktop/Binning_refiner_R_wd2')
input_bin_folder = 'input_bin_folder'
input_bin_folder


file_out = 'contigs2bin.tab'

# get input_bin_subfolders
input_bin_subfolders = list.dirs(input_bin_folder, recursive=FALSE)
input_bin_subfolders


# get contig to bin information
file_out_handle = file(file_out, open = "w")

# initialize a dataframe
ctg_to_bin_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(ctg_to_bin_df)  = c("MetaBAT", "MyCC", "CONCOCT")


# read in contig to bin information
for (bin_subfolder in input_bin_subfolders){
  bin_file_list = list.files(bin_subfolder)
  for (each_bin in bin_file_list){
    pwd_each_bin = paste(bin_subfolder, each_bin, sep = '/')
    for (each_seq in read.fasta(file = pwd_each_bin, seqtype = 'DNA', as.string = 1, forceDNAtolower = 0)){
      
      # get sequence id
      each_seq_id = getName(each_seq)
      
      # remove path from bin_subfolder
      bin_subfolder_only_name = strsplit(bin_subfolder, '/')[[1]][-1]
      
      # turn what to write out into dataframe
      for_write_df = data.frame(each_seq_id, bin_subfolder_only_name, each_bin)
      
      # write out into file
      write.table(for_write_df, file_out_handle, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
      
    }
  }
}

# close file
close(file_out_handle)

ctg_to_bin_df


#read in contig to bin file 

# contig2bin_matrix = read.csv(file_out)
# 
# print(contig2bin_matrix)
# 




# ##########################################################################################
# #                                        FOR TEST                                        #
# ##########################################################################################
# 
# for_test = list(a=c(1,'b',3), b=NULL)
# for_test
# 
# if (exists('c', where=for_test) == TRUE){
#   for_test$a = c(for_test$a, 555)
# } else {
#   for_test = list.append(for_test, c=c('haha'))
# }
# 
# for_test
# 
# ##########################################################################################
# ##########################################################################################
# 
# 
# 
# #bin_folders = list.dirs(path = ".", all.files = FALSE)
# bin_folders = c('MetaBAT', 'MyCC')
# 
# file_out = file('bin_contigs.tab', open = "w")
# 
# for (each_folder in bin_folders){
#   
# #  # create tmp folder
# #  each_folder_tmp = paste(each_folder, '_tmp', sep = '')
# #  if(file.exists(each_folder_tmp)){
# #    unlink(each_folder_tmp, recursive = TRUE)
# #    dir.create(each_folder_tmp)
# #  }else{
# #    dir.create(each_folder_tmp)
# #  }
#   
#   # get bin file list
#   current_bins = list.files(each_folder)
#   
#   # read in sequences
#   for(each_bin in current_bins){
#     pwd_each_bin = paste(each_folder, each_bin, sep = '/')
#     each_bin_seqs = read.fasta(file = pwd_each_bin, seqtype = 'DNA', as.string = 1, forceDNAtolower = 0)
#     current_bin_seq_list = c(each_bin)
#     print(typeof(current_bin_seq_list))
#     for (each_seq in each_bin_seqs){
#       seq_id = getName(each_seq)
#       current_bin_seq_list = c(current_bin_seq_list, seq_id)
#     }
#     print(current_bin_seq_list)
#     print(typeof(current_bin_seq_list))
#     
#     to_write = paste(toString(current_bin_seq_list), sep = ',')
#     print(to_write)
#     
#     write.csv(to_write, file = file_out)
#     
#   }
# }
# 
# close(file_out)
# 
# 
# 
# 
# # if (exists('each_seq_id', where=seq_lol) == FALSE) {
# #   seq_lol = list.append(seq_lol, each_seq_id=c(each_bin))
# # } else {
# #   seq_lol$each_seq_id = c(seq_lol$each_seq_id, each_bin)
# # }