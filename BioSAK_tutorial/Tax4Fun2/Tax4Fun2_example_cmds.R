
################################################################################
############################### install Tax4Fun2 ###############################
################################################################################

# Copy the source code Tax4Fun2_1.1.5.tar.gz to the folder where you want Tax4Fun2 to be installed.
# Here we are going to install it to the "Software" folder under my home directory
# Forward to the "Software" folder
cd /Users/songweizhi/Software

# option 1: paste the following two lines in your terminal to install Tax4Fun2
R
install.packages(pkgs = "Tax4Fun2_1.1.5.tar.gz", repos = NULL, source = TRUE)

# option 2:
open RStudio > Tools > install packages > select "Package Archive File" from "Install from" > Install

# you'll need to have Tax4Fun2's default database (Tax4Fun2_ReferenceData_v2) to perform functional prediction
# it is in the folder: torsten/Weizhi_Song/Tax4Fun2/Tax4Fun2_ReferenceData_v2

################################################################################
######## Making functional predictions using Tax4Fun2's default database #######
################################################################################

library(Tax4Fun2)

# modify the following 8 lines
pwd_op_folder   = '/Users/songweizhi/Desktop/Tax4Fun2_outputs'  # specify output folder
query_otu_seq   = 'path/to/your_OTU_sequence.fa'  # OTU sequence file
query_otu_table = 'path/to/your_OTU_table.txt'  # OTU table
norm_by_cn      = TRUE  # normalize_by_copy_number (TRUE or FALSE)
norm_path       = TRUE  # normalize_pathways (TRUE or FALSE)
iden            = 0.97  # min_identity_to_reference
num_of_threads  = 6     # number of CPU cores to use 
pwd_ref_data    = '/Users/songweizhi/Desktop/Tax4Fun2_wd/batch_0/Tax4Fun2_ReferenceData_v2'  # path to Tax4Fun2's default database

# Run the reference blast
runRefBlast(path_to_otus = query_otu_seq, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", use_force = T, num_threads = num_of_threads)

# Predicting functional profiles
makeFunctionalPrediction(path_to_otu_table = query_otu_table, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", normalize_by_copy_number = norm_by_cn, min_identity_to_reference = iden, normalize_pathways = norm_path)


############################################################################################################
########## Making functional predictions using the default database and a user-generated database ##########
############################################################################################################

################## first, use the following lines to Generate your own reference database ##################

library(Tax4Fun2)

# modify the following 4 lines
pwd_ref_data    = '/Users/songweizhi/Tax4Fun2_ReferenceData_v2'  # path to Tax4Fun2's default database
pwd_user_data   = '/Users/songweizhi/Sponge_associated_MAGs_with_16S'  # path to the folder holds the genome files 
name_user_data  = 'Sponge_associated_Tax4Fun2_db'  # specify the name of the to-be generated database
gnm_ext         = 'fna'  # extension of the genome files, normally in fasta, fna or fa

# 1. Extracting SSU sequences from genomes
extractSSU(genome_folder = pwd_user_data, file_extension = gnm_ext, path_to_reference_data = pwd_ref_data)

# 2. Assigning function to genomes
assignFunction(genome_folder = pwd_user_data, file_extension = gnm_ext, path_to_reference_data = pwd_ref_data, num_of_threads = num_threads, fast = TRUE)

# 3. Generate user-defined reference data
generateUserData(path_to_reference_data = pwd_ref_data, path_to_user_data = pwd_user_data, name_of_user_data = name_user_data, SSU_file_extension = "_16SrRNA.ffn", KEGG_file_extension = "_funPro.txt")


###### then, making functional predictions using the default database and the user-generated database ######

# modify the following 10 lines
pwd_op_folder   = '/Users/songweizhi/Desktop/Tax4Fun2_outputs'  # specify output folder
query_otu_seq   = 'path/to/your_OTU_sequence.fa'  # OTU sequence file
query_otu_table = 'path/to/your_OTU_table.txt'  # OTU table
norm_by_cn      = TRUE  # normalize_by_copy_number (TRUE or FALSE)
norm_path       = TRUE  # normalize_pathways (TRUE or FALSE)
iden            = 0.97  # min_identity_to_reference
num_of_threads  = 6   # number of CPU cores to use 
pwd_ref_data    = '/Users/songweizhi/Tax4Fun2_ReferenceData_v2'  # need to be the same as above
pwd_user_data   = '/Users/songweizhi/Sponge_associated_MAGs_with_16S'  # need to be the same as above
name_user_data  = 'Sponge_associated_Tax4Fun2_db'  # need to be the same as above

# 1. Run the reference blast with include_user_data = TRUE and specifiy the path to the user data
runRefBlast(path_to_otus = query_otu_seq, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", use_force = T, num_threads = num_of_threads, include_user_data = T, path_to_user_data = pwd_user_data, name_of_user_data = name_user_data)

# 2. Making the prediction with your data included
makeFunctionalPrediction(path_to_otu_table = query_otu_table, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", normalize_by_copy_number = norm_by_cn, min_identity_to_reference = iden, normalize_pathways = norm_path, include_user_data = T, path_to_user_data = pwd_user_data, name_of_user_data = name_user_data)

