
#######################################################################################################################
############################################# Section 1: install Tax4Fun2 #############################################
#######################################################################################################################

# there are two ways to install Tax4Fun2 on you computer

# option 1:
# open RStudio > Tools > install packages > select "Package Archive File" from "Install from" > choose Tax4Fun2_1.1.5.tar.gz > Install

# option 2:
# Copy the source code Tax4Fun2_1.1.5.tar.gz to the folder where you want Tax4Fun2 to be installed.
# Here we are going to install it to the "Software" folder under my home directory
# Forward to the "Software" folder
cd /Users/songweizhi/Software
# paste the following two lines in your terminal to install Tax4Fun2
R
install.packages(pkgs = "Tax4Fun2_1.1.5.tar.gz", repos = NULL, source = TRUE)

# you'll need to have Tax4Fun2's default database (Tax4Fun2_ReferenceData_v2) to perform functional prediction


#######################################################################################################################
######################### Section 2: Making functional predictions using the default database #########################
#######################################################################################################################

library(Tax4Fun2)

# modify the following 8 lines
pwd_op_folder   = 'path/to/Tax4Fun2_outputs'            # specify output folder
query_otu_seq   = 'path/to/your/OTU_sequence.fa'        # OTU sequence file
query_otu_table = 'path/to/your/OTU_table.txt'          # OTU table
pwd_ref_data    = 'path/to/Tax4Fun2_ReferenceData_v2'   # path to Tax4Fun2's default database (i.e., Tax4Fun2_ReferenceData_v2), need to be decompressed before use
norm_by_cn      = TRUE                                  # normalize_by_copy_number (TRUE or FALSE)
norm_path       = TRUE                                  # normalize_pathways (TRUE or FALSE)
iden            = 0.97                                  # min_identity_to_reference
num_of_threads  = 6                                     # number of CPU cores to use 

# Run the reference blast
runRefBlast(path_to_otus = query_otu_seq, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", use_force = T, num_threads = num_of_threads)

# Predicting functional profiles
makeFunctionalPrediction(path_to_otu_table = query_otu_table, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", normalize_by_copy_number = norm_by_cn, min_identity_to_reference = iden, normalize_pathways = norm_path)


#######################################################################################################################
########## Section 3: Making functional predictions using the default database and a user-generated database ##########
#######################################################################################################################

library(Tax4Fun2)

################## 3.1 use the following lines to Generate your own reference database ##################

# modify the following 4 lines
pwd_ref_data    = 'path/to/Tax4Fun2_ReferenceData_v2'   # path to Tax4Fun2's default database, need to be decompressed before use
#pwd_user_data   = 'path/to/Ulva_MAGs_w16S'             # path to the folder holds the genome files 
#name_user_data  = 'Ulva_MAGs_db'                       # specify the name of the generated database, specify only the name, do not include path here!
pwd_user_data   = 'path/to/Ulva_MAGs_w16S'              # path to the folder holds the genome files 
name_user_data  = 'Ulva_MAGs_db'                        # specify the name of the generated database, specify only the name, do not include path here!
gnm_ext         = 'fna'                                 # extension of the genome files, normally in fasta, fna or fa

extractSSU(genome_folder = pwd_user_data, file_extension = gnm_ext, path_to_reference_data = pwd_ref_data)
assignFunction(genome_folder = pwd_user_data, file_extension = gnm_ext, path_to_reference_data = pwd_ref_data, num_of_threads = num_threads, fast = TRUE)
generateUserData(path_to_reference_data = pwd_ref_data, path_to_user_data = pwd_user_data, name_of_user_data = name_user_data, SSU_file_extension = "_16SrRNA.ffn", KEGG_file_extension = "_funPro.txt")

###### 3.2 making functional predictions using the default database and the user-generated database ######

# modify the following 10 lines
pwd_op_folder   = 'path/to/Tax4Fun2_outputs'            # specify output folder
query_otu_seq   = 'path/to/your/OTU_sequence.fa'        # OTU sequence file
query_otu_table = 'path/to/your/OTU_table.txt'          # OTU table
pwd_ref_data    = 'path/to/Tax4Fun2_ReferenceData_v2'   # need to be the same as in 3.1
#pwd_user_data   = 'path/to/Ulva_MAGs_w16S'             # need to be the same as in 3.1
#name_user_data  = 'Ulva_MAGs_db'                       # need to be the same as in 3.1 , specify only the name, do not include path here 
pwd_user_data   = 'path/to/Ulva_MAGs_w16S'              # need to be the same as in 3.1
name_user_data  = 'Ulva_MAGs_db'                        # need to be the same as in 3.1 , specify only the name, do not include path here 
norm_by_cn      = TRUE                                  # normalize_by_copy_number (TRUE or FALSE)
norm_path       = TRUE                                  # normalize_pathways (TRUE or FALSE)
iden            = 0.97                                  # min_identity_to_reference
num_of_threads  = 6                                     # number of CPU cores to use 

# 1. Run the reference blast with include_user_data = TRUE and specifiy the path to the user data
runRefBlast(path_to_otus = query_otu_seq, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", use_force = T, num_threads = num_of_threads, include_user_data = T, path_to_user_data = pwd_user_data, name_of_user_data = name_user_data)

# 2. Making the prediction with your data included
makeFunctionalPrediction(path_to_otu_table = query_otu_table, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", normalize_by_copy_number = norm_by_cn, min_identity_to_reference = iden, normalize_pathways = norm_path, include_user_data = T, path_to_user_data = pwd_user_data, name_of_user_data = name_user_data)

