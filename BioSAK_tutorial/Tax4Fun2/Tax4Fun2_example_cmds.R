
# Tax4Fun2 was developed by Dr Bernd Wemheuer, this is only a short script showing how it can be installed and executed for functional prediction.

#######################################################################################################################
############################################# Section 1: install Tax4Fun2 #############################################
#######################################################################################################################

# source code of Tax4Fun2 (Tax4Fun2_1.1.5.tar.gz) is available here:
# https://github.com/songweizhi/BioSAK/tree/master/BioSAK_tutorial/Tax4Fun2

# There are two approaches to install Tax4Fun2 on you computer

# Approach 1:
# open RStudio > Tools > install packages > from "Install from" select "Package Archive File" > choose Tax4Fun2_1.1.5.tar.gz > Install

# Approach 2:
# Copy the source code Tax4Fun2_1.1.5.tar.gz to the folder where you want Tax4Fun2 to be installed.
# Here we are going to install it to the "Software" folder under my home directory
# Forward to the "Software" folder
cd /Users/songweizhi/Software
# paste the following two lines in your terminal to install Tax4Fun2
R
install.packages(pkgs = "Tax4Fun2_1.1.5.tar.gz", repos = NULL, source = TRUE)


#######################################################################################################################
############################################# Section 2: before you start #############################################
#######################################################################################################################

# You'll need to have Tax4Fun2's default database (Tax4Fun2_ReferenceData_v2) to perform functional prediction
# It's can be found in torsten/Weizhi_Song/Software_Scripts_Databases/Tax4Fun2/Tax4Fun2_ReferenceData_v2.tar.gz in Torsten's storage

# There are two ways of running Tax4Fun2
# refers to section 3 if you want to make functional predictions with Tax4Fun2's default database
# refers to section 4 if you want to make functional predictions with Tax4Fun2's default database together with a user-generated database


#######################################################################################################################
######################### Section 3: Making functional predictions using the default database #########################
#######################################################################################################################

library(Tax4Fun2)

# modify the following 8 lines
pwd_op_folder   = 'path/to/Tax4Fun2_output_folder'      # specify output folder
query_otu_seq   = 'path/to/demo_OTU.fasta'              # OTU sequence file
query_otu_table = 'path/to/demo_OTU_table.csv'          # OTU table
pwd_ref_data    = 'path/to/Tax4Fun2_ReferenceData_v2'   # path to Tax4Fun2's default database (i.e., Tax4Fun2_ReferenceData_v2), need to be decompressed before use
norm_by_cn      = TRUE                                  # normalize_by_copy_number (TRUE or FALSE)
norm_path       = TRUE                                  # normalize_pathways (TRUE or FALSE)
iden            = 0.97                                  # min_identity_to_reference
num_of_threads  = 6                                     # number of CPU cores to use 

# Run Tax4Fun2
runRefBlast(path_to_otus = query_otu_seq, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", use_force = T, num_threads = num_of_threads)
makeFunctionalPrediction(path_to_otu_table = query_otu_table, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", normalize_by_copy_number = norm_by_cn, min_identity_to_reference = iden, normalize_pathways = norm_path)


#######################################################################################################################
########## Section 4: Making functional predictions using the default database and a user-generated database ##########
#######################################################################################################################

# Please refers to steps in section 4.1 to generate your own database.
# You'll need to store your reference genomes/MAGs in a folder and provide full path to that folder to Tax4Fun2 with "pwd_user_data".
# You also need to specify a name to the to-be generated database with "name_user_data".
# !!! Reference genomes/MAGs without 16S rRNA genes will be ignored in database building.
# !!! You need to know the copy number of 16S rRNA genes in your reference genomes/MAGs, MarkerMAG might be able to help with this.

# I have created a user database with MAGs derived from sponge-associated microbial community.
# located at torsten/Weizhi_Song/Software_Scripts_Databases/Tax4Fun2/Sponge_associated_MAGs_with_added_16S.tar.gz
# To use it, replace the corresponding lines in section 4.2 with
# pwd_user_data  = 'path/to/Sponge_associated_MAGs_with_added_16S'
# name_user_data = 'Sponge_associated_MAGs_Tax4Fun2_db'

# Tahsin Khan has created a user database with Ecklonia radiata-associated MAGs,
# located at torsten/Weizhi_Song/Software_Scripts_Databases/Tax4Fun2/Ulva_MAGs_w16S.tar.gz
# To use it, replace the corresponding lines in section 4.2 with
# pwd_user_data  = 'path/to/Ulva_MAGs_w16S'
# name_user_data = 'Ulva_MAGs_db'

# Database files need to be decompressed before use
# tar -xzvf Sponge_associated_MAGs_with_added_16S.tar.gz
# tar -xzvf Ulva_MAGs_w16S.tar.gz

library(Tax4Fun2)

################## 4.1 use the following lines to Generate your own reference database ##################

# modify the following 4 lines
pwd_ref_data    = 'path/to/Tax4Fun2_ReferenceData_v2'   # path to Tax4Fun2's default database, need to be decompressed before use
pwd_user_data   = 'path/to/your/genome_folder'          # path to the folder holds the reference genome/MAG files
name_user_data  = 'name_of_user_database'               # specify the name of the generated database, specify only the name, do not include path here!
gnm_ext         = 'fna'                                 # extension of the genome files, normally in fasta, fna or fa

# Generate your own database
extractSSU(genome_folder = pwd_user_data, file_extension = gnm_ext, path_to_reference_data = pwd_ref_data)
assignFunction(genome_folder = pwd_user_data, file_extension = gnm_ext, path_to_reference_data = pwd_ref_data, num_of_threads = num_threads, fast = TRUE)
generateUserData(path_to_reference_data = pwd_ref_data, path_to_user_data = pwd_user_data, name_of_user_data = name_user_data, SSU_file_extension = "_16SrRNA.ffn", KEGG_file_extension = "_funPro.txt")

###### 4.2 making functional predictions using the default database and the user-generated database ######

# modify the following 10 lines
pwd_op_folder   = 'path/to/Tax4Fun2_output_folder'      # specify output folder
query_otu_seq   = 'path/to/demo_OTU.fasta'              # OTU sequence file
query_otu_table = 'path/to/demo_OTU_table.csv'          # OTU table
pwd_ref_data    = 'path/to/Tax4Fun2_ReferenceData_v2'   # need to be the same as in 4.1
pwd_user_data   = 'path/to/your/genome_folder'          # need to be the same as in 4.1
name_user_data  = 'name_of_user_database'               # need to be the same as in 4.1 , specify only the name, do not include path here
norm_by_cn      = TRUE                                  # normalize_by_copy_number (TRUE or FALSE)
norm_path       = TRUE                                  # normalize_pathways (TRUE or FALSE)
iden            = 0.97                                  # min_identity_to_reference
num_of_threads  = 6                                     # number of CPU cores to use 

# Run Tax4Fun2 with your data included
runRefBlast(path_to_otus = query_otu_seq, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", use_force = T, num_threads = num_of_threads, include_user_data = T, path_to_user_data = pwd_user_data, name_of_user_data = name_user_data)
makeFunctionalPrediction(path_to_otu_table = query_otu_table, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", normalize_by_copy_number = norm_by_cn, min_identity_to_reference = iden, normalize_pathways = norm_path, include_user_data = T, path_to_user_data = pwd_user_data, name_of_user_data = name_user_data)

