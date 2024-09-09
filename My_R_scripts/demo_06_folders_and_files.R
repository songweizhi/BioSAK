
file_list = dir(path = '/Users/songweizhi/Desktop/ccc', pattern = '.fasta')

print(file_list)


library(tools)

pwd_file = '/Users/songweizhi/Desktop/grouping_out.purity.txt'

file_path = dirname(pwd_file)
file_path
# "/Users/songweizhi/Desktop"

file_path2 = dirname('grouping_out.purity.txt')
file_path2
# "."

file_name = basename(pwd_file)
file_name
# "grouping_out.purity.txt"

file_name_without_extension = file_path_sans_ext(file_name)
file_name_without_extension
# "grouping_out.purity"

file_name_without_extension2 = file_path_sans_ext(basename(pwd_file))
file_name_without_extension2
# "grouping_out.purity"

file_extension = file_ext(file_name)
file_extension
# "txt"

file_extension2 = file_ext(pwd_file)
file_extension2
# "txt"

# remove path from bin_subfolder
bin_subfolder_only_name = strsplit(bin_subfolder, '/')[[1]][-1]
