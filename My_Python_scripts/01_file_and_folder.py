
######################## get file list ########################

import glob

file_re = '/Users/songweizhi/Desktop/test/*.fa'
file_list = glob.glob(file_re)
print(file_list)


####################### get folder list #######################

import os

dir = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/op_qualified_OGs_gene_tree_dir'
sub_dir_list = next(os.walk(dir))[1]
print(sub_dir_list)


#################### sep_path_basename_ext ####################

import os

def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


file_demo = '/Users/weizhisong/Softwares/artemis/act.txt'
file_path, file_basename, file_ext = sep_path_basename_ext(file_demo)


#################### whether folder/file exist #####################

import os

folder_1 = '/Users/songweizhi/Desktop/folder_1'
file_1 = '/Users/songweizhi/Desktop/folder_1/file.txt'

print(os.path.isdir(folder_1))
print(os.path.isfile(file_1))


######################## get file size ########################

text_file = '/Users/songweizhi/Desktop/a.txt'
text_file_size = os.stat(text_file).st_size


###############################################################

# Sometimes when we only care about certain parts of a tuple when unpacking, a dummy variable like _ is used as placeholder
___, filename = os.path.split('/home/luciano/.ssh/idrsa.pub')
print(filename)     # idrsa.pub


#################### remove folders or files #####################

def rm_folder_file(target_re):

    target_list = glob.glob(target_re)

    for target in target_list:

        if os.path.isdir(target) is True:
            os.system('rm -r %s' % target)
            
        elif os.path.isfile(target) is True:
            os.system('rm %s' % target)

