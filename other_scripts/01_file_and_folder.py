
########################################################################################################################
######################################################### file #########################################################
########################################################################################################################

#################### sep_path_basename_ext ####################

import os

def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


file_1 = 'input_sequence.fasta'
file_2 = '/Users/weizhisong/Softwares/artemis/act.txt'
file_path_1, file_basename_1, file_extension_1 = sep_path_basename_ext(file_1)
file_path_2, file_basename_2, file_extension_2 = sep_path_basename_ext(file_2)

print(file_path_1)          # .
print(file_basename_1)      # input_sequence
print(file_extension_1)     # .fasta
print(file_path_2)          # /Users/weizhisong/Softwares/artemis
print(file_basename_2)      # act
print(file_extension_2)     # .txt


######################## get file size ########################

text_file = '/Users/songweizhi/Desktop/a.txt'
text_file_size = os.stat(text_file).st_size


###############################################################

# Sometimes when we only care about certain parts of a tuple when unpacking, a dummy variable like _ is used as placeholder
___, filename = os.path.split('/home/luciano/.ssh/idrsa.pub')
print(filename)     # idrsa.pub


########################################################################################################################
######################################################## folder ########################################################
########################################################################################################################

##################### force_create_folder #####################

import os
import shutil

def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


################## get_no_hidden_folder_list ##################

import os
import glob

def get_no_hidden_folder_list(wd):
    folder_list = []
    for each_folder in os.listdir(wd):
        if not each_folder.startswith('.'):
            folder_list.append(each_folder)
    return folder_list


####################### get folder list #######################

wd = '/Users/songweizhi/Desktop/test_'
print(get_no_hidden_folder_list(wd))


######################## get file list ########################

file_re = '/Users/songweizhi/Desktop/test/*.fa'
file_list = [os.path.basename(file_name) for file_name in glob.glob(file_re)]
print(file_list)


#################### whether folder exist #####################

import os

# check whether folder exist
folder_1 = '/Users/songweizhi/Desktop/folder_1'
folder_2 = '/Users/songweizhi/Desktop/folder_2'
file_1 = '/Users/songweizhi/Desktop/folder_1/file.txt'
os.mkdir(folder_2)

print(os.path.isdir(folder_1))  # False
print(os.path.isdir(folder_2))  # True
print(os.path.isfile(file_1))   # False


#################### remove folders or files #####################

def rm_folder_file(target_re):

    target_list = glob.glob(target_re)

    for target in target_list:

        if os.path.isdir(target) is True:
            os.system('rm -r %s' % target)
            
        elif os.path.isfile(target) is True:
            os.system('rm %s' % target)

