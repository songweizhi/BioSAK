
######################## get file list ########################

import glob

file_dir = '/Users/songweizhi/Desktop/DateArTree'
file_ext = 'txt'

file_re = '%s/*.%s' % (file_dir, file_ext)

file_list = glob.glob(file_re)
file_list = [os.path.basename(i) for i in glob.glob(file_re)]


print(file_list)


####################### get folder list #######################

import os

dir = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/op_qualified_OGs_gene_tree_dir'
sub_dir_list = next(os.walk(dir))[1]
print(sub_dir_list)


############# get_time_since_last_modification #############

import os
import time

def get_time_since_last_modification(target_folder_or_file):

    last_modified_time = os.path.getmtime(target_folder_or_file)
    current_time = time.time()

    tslm_sec  = current_time - last_modified_time
    tslm_min  = tslm_sec / 60
    tslm_hour = tslm_sec / (60 * 60)
    tslm_day  = tslm_sec / (60 * 60 * 24)

    tslm_sec  = float("{0:.2f}".format(tslm_sec))
    tslm_min  = float("{0:.2f}".format(tslm_min))
    tslm_hour = float("{0:.2f}".format(tslm_hour))
    tslm_day  = float("{0:.2f}".format(tslm_day))

    return tslm_sec, tslm_min, tslm_hour, tslm_day


target_file = '/Users/songweizhi/Desktop/foo.txt'
tslm_sec, tslm_min, tslm_hour, tslm_day = get_time_since_last_modification(target_file)


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

