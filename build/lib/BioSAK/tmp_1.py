
import os

dir = '/Users/songweizhi/Desktop/demo'
sub_dir_list = next(os.walk(dir))[1]
print(sub_dir_list)
