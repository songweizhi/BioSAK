import os

text_file = '/Users/songweizhi/Desktop/Sponge_r226/sponge_phylogeny/18S/RefSeq_18S_wd/tmp/AJ633926.1.gbk'
text_file_size = os.stat(text_file).st_size
print(text_file_size)