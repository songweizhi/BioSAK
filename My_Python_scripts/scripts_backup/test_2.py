import os
import glob
from Bio import SeqIO

file_re = '/Users/songweizhi/Desktop/GoodBins_0.5_0.05/*.fasta'
#file_re = '/srv/scratch/z5039045/get_bin_depth_wd/GoodBins_0.5_0.05/*.fasta'
file_list = [os.path.basename(file_name) for file_name in glob.glob(file_re)]
file_list = glob.glob(file_re)

n = 0
for each_file in file_list:

    for each_seq in SeqIO.parse(each_file, 'fasta'):
        n += 1
print(n)


print(file_list)
