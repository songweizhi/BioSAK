import os
import glob
import pandas
import argparse
import numpy as np


ko00001_keg = '/Users/songweizhi/DB/KEGG_DB/ko00001.keg'


for each_line in open(ko00001_keg):
    if each_line[0] in ['A', 'B', 'C']:
        print(each_line.strip())


remaining_ABC = '''


'''
