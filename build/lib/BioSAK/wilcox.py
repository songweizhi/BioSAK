import argparse
import numpy as np
from scipy.stats import wilcoxon


wilcox_usage = '''
===================== wilcox example commands ===================== 

BioSAK wilcox -1 paired_treatment_1.txt -2 paired_treatment_2.txt

===================================================================
'''


def wilcox(args):

    s1_txt = args['1']
    s2_txt = args['2']

    s1_list = []
    for e1 in open(s1_txt):
        s1_list.append(float(e1.strip()))

    s2_list = []
    for e2 in open(s2_txt):
        s2_list.append(float(e2.strip()))

    s1_array = np.array(s1_list)
    s2_array = np.array(s2_list)
    stats_result = wilcoxon(s1_array, s2_array)

    print(stats_result)


if __name__ == '__main__':

    wilcox_parser = argparse.ArgumentParser()
    wilcox_parser.add_argument('-1', required=True, help='treatment 1')
    wilcox_parser.add_argument('-2', required=True, help='treatment 2')
    args = vars(wilcox_parser.parse_args())
    wilcox(args)
