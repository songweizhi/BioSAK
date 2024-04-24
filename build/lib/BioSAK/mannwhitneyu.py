import argparse
from scipy.stats import mannwhitneyu


mannwhitneyu_usage = '''
============== mannwhitneyu example commands ============== 

BioSAK mannwhitneyu -1 treatment_1.txt -2 treatment_2.txt

===========================================================
'''


def mannwhitneyu(args):
    s1_txt = args['1']
    s2_txt = args['2']

    s1_list = []
    for e1 in open(s1_txt):
        s1_list.append(float(e1.strip()))

    s2_list = []
    for e2 in open(s2_txt):
        s2_list.append(float(e2.strip()))

    _, p_value = mannwhitneyu(s1_list, s2_list)

    print('P value: %s' % p_value)


if __name__ == '__main__':

    mannwhitneyu_parser = argparse.ArgumentParser()
    mannwhitneyu_parser.add_argument('-1', required=True, help='sample 1')
    mannwhitneyu_parser.add_argument('-2', required=True, help='sample 2')
    args = vars(mannwhitneyu_parser.parse_args())
    mannwhitneyu(args)
