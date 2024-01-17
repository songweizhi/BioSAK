import os
import argparse


PhyloBiAssoc_usage = '''
=========== PhyloBiAssoc example commands ===========

BioSAK PhyloBiAssoc -t demo.tre -d demo.txt

# Note: the header for the first two columns has to be "ID" and "cate".

# It will perform:
# 1) binaryPGLMM test if phylosig p-value <= 0.05 (significant phylogenetic signal)
# 2) chi-squared test if phylosig p-value > 0.05  (no phylogenetic signal)
# 3) do nothing if phylosig returns NaN (might due to the same value across all genomes)

# https://www.rdocumentation.org/packages/ape/versions/5.7-1/topics/binaryPGLMM

=====================================================
'''


def PhyloBiAssoc(args):

    tree_file        = args['t']
    data_file        = args['d']
    pwd_current_script  = os.path.realpath(__file__)
    current_script_path = '/'.join(pwd_current_script.split('/')[:-1])
    PhyloBiAssoc_R   = '%s/PhyloBiAssoc.R' % current_script_path
    PhyloBiAssoc_cmd = 'Rscript %s -t %s -d %s' % (PhyloBiAssoc_R, tree_file, data_file)
    os.system(PhyloBiAssoc_cmd)


if __name__ == "__main__":
    PhyloBiAssoc_parser = argparse.ArgumentParser(usage=PhyloBiAssoc_usage)
    PhyloBiAssoc_parser.add_argument('-t', required=True, help='tree file')
    PhyloBiAssoc_parser.add_argument('-d', required=True, help='data file')
    args = vars(PhyloBiAssoc_parser.parse_args())
    PhyloBiAssoc(args)
