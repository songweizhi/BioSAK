import os
import argparse


PhyloBiAssoc_usage = '''
=========== PhyloBiAssoc example commands ===========

BioSAK PhyloBiAssoc -t demo.tre -d demo.txt

# Note: the header for the first two columns has to be "ID" and "cate".

# Output example
# gene_1  phylosig  7.973475e-26  binaryPGLMM  0.03255813
# gene_2  phylosig  1             chisq.test   0.7183411
# gene_3  phylosig  2.66378e-08   binaryPGLMM  0.5969282
# gene_4  phylosig  7.169588e-08  binaryPGLMM  0.999338

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
    PhyloBiAssoc_parser = argparse.ArgumentParser()
    PhyloBiAssoc_parser.add_argument('-t', required=True, help='tree file')
    PhyloBiAssoc_parser.add_argument('-d', required=True, help='data file')
    args = vars(PhyloBiAssoc_parser.parse_args())
    PhyloBiAssoc(args)
