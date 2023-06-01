import os
import argparse
from BioSAK.BioSAK_config import config_dict


PhyloBiAssoc_usage = '''
=========== PhyloBiAssoc example commands ===========

BioSAK PhyloBiAssoc -tree demo.tre -data demo.txt

=====================================================
'''


def PhyloBiAssoc(args):

    tree_file       = args['tree']
    data_file       = args['data']
    PhyloBiAssoc_R  = config_dict['PhyloBiAssoc_R']

    PhyloBiAssoc_cmd = 'Rscript %s -t %s -d %s' % (PhyloBiAssoc_R, tree_file, data_file)
    os.system(PhyloBiAssoc_cmd)


if __name__ == "__main__":
    PhyloBiAssoc_parser = argparse.ArgumentParser()
    PhyloBiAssoc_parser.add_argument('-tree', required=True, help='tree file')
    PhyloBiAssoc_parser.add_argument('-data', required=True, help='data file')
    args = vars(PhyloBiAssoc_parser.parse_args())
    PhyloBiAssoc(args)
