import os
import argparse


VisBlastOp_usage = '''
====================== VisBlastOp example commands ======================

BioSAK VisBlastOp -h

=========================================================================
'''


def VisBlastOp(args):
    print('To be added!')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-gene', required=True, help='gene id')
    args = vars(parser.parse_args())
    VisBlastOp(args)
