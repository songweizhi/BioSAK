import argparse


DendroDF_usage = '''
================== DendroDF example commands ==================

BioSAK DendroDF -i matrix.txt -o Dendrogram.pdf

===============================================================
'''


def DendroDF(args):

    input_matrix_txt = args['i']
    output_pdf_txt   = args['o']


if __name__ == '__main__':

    DendroDF_parser = argparse.ArgumentParser(usage=DendroDF_usage)
    DendroDF_parser.add_argument('-i', required=True, help='input data matrix')
    DendroDF_parser.add_argument('-o', required=True, help='output plot, in pdf format')
    args = vars(DendroDF_parser.parse_args())
    DendroDF(args)
