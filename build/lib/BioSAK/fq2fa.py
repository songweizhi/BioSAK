import os
import argparse


fq2fa_usage = '''
========== fq2fa example commands ==========

BioSAK fq2fa -i reads_R1.fq -o reads_R1.fa
BioSAK fq2fa -i reads_R2.fq -o reads_R2.fa

============================================
'''


def fq2fa(args):
    fq_in   = args['i']
    fa_out  = args['o']
    sed_cmd = "sed -n '1~4s/^@/>/p;2~4p' %s > %s" % (fq_in, fa_out)
    os.system(sed_cmd)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, help='input fastq')
    parser.add_argument('-o', required=True, help='output fasta')
    args = vars(parser.parse_args())
    fq2fa(args)
