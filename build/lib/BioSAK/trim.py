import os
import argparse


trim_usage = '''
========================== trim example commands ==========================

BioSAK trim -qci -qco -a TruSeq3-PE-2.fa -1 R1.fastq -2 R2.fastq -x fastq
BioSAK trim -qci -qco -a TruSeq3-PE-2.fa -1 R1.fq.gz -2 R2.fq.gz -x fq.gz

===========================================================================
'''


def trim(args):

    fastq_r1        = args['1']
    fastq_r2        = args['2']
    file_ext        = args['x']
    adapter_file    = args['a']
    leading         = args['leading']
    trailing        = args['trailing']
    crop            = args['crop']
    headcrop        = args['headcrop']
    swl             = args['swl']
    swq             = args['swq']
    minlen          = args['minlen']
    qc_fastq_in     = args['qci']
    qc_fastq_out    = args['qco']

    r1_base = fastq_r1[:-(len(file_ext) + 1)]
    r2_base = fastq_r2[:-(len(file_ext) + 1)]

    r1_p  = '%s_P.%s'   % (r1_base, file_ext)
    r1_up = '%s_UP.%s'  % (r1_base, file_ext)
    r2_p  = '%s_P.%s'   % (r2_base, file_ext)
    r2_up = '%s_UP.%s'  % (r2_base, file_ext)

    qc_cmd_in   = 'fastqc %s %s' % (fastq_r1, fastq_r2)
    trim_cmd    = 'trimmomatic PE %s %s %s %s %s %s ILLUMINACLIP:%s:2:30:10 LEADING:%s TRAILING:%s CROP:%s HEADCROP:%s SLIDINGWINDOW:%s:%s MINLEN:%s' % (fastq_r1, fastq_r2, r1_p, r1_up, r2_p, r2_up, adapter_file, leading, trailing, crop, headcrop, swl, swq, minlen)
    qc_cmd_out  = 'fastqc %s %s' % (r1_p, r2_p)

    # run trimmomatic
    print(trim_cmd)
    os.system(trim_cmd)

    # run fastqc
    if qc_fastq_in is True:
        print(qc_cmd_in)
        os.system(qc_cmd_in)

    if qc_fastq_out is True:
        print(qc_cmd_out)
        os.system(qc_cmd_out)

    print('Done!')


if __name__ == '__main__':
    trim_parser = argparse.ArgumentParser(usage=trim_usage)
    trim_parser.add_argument('-1',              required=True,                          help='fastq R1')
    trim_parser.add_argument('-2',              required=True,                          help='fastq R2')
    trim_parser.add_argument('-x',              required=True,                          help='file extension, e.g., fastq or fastq.gz')
    trim_parser.add_argument('-a',              required=True,                          help='adapter file, e.g., TruSeq3-PE-2.fa')
    trim_parser.add_argument('-leading',        required=False, type=int, default=25,   help='leading, default is 25')
    trim_parser.add_argument('-trailing',       required=False, type=int, default=25,   help='trailing, default is 25')
    trim_parser.add_argument('-crop',           required=False, type=int, default=145,  help='crop, default is 145')
    trim_parser.add_argument('-headcrop',       required=False, type=int, default=5,    help='headcrop, default is 5')
    trim_parser.add_argument('-swl',            required=False, type=int, default=5,    help='slidingwindow length, default is 5')
    trim_parser.add_argument('-swq',            required=False, type=int, default=25,   help='slidingwindow q-value, default is 25')
    trim_parser.add_argument('-minlen',         required=False, type=int, default=36,   help='minlen, default is 36')
    trim_parser.add_argument('-qci',            required=False, action="store_true",    help='specify to run fastqc for input files')
    trim_parser.add_argument('-qco',            required=False, action="store_true",    help='specify to run fastqc for output files')
    args = vars(trim_parser.parse_args())
    trim(args)
