import os
import argparse
from time import sleep


get_abd2_mapping_usage = '''
================================= get_abd2_mapping example commands =================================

python3 get_abd2_mapping.py -r masked_ref.fa -t 36 -p Cliona1 -1 Cliona1_1.fastq -2 Cliona1_2.fastq
python3 get_abd2_mapping.py -r masked_ref.fa -t 36 -p Cliona2 -1 Cliona2_1.fastq -2 Cliona2_2.fastq
python3 get_abd2_mapping.py -r masked_ref.fa -t 36 -p Cliona3 -1 Cliona3_1.fastq -2 Cliona3_2.fastq
python3 get_abd2_mapping.py -r masked_ref.fa -t 36 -p Cliona4 -1 Cliona4_1.fastq -2 Cliona4_2.fastq

# Output files from this step: cov, rpkm and stat files

=====================================================================================================
'''


def mapping(args):

    fq_r1             = args['1']
    fq_r2             = args['2']
    ref_seq           = args['r']
    op_prefix         = args['p']
    num_threads       = args['t']

    bwa_cmd           = 'bwa mem -5SP -t %s %s %s %s | samblaster > %s.sam'                                                                                        % (num_threads, ref_seq, fq_r1, fq_r2, op_prefix)
    samtools_view_cmd = 'samtools view -@ 32 -bS -h -b %s.sam > %s.bam'                                                                                            % (op_prefix, op_prefix)
    samtools_sort_cmd = 'samtools sort -@ 32 %s.bam -o %s.sorted.bam'                                                                                              % (op_prefix, op_prefix)
    coverm_filter_cmd = 'coverm filter -b %s.sorted.bam --min-read-aligned-percent 0.9 --min-read-percent-identity 0.99 --output-bam-files %s.sorted_filtered.bam' % (op_prefix, op_prefix)
    pileup_sh_cmd     = 'pileup.sh in=%s.sorted_filtered.bam out=%s.sorted_filtered.cov rpkm=%s.sorted_filtered.rpkm overwrite=true'                               % (op_prefix, op_prefix, op_prefix)
    seqkit_stat_cmd   = 'seqkit stat %s > %s.stat'                                                                                                                 % (fq_r1, op_prefix)

    print(bwa_cmd)
    os.system(bwa_cmd)

    print(samtools_view_cmd)
    os.system(samtools_view_cmd)
    sleep(1)
    os.system('rm %s.sam' % op_prefix)

    print(samtools_sort_cmd)
    os.system(samtools_sort_cmd)
    sleep(1)
    os.system('rm %s.bam' % op_prefix)

    print(coverm_filter_cmd)
    os.system(coverm_filter_cmd)
    sleep(1)
    os.system('rm %s.sorted.bam' % op_prefix)

    print(pileup_sh_cmd)
    os.system(pileup_sh_cmd)
    sleep(1)
    os.system('rm %s.sorted_filtered.bam' % op_prefix)

    print(seqkit_stat_cmd)
    os.system(seqkit_stat_cmd)


if __name__ == '__main__':

    abund_parser = argparse.ArgumentParser()
    abund_parser.add_argument('-1', required=True,                          help='fastq r1')
    abund_parser.add_argument('-2', required=True,                          help='fastq r2')
    abund_parser.add_argument('-r', required=True,                          help='masked reference sequence')
    abund_parser.add_argument('-p', required=True,                          help='output prefix')
    abund_parser.add_argument('-t', required=False, type=int, default=1,    help='number of threads')
    args = vars(abund_parser.parse_args())
    mapping(args)
