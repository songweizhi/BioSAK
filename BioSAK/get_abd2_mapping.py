import os
import argparse
from time import sleep


get_abd2_mapping_usage = '''
================================= get_abd2_mapping example commands =================================

BioSAK get_abd2_mapping -r masked_ref.fa -t 36 -p Cliona1 -1 Cliona1_1.fastq -2 Cliona1_2.fastq
BioSAK get_abd2_mapping -r masked_ref.fa -t 36 -p Cliona2 -1 Cliona2_1.fastq -2 Cliona2_2.fastq
BioSAK get_abd2_mapping -r masked_ref.fa -t 36 -p Cliona3 -1 Cliona3_1.fastq -2 Cliona3_2.fastq
BioSAK get_abd2_mapping -r masked_ref.fa -t 36 -p Cliona4 -1 Cliona4_1.fastq -2 Cliona4_2.fastq

# Output files from this step: cov, rpkm and stat files

=====================================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def get_abd2_mapping(args):

    fq_r1             = args['1']
    fq_r2             = args['2']
    ref_seq           = args['r']
    op_dir            = args['o']
    op_prefix         = args['p']
    num_threads       = args['t']
    force_overwrite   = args['f']

    # check input files
    if os.path.isfile(ref_seq) is False:
        print('%s not found, program exited!' % ref_seq)
        exit()
    if os.path.isfile(fq_r1) is False:
        print('%s not found, program exited!' % fq_r1)
        exit()
    if os.path.isfile(fq_r2) is False:
        print('%s not found, program exited!' % fq_r2)
        exit()

    # create op dir
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output directory detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    fq1_name, _, _, _ = sep_path_basename_ext(fq_r1)
    fq2_name, _, _, _ = sep_path_basename_ext(fq_r2)
    ref_name, _, _, _ = sep_path_basename_ext(ref_seq)

    perform_decompress = False
    fq_r1_decompressed = fq_r1
    fq_r2_decompressed = fq_r2
    if fq_r1.endswith('.gz'):
        fq_r1_decompressed  = '%s/%s'               % (op_dir, fq1_name[:-3])
        fq_r2_decompressed  = '%s/%s'               % (op_dir, fq2_name[:-3])
        gunzip_cmd_r1       = 'gunzip -c %s > %s'   % (fq_r1, fq_r1_decompressed)
        gunzip_cmd_r2       = 'gunzip -c %s > %s'   % (fq_r2, fq_r2_decompressed)
        perform_decompress  = True
        os.system(gunzip_cmd_r1)
        os.system(gunzip_cmd_r2)

    elif fq_r1.endswith('.tar.gz'):
        tar_cmd_r1 = ''
        tar_cmd_r2 = ''
        print('tar to be added, program exited!')
        exit()

    # index reference sequences
    index_ref_cmd     = 'cp %s %s/; bwa index %s/%s'                                                                                                                        % (ref_seq, op_dir, op_dir, ref_name)
    bwa_cmd           = 'bwa mem -5SP -t %s %s/%s %s %s | samblaster > %s/%s.sam'                                                                                           % (num_threads, op_dir, ref_name, fq_r1_decompressed, fq_r2_decompressed, op_dir, op_prefix)
    samtools_view_cmd = 'samtools view -@ 32 -bS -h -b %s/%s.sam > %s/%s.bam'                                                                                               % (op_dir, op_prefix, op_dir, op_prefix)
    samtools_sort_cmd = 'samtools sort -@ 32 %s/%s.bam -o %s/%s.sorted.bam'                                                                                                 % (op_dir, op_prefix, op_dir, op_prefix)
    coverm_filter_cmd = 'coverm filter -b %s/%s.sorted.bam --min-read-aligned-percent 0.9 --min-read-percent-identity 0.99 --output-bam-files %s/%s.sorted_filtered.bam'    % (op_dir, op_prefix, op_dir, op_prefix)
    pileup_sh_cmd     = 'pileup.sh in=%s/%s.sorted_filtered.bam out=%s/%s.sorted_filtered.cov rpkm=%s/%s.sorted_filtered.rpkm overwrite=true'                               % (op_dir, op_prefix, op_dir, op_prefix, op_dir, op_prefix)
    seqkit_stat_cmd   = 'seqkit stat %s > %s/%s.stat'                                                                                                                       % (fq_r1_decompressed, op_dir, op_prefix)

    print(seqkit_stat_cmd)
    os.system(seqkit_stat_cmd)

    print(index_ref_cmd)
    os.system(index_ref_cmd)

    print(bwa_cmd)
    os.system(bwa_cmd)
    if perform_decompress is True:
        os.system('rm %s/%s' % (op_dir, fq1_name[:-3]))
        os.system('rm %s/%s' % (op_dir, fq2_name[:-3]))
    os.system('rm %s/%s'   % (op_dir, ref_name))
    os.system('rm %s/%s.*' % (op_dir, ref_name))

    print(samtools_view_cmd)
    os.system(samtools_view_cmd)
    sleep(1)
    os.system('rm %s/%s.sam' % (op_dir, op_prefix))

    print(samtools_sort_cmd)
    os.system(samtools_sort_cmd)
    sleep(1)
    os.system('rm %s/%s.bam' % (op_dir, op_prefix))

    print(coverm_filter_cmd)
    os.system(coverm_filter_cmd)
    sleep(1)
    os.system('rm %s/%s.sorted.bam' % (op_dir, op_prefix))

    print(pileup_sh_cmd)
    os.system(pileup_sh_cmd)
    sleep(1)
    os.system('rm %s/%s.sorted_filtered.bam' % (op_dir, op_prefix))

    print('Done!')


if __name__ == '__main__':

    abund_parser = argparse.ArgumentParser()
    abund_parser.add_argument('-1', required=True,                          help='fastq r1')
    abund_parser.add_argument('-2', required=True,                          help='fastq r2')
    abund_parser.add_argument('-r', required=True,                          help='masked reference sequence')
    abund_parser.add_argument('-o', required=True,                          help='output directory')
    abund_parser.add_argument('-p', required=True,                          help='output prefix')
    abund_parser.add_argument('-t', required=False, type=int, default=1,    help='number of threads, default is 1')
    abund_parser.add_argument('-f', required=False, action="store_true",    help='force overwrite')
    args = vars(abund_parser.parse_args())
    get_abd2_mapping(args)
