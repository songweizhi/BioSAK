import os
import argparse


abd_16s_usage = '''
======================== abd_16s example commands ========================

BioSAK abd_16s -i 16S.fa -r1 R1.fa -r2 R2.fa -o op_dir -p demo -t 12 -f 

==========================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def abd_16s(args):

    input_16s_seq           = args['i']
    reads_r1_fasta          = args['r1']
    reads_r2_fasta          = args['r2']
    op_dir                  = args['o']
    op_prefix               = args['p']
    num_threads             = args['t']
    force_overwrite         = args['f']

    # create op dir
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output directory detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    ########################################################################################################################

    _, _, input_16s_base, _ = sep_path_basename_ext(input_16s_seq)

    bowtie_build_cmd        = 'bowtie2-build --quiet --threads %s -f %s %s'                                                                                                     % (num_threads, input_16s_seq, input_16s_base)
    bowtie_read_to_16s_cmd  = 'bowtie2 -x %s -U %s,%s -S %s/%s.sam -p %s -f --xeq --local --all --no-unal -N 1 -L 30 --quiet'                                                   % (input_16s_base, reads_r1_fasta, reads_r2_fasta, op_dir, op_prefix, num_threads)
    samtools_view_cmd       = 'samtools view -@ 32 -bS -h -b %s/%s.sam > %s/%s.bam'                                                                                             % (op_dir, op_prefix, op_dir, op_prefix)
    samtools_sort_cmd       = 'samtools sort -@ 32 %s/%s.bam -o %s/%s.sorted.bam'                                                                                               % (op_dir, op_prefix, op_dir, op_prefix)
    coverm_filter_cmd       = 'coverm filter -b %s/%s.sorted.bam --min-read-aligned-percent 0.9 --min-read-percent-identity 0.99 --output-bam-files %s/%s.sorted_filtered.bam'  % (op_dir, op_prefix, op_dir, op_prefix)
    pileup_sh_cmd           = 'pileup.sh in=%s/%s.sorted_filtered.bam out=%s/%s.sorted_filtered.cov rpkm=%s/%s.sorted_filtered.rpkm overwrite=true'                             % (op_dir, op_prefix, op_dir, op_prefix, op_dir, op_prefix)
    seqkit_stat_cmd         = 'seqkit stat %s > %s/%s.stat'                                                                                                                     % (reads_r1_fasta, op_dir, op_prefix)

    os.system(bowtie_build_cmd)
    os.system(bowtie_read_to_16s_cmd)
    os.system(samtools_view_cmd)
    os.system(samtools_sort_cmd)
    os.system(coverm_filter_cmd)
    os.system(pileup_sh_cmd)
    os.system(seqkit_stat_cmd)


if __name__ == '__main__':

    abd_16s_parser = argparse.ArgumentParser()
    abd_16s_parser.add_argument('-i',   required=True,                          help='16S sequences')
    abd_16s_parser.add_argument('-r1',  required=True,                          help='reads r1')
    abd_16s_parser.add_argument('-r2',  required=True,                          help='reads r2')
    abd_16s_parser.add_argument('-o',   required=True,                          help='output directory')
    abd_16s_parser.add_argument('-p',   required=True,                          help='output prefix')
    abd_16s_parser.add_argument('-t',   required=False, type=int, default=1,    help='number of threads, default is 1')
    abd_16s_parser.add_argument('-f',   required=False, action="store_true",    help='force overwrite')
    args = vars(abd_16s_parser.parse_args())
    abd_16s(args)
