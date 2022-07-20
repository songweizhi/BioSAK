import os
import argparse

GTDB_tree_r207_usage = '''
============== GTDB_tree_r207 example command ==============

module load python/3.7.3
source ~/mypython3env/bin/activate
module load perl/5.28.0
module load prodigal/2.6.3
module load pplacer/1.1.alpha19
module load hmmer/3.3
module load fasttree/2.1.11
module load gcc/6.2.0
module load gsl/2.6
module load fastani/1.32
module load R/3.5.3
export GTDBTK_DATA_PATH=/srv/scratch/z5039045/DB/GTDB_r207/release207

BioSAK GTDB_tree_r207 -p Demo -i gnm_folder -x fa -t 12

============================================================
'''

def GTDB_tree_r207(args):

    input_gnm_dir       = args['i']
    output_prefix       = args['p']
    file_extension      = args['x']
    num_threads         = args['t']

    output_dir          = '%s_GTDB_tree_r207_wd'                                                    % output_prefix
    msa_bac120_gz       = '%s/align/gtdbtk.bac120.user_msa.fasta.gz'                                % output_dir
    msa_bac120          = '%s/align/gtdbtk.bac120.user_msa.fasta'                                   % output_dir
    msa_ar53_gz         = '%s/align/gtdbtk.ar53.user_msa.fasta.gz'                                  % output_dir
    msa_ar53            = '%s/align/gtdbtk.ar53.user_msa.fasta'                                     % output_dir

    cmd_identify        = 'gtdbtk identify --genome_dir %s -x %s --out_dir %s --cpus %s'            % (input_gnm_dir, file_extension, output_dir, num_threads)
    cmd_align           = 'gtdbtk align --identify_dir %s --out_dir %s --cpus %s'                   % (output_dir, output_dir, num_threads)
    cmd_gunzip_bac120   = 'gunzip %s'                                                               % msa_bac120_gz
    cmd_gunzip_ar53     = 'gunzip %s'                                                               % msa_ar53_gz
    cmd_infer_bac120    = 'gtdbtk infer --msa_file %s --out_dir %s --cpus %s --prefix %s_bac120'    % (msa_bac120, output_dir, num_threads, output_prefix)
    cmd_infer_ar53      = 'gtdbtk infer --msa_file %s --out_dir %s --cpus %s --prefix %s_ar53'      % (msa_ar53, output_dir, num_threads, output_prefix)

    print(cmd_identify)
    os.system(cmd_identify)
    print(cmd_align)
    os.system(cmd_align)

    if os.path.isfile(msa_bac120_gz):
        print(cmd_gunzip_bac120)
        os.system(cmd_gunzip_bac120)
        print(cmd_infer_bac120)
        os.system(cmd_infer_bac120)

    if os.path.isfile(msa_ar53_gz):
        print(cmd_gunzip_ar53)
        os.system(cmd_gunzip_ar53)
        print(cmd_infer_ar53)
        os.system(cmd_infer_ar53)

    inferred_bac120_tree = '%s/%s_bac120.unrooted.tree' % (output_dir, output_prefix)
    inferred_ar53_tree   = '%s/%s_ar53.unrooted.tree'   % (output_dir, output_prefix)

    if os.path.isfile(inferred_bac120_tree):
        print('Inferred bacterial tree:\t%s' % inferred_bac120_tree)
    if os.path.isfile(inferred_ar53_tree):
        print('Inferred archaeal tree:\t%s' % inferred_ar53_tree)

    print('Done!')


if __name__ == '__main__':

    GTDB_tree_r207_parser = argparse.ArgumentParser(usage=GTDB_tree_r207_usage)
    GTDB_tree_r207_parser.add_argument('-p',               required=True,                            help='output prefix')
    GTDB_tree_r207_parser.add_argument('-i',               required=True,                            help='genome folder')
    GTDB_tree_r207_parser.add_argument('-x',               required=True,                            help='genome file extension')
    GTDB_tree_r207_parser.add_argument('-t',               required=False, type=int, default=1,      help='number of threads')
    args = vars(GTDB_tree_r207_parser.parse_args())
    GTDB_tree_r207(args)


'''
Test2_ar53.unrooted.tree
Test2_bac120.unrooted.tree
'''