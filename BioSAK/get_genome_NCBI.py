import os
import glob
import argparse
import multiprocessing as mp
from distutils.spawn import find_executable


get_genome_NCBI_parser_usage = '''
===================================== get_genome_NCBI example commands =====================================

# Dependency: datasets (https://www.ncbi.nlm.nih.gov/datasets/docs/v1/download-and-install/)

# example command
BioSAK get_genome_NCBI -i GCA_000006945
BioSAK get_genome_NCBI -i GCA_009840555.1
BioSAK get_genome_NCBI -i genome_id.txt -o downloaded_gnm_dir -t 8 -f

# Format of genome_id.txt
GCA_009837245
GCA_009840555.1
GCA_000006945.2

# You can use this link to get the id of genomes from a specific taxon , IDs are in the "Assembly" column.
https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/refseq_category:reference

============================================================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def datasets_worker(arg_list):

    gnm_id                      = arg_list[0]
    op_dir                      = arg_list[1]
    downloaded_file             = '%s/%s.zip'                                                                               % (op_dir, gnm_id)
    downloaded_file_unzip_dir   = '%s/%s'                                                                                   % (op_dir, gnm_id)
    datasets_cmd                = 'datasets download genome accession %s --include genome --no-progressbar --filename %s'   % (gnm_id, downloaded_file)
    os.system(datasets_cmd)

    # datasets download genome accession %s --include genome --no-progressbar --filename %s
    # datasets download genome accession GCA_014843695.1 --include gff3,rna,cds,protein,genome,seq-report

    if os.path.isfile(downloaded_file) is True:

        if os.path.isdir(downloaded_file_unzip_dir) is True:
            os.system('rm -r %s' % downloaded_file_unzip_dir)

        unzip_cmd = 'unzip -q %s -d %s' % (downloaded_file, downloaded_file_unzip_dir)
        os.system(unzip_cmd)

        fna_file_re  = '%s/ncbi_dataset/data/*.*/*.*_genomic.fna' % (downloaded_file_unzip_dir)
        fna_file_list = glob.glob(fna_file_re)

        if len(fna_file_list) == 1:
            pwd_fna_file_unzipped = fna_file_list[0]

            _, f_base, _ = sep_path_basename_ext(pwd_fna_file_unzipped)
            pwd_fna_file = '%s/%s.fna' % (op_dir, '_'.join(f_base.split('_')[:2]))

            mv_cmd = 'mv %s %s'  % (pwd_fna_file_unzipped, pwd_fna_file)
            os.system(mv_cmd)
            os.system('rm %s'    % downloaded_file)
            os.system('rm -r %s' % downloaded_file_unzip_dir)
        else:
            print('Download failed: %s!' % gnm_id)
    else:
        print('Download failed: %s!' % gnm_id)


def download_GenBank_genome(args):

    assembly_id     = args['i']
    output_folder   = args['o']
    num_threads     = args['t']
    force_overwrite = args['f']

    if find_executable('datasets') is None:
        print('datasets not detected, please install it first')
        print('get help with: BioSAK get_genome_NCBI -h')
        exit()

    if os.path.isfile(assembly_id) is False:
        print('Download single genome mode')
        if output_folder is None:
            output_folder = '.'
        datasets_worker([assembly_id, output_folder])

    else:
        # check output dir
        if output_folder is None:
            print('An output directory need to be provided, program exited!')
            exit()
        else:
            if os.path.isdir(output_folder) is True:
                if force_overwrite is True:
                    os.system('rm -r %s' % output_folder)
                else:
                    print('Specified output directory detected, program exited!')
                    exit()
            os.mkdir(output_folder)

        report_file = '%s_report.txt'        % output_folder
        failed_file = '%s_report_failed.txt' % output_folder

        datasets_worker_arg_lol = []
        for each_gnm in open(assembly_id):
            gnm_id = each_gnm.strip()
            datasets_worker_arg_lol.append([gnm_id, output_folder])

        # run with multiprocessing
        pool = mp.Pool(processes=num_threads)
        pool.map(datasets_worker, datasets_worker_arg_lol)
        pool.close()
        pool.join()

        # report
        failed_gnm_set = set()
        id_to_file_dict = dict()
        report_file_handle = open(report_file, 'w')
        for each_gnm in open(assembly_id):
            gnm_id = each_gnm.strip()
            current_gnm_re = '%s/%s*.fna' % (output_folder, gnm_id)
            current_gnm_list = glob.glob(current_gnm_re)

            if len(current_gnm_list) == 1:
                pwd_current_gnm = current_gnm_list[0]
                id_to_file_dict[gnm_id] = pwd_current_gnm.split('/')[-1]
                report_file_handle.write('%s\t%s\n' % (gnm_id, pwd_current_gnm.split('/')[-1]))
            else:
                failed_gnm_set.add(gnm_id)
        report_file_handle.close()

        # write out failed files
        if len(failed_gnm_set) > 0:
            failed_file_handle = open(failed_file, 'w')
            failed_file_handle.write('\n'.join(failed_gnm_set))
            failed_file_handle.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage=get_genome_NCBI_parser_usage)
    parser.add_argument('-i',   required=True,                       help='IDs of genomes to download')
    parser.add_argument('-o',   required=False, default=None,        help='output folder')
    parser.add_argument('-t',   required=False, default=1, type=int, help='number of threads')
    parser.add_argument('-f',   required=False, action="store_true", help='force overwrite existing genome directory')
    args = vars(parser.parse_args())
    download_GenBank_genome(args)
