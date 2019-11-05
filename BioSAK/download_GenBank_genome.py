import os
import shutil
import argparse
from datetime import datetime
import multiprocessing as mp


download_GenBank_genome_parser_usage = '''
================================== dwnld_GenBank_genome example commands ==================================

# Usage:
# 1. Go to https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/
# 2. Search genomes you want to download (e.g. prokaryotes, proteobacteria or psychrobacter)
# 3. Click "Download" on the right side
# 4. provide the downloaded csv file with '-csv'

# Download all genomes in file prokaryotes.csv
BioSAK dwnld_GenBank_genome -csv prokaryotes.csv

# Only download genomes with provided accessions in accessions.txt
BioSAK dwnld_GenBank_genome -csv prokaryotes.csv -id accessions.txt

# accessions.txt file format 
# one accession per line, accessions can be found from the 6th column of the prokaryotes.csv file.
GCA_000265385.1
GCA_900176135
GCA_900107535.1

===========================================================================================================
'''


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def genome_download_worker(argument_list):

    genome_record_split = argument_list[0]
    downloaded_genome_folder = argument_list[1]
    GenBank_FTP = genome_record_split[-2][1:-1]
    GenBank_FTP_id = GenBank_FTP.strip().split('/')[-1]
    assembly_id = genome_record_split[5]

    # download
    wget_cmd = 'wget %s/%s_genomic.fna.gz -P %s -q' % (GenBank_FTP, GenBank_FTP_id, downloaded_genome_folder)
    os.system(wget_cmd)

    # decompress
    os.system('gunzip %s/%s_genomic.fna.gz' % (downloaded_genome_folder, GenBank_FTP_id))

    # rename
    os.system('mv %s/%s_genomic.fna %s/%s.fna' % (downloaded_genome_folder, GenBank_FTP_id, downloaded_genome_folder, assembly_id))


def download_GenBank_genome(args):

    assembly_id_file = args['id']
    csv_file = args['csv']
    num_threads = args['t']
    time_format = '[%Y-%m-%d %H:%M:%S] '

    in_file_path, in_file_basename, in_file_extension = sep_path_basename_ext(csv_file)
    downloaded_genome_folder = '%s_genomes' % in_file_basename
    force_create_folder(downloaded_genome_folder)

    genomes_to_download = []
    if (assembly_id_file is not None) and (os.path.isfile(assembly_id_file) is True):

        # get id list to download
        for each_id in open(assembly_id_file):
            each_id = each_id.strip()
            each_id_no_version = each_id
            if '.' in each_id:
                each_id_no_version = each_id.split('.')[0]
            genomes_to_download.append(each_id_no_version)

    # report
    if (assembly_id_file is not None) and (os.path.isfile(assembly_id_file) is True):
        print(datetime.now().strftime(time_format) + 'Downloading %s genomes with %s cores' % (len(genomes_to_download), num_threads))
    else:
        print(datetime.now().strftime(time_format) + 'Downloading genomes with %s cores' % (num_threads))

    # download genome with multiprocessing
    list_for_multiple_arguments_download_worker = []
    for genome_record in open(csv_file):

        if not genome_record.startswith('#Organism Name'):
            genome_record_split = genome_record.strip().split(',')
            assembly_id = genome_record_split[5]
            assembly_id_no_version = assembly_id.split('.')[0][1:]

            if (assembly_id_file is not None) and (os.path.isfile(assembly_id_file) is True):
                if assembly_id_no_version in genomes_to_download:
                    list_for_multiple_arguments_download_worker.append([genome_record_split, downloaded_genome_folder])
            else:
                list_for_multiple_arguments_download_worker.append([genome_record_split, downloaded_genome_folder])

    # run COG annotaion files with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(genome_download_worker, list_for_multiple_arguments_download_worker)
    pool.close()
    pool.join()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-csv', required=True,                          help='csv file from NCBI genome_browse')
    parser.add_argument('-id',  required=False, default=None,           help='assembly accessions, one per line')
    parser.add_argument('-t',   required=False, default=1, type=int,    help='number of threads')

    args = vars(parser.parse_args())

    download_GenBank_genome(args)
