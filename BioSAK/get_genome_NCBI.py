import os
import glob
import argparse
from datetime import datetime
import multiprocessing as mp


get_genome_NCBI_parser_usage = '''
===================================== get_genome_NCBI example commands =====================================

# Usage:
Before you start, you need to download https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt,
and provide it to get_genome_NCBI with '-db'.

# Download all genomes in prokaryotes.csv
BioSAK get_genome_NCBI -db prokaryotes.csv -out gnm_dir -fna -name -t 8

# Only download genomes in gnm_to_download.txt
BioSAK get_genome_NCBI -db prokaryotes.csv -out gnm_dir -fna -name -t 8 -id gnm_to_download.txt

# Format of genome id file (specify with "-id")
# Genome ID can be found in the "Assembly Accession" column of the prokaryotes.txt, one id per line
GCA_009840555.1
GCA_009840575.1

# You may want to get the id of genomes from a specific taxon using:
https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/refseq_category:reference

============================================================================================================
'''


def genome_download_worker(argument_list):

    genome_record_split =      argument_list[0]
    col_index =                argument_list[1]
    downloaded_genome_folder = argument_list[2]
    get_fna =                  argument_list[3]
    get_faa =                  argument_list[4]
    get_gbff =                 argument_list[5]
    with_name =                argument_list[6]

    genome_name = genome_record_split[0]
    genome_name_no_space = '_'.join(genome_name.split(' '))
    if ('(' in genome_name_no_space) and (')' in genome_name_no_space):
        genome_name_no_space_no_parenthesis_tmp = genome_name_no_space.replace('(', '_')
        genome_name_no_space = genome_name_no_space_no_parenthesis_tmp.replace(')', '')

    GenBank_FTP = genome_record_split[col_index['FTP Path']]
    assembly_id = genome_record_split[col_index['Assembly Accession']]
    GenBank_FTP_id = GenBank_FTP.strip().split('/')[-1]

    # prepare cmds
    fna_file =                  '%s_genomic.fna.gz'                     % (GenBank_FTP_id)
    faa_file =                  '%s_protein.faa.gz'                     % (GenBank_FTP_id)
    gbff_file =                 '%s_genomic.gbff.gz'                    % (GenBank_FTP_id)
    pwd_fna_file =              '%s/%s'                                 % (downloaded_genome_folder, fna_file)
    pwd_faa_file =              '%s/%s'                                 % (downloaded_genome_folder, faa_file)
    pwd_gbff_file =             '%s/%s'                                 % (downloaded_genome_folder, gbff_file)
    ftp_fna_file =              '%s/%s'                                 % (GenBank_FTP, fna_file)
    ftp_faa_file =              '%s/%s'                                 % (GenBank_FTP, faa_file)
    ftp_gbff_file =             '%s/%s'                                 % (GenBank_FTP, gbff_file)
    wget_fna_cmd =              'wget %s -P %s -q'                      % (ftp_fna_file, downloaded_genome_folder)
    wget_faa_cmd =              'wget %s -P %s -q'                      % (ftp_faa_file, downloaded_genome_folder)
    wget_gbff_cmd =             'wget %s -P %s -q'                      % (ftp_gbff_file, downloaded_genome_folder)
    gunzip_fna_cmd =            'gunzip %s'                             % (pwd_fna_file)
    gunzip_faa_cmd =            'gunzip %s'                             % (pwd_faa_file)
    gunzip_gbff_cmd =           'gunzip %s'                             % (pwd_gbff_file)
    rename_fna_cmd =            'mv %s/%s_genomic.fna %s/%s.fna'        % (downloaded_genome_folder, GenBank_FTP_id, downloaded_genome_folder, assembly_id)
    rename_faa_cmd =            'mv %s/%s_protein.faa %s/%s.faa'        % (downloaded_genome_folder, GenBank_FTP_id, downloaded_genome_folder, assembly_id)
    rename_gbff_cmd =           'mv %s/%s_genomic.gbff %s/%s.gbff'      % (downloaded_genome_folder, GenBank_FTP_id, downloaded_genome_folder, assembly_id)
    rename_fna_cmd_with_name =  'mv %s/%s_genomic.fna %s/%s_%s.fna'     % (downloaded_genome_folder, GenBank_FTP_id, downloaded_genome_folder, assembly_id, genome_name_no_space)
    rename_faa_cmd_with_name =  'mv %s/%s_protein.faa %s/%s_%s.faa'     % (downloaded_genome_folder, GenBank_FTP_id, downloaded_genome_folder, assembly_id, genome_name_no_space)
    rename_gbff_cmd_with_name = 'mv %s/%s_genomic.gbff %s/%s_%s.gbff'   % (downloaded_genome_folder, GenBank_FTP_id, downloaded_genome_folder, assembly_id, genome_name_no_space)

    # download, decompress and rename
    if get_fna is True:
        os.system(wget_fna_cmd)
        if os.path.isfile(pwd_fna_file) is True:
            os.system(gunzip_fna_cmd)
            if with_name is False:
                os.system(rename_fna_cmd)
            else:
                os.system(rename_fna_cmd_with_name)
        else:
            print('fna file not found in: %s/' % (GenBank_FTP))

    if get_faa is True:
        os.system(wget_faa_cmd)
        if os.path.isfile(pwd_faa_file) is True:
            os.system(gunzip_faa_cmd)
            if with_name is False:
                os.system(rename_faa_cmd)
            else:
                os.system(rename_faa_cmd_with_name)
        else:
            print('faa file not found in: %s/' % (GenBank_FTP))

    if get_gbff is True:
        os.system(wget_gbff_cmd)
        if os.path.isfile(pwd_gbff_file) is True:
            os.system(gunzip_gbff_cmd)
            if with_name is False:
                os.system(rename_gbff_cmd)
            else:
                os.system(rename_gbff_cmd_with_name)
        else:
            print('gbff file not found in: %s/' % (GenBank_FTP))


def download_GenBank_genome(args):

    csv_file =          args['db']
    output_folder =     args['out']
    assembly_id_file =  args['id']
    get_fna =           args['fna']
    get_faa =           args['faa']
    get_gbff =          args['gbff']
    with_name =         args['name']
    num_threads =       args['t']
    force_overwrite =   args['force']

    time_format = '[%Y-%m-%d %H:%M:%S] '

    if (get_fna is False) and (get_faa is False) and (get_gbff is False):
        print(datetime.now().strftime(time_format) + 'Please specify at least one file type to download, program exited')
        exit()

    if os.path.isdir(output_folder) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % output_folder)
        else:
            print(datetime.now().strftime(time_format) + 'Specified genome directory detected, program exited!')
            exit()

    os.mkdir(output_folder)

    # report
    print(datetime.now().strftime(time_format) + 'Downloading genomes with %s cores' % (num_threads))

    provided_version_dict = {}
    assembly_id_set = set()
    assembly_id_set_no_version = set()
    if assembly_id_file is not None:
        for each_id in open(assembly_id_file):
            assembly_id_set.add(each_id.strip())
            assembly_id_set_no_version.add(each_id.strip().split('.')[0])
            provided_version_dict[each_id.strip().split('.')[0]] = each_id.strip()

    # download genome with multiprocessing
    different_version_found_dict = {}
    genomes_in_csv = set()
    genomes_in_csv_no_version = set()
    list_for_multiple_arguments_download_worker = []
    col_index = {}
    for genome_record in open(csv_file):
        genome_record_split = genome_record.strip().split('\t')
        if genome_record.startswith('#Organism/Name'):
            col_index = {key: i for i, key in enumerate(genome_record_split)}
        else:
            assembly_id = genome_record_split[col_index['Assembly Accession']]
            assembly_id_no_version = assembly_id.split('.')[0]
            if assembly_id_file is None:
                list_for_multiple_arguments_download_worker.append([genome_record_split, col_index, output_folder, get_fna, get_faa, get_gbff, with_name])
                genomes_in_csv.add(assembly_id)
                genomes_in_csv_no_version.add(assembly_id.split('.')[0])
            else:
                if assembly_id in assembly_id_set:
                    list_for_multiple_arguments_download_worker.append([genome_record_split, col_index, output_folder, get_fna, get_faa, get_gbff, with_name])
                elif assembly_id_no_version in assembly_id_set_no_version:
                    list_for_multiple_arguments_download_worker.append([genome_record_split, col_index, output_folder, get_fna, get_faa, get_gbff, with_name])
                    different_version_found_dict[assembly_id_no_version] = assembly_id

    # run with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(genome_download_worker, list_for_multiple_arguments_download_worker)
    pool.close()
    pool.join()

    # remove unzip broken files
    os.system('rm %s/*.gz' % output_folder)

    # get stat
    downloaded_file_re = '%s/*_*' % output_folder
    downloaded_file_list = [os.path.basename(file_name) for file_name in glob.glob(downloaded_file_re)]
    downloaded_file_list_no_version = [i.split('.')[0] for i in downloaded_file_list]
    failed_gnm_set = set()
    if assembly_id_file is None:
        for each_gnm in genomes_in_csv:
            each_gnm_no_version = each_gnm.split('.')[0]
            if each_gnm_no_version not in downloaded_file_list_no_version:
                failed_gnm_set.add(each_gnm)
    else:
        for each_provided_gnm in assembly_id_set:
            provided_gnm_no_version = each_provided_gnm.split('.')[0]
            if provided_gnm_no_version not in downloaded_file_list_no_version:
                failed_gnm_set.add(each_provided_gnm)

    if len(failed_gnm_set) > 0:
        print(datetime.now().strftime(time_format) + 'Failed genomes number: %s' % len(failed_gnm_set))
        failed_gnm_txt = '%s_failed.txt'% output_folder
        failed_gnm_txt_handle = open(failed_gnm_txt, 'w')
        for each_failed_gnm in failed_gnm_set:
            failed_gnm_txt_handle.write('%s\n' % each_failed_gnm)
        failed_gnm_txt_handle.close()
        print(datetime.now().strftime(time_format) + 'See details in: %s' % failed_gnm_txt)

    if len(different_version_found_dict) > 0:
        print(datetime.now().strftime(time_format) + 'Different version found for %s genomes' % len(different_version_found_dict))
        diff_version_gnm_txt = '%s_diff_version.txt'% output_folder
        diff_version_gnm_txt_handle = open(diff_version_gnm_txt, 'w')
        diff_version_gnm_txt_handle.write('Provided\tDownloaded\n')
        for each_diff_version_gnm in different_version_found_dict:
            diff_version_gnm_txt_handle.write('%s\t%s\n' % (provided_version_dict[each_diff_version_gnm], different_version_found_dict[each_diff_version_gnm]))
        diff_version_gnm_txt_handle.close()
        print(datetime.now().strftime(time_format) + 'See details in: %s' % diff_version_gnm_txt)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-db',          required=True,                       help='local path to prokaryotes.txt')
    parser.add_argument('-out',         required=True,                       help='output folder')
    parser.add_argument('-id',          required=False, default=None,        help='IDs of genomes to download')
    parser.add_argument('-fna',         required=False, action="store_true", help='download fna file')
    parser.add_argument('-faa',         required=False, action="store_true", help='download faa file')
    parser.add_argument('-gbff',        required=False, action="store_true", help='download gbff file')
    parser.add_argument('-name',        required=False, action="store_true", help='include genome name in the downloaded files')
    parser.add_argument('-force',       required=False, action="store_true", help='force overwrite existing genome directory')
    parser.add_argument('-t',           required=False, default=1, type=int, help='number of threads')
    args = vars(parser.parse_args())
    download_GenBank_genome(args)

'''
BioSAK get_genome_NCBI -csv /Users/songweizhi/DB/NCBI/prokaryotes-2022-05-18.txt -out gnm_dir -fna -t 10
python3 ~/PycharmProjects/BioSAK/BioSAK/get_genome_NCBI.py -csv /Users/songweizhi/DB/NCBI/prokaryotes-2022-05-18.txt -out gnm_dir -fna -t 10

python3 ~/PycharmProjects/BioSAK/BioSAK/get_genome_NCBI.py -db /Users/songweizhi/DB/NCBI/prokaryotes-2022-05-18.txt -id prokaryotes_to_download.txt -out gnm_dir -fna -t 10

'''