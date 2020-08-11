import os
import shutil
import argparse
from requests import get
import multiprocessing as mp
from datetime import datetime
from bs4 import BeautifulSoup


''' 
cd /Users/songweizhi/Desktop
python3 ~/PycharmProjects/BioSAK/BioSAK/get_ko_gene_seqs.py -ko interested_ko.txt -o test_test -force -t 4
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


def download_ko_genes(argument_list):

    ko              = argument_list[0]
    fasta_file_out  = argument_list[1]


    ko_url = 'https://www.genome.jp/dbget-bin/www_bget?ko:%s' % ko
    web_address_K18831_soup = BeautifulSoup(get(ko_url).text, features="lxml")

    # get total number of genes
    total_gene_num = 0
    for each_row in web_address_K18831_soup.find_all('tr'):
        try:
            row_name = each_row.th.text
            if row_name == 'Genes':
                gene_list = each_row.find_all('a')
                for each_gene in gene_list:
                    if '/dbget-bin/www_bget?' in each_gene.get('href'):
                        total_gene_num += 1
        except:
            pass

    time_format = '[%Y-%m-%d %H:%M:%S]'
    print('%s Downloading %s genes for %s.' % ((datetime.now().strftime(time_format)), total_gene_num, ko))

    # write out sequences
    file_out_handle = open(fasta_file_out, 'w')
    for each_row in web_address_K18831_soup.find_all('tr'):

        try:
            row_name = each_row.th.text
            if row_name == 'Genes':
                gene_list = each_row.find_all('a')

                for each_gene in gene_list:

                    if '/dbget-bin/www_bget?' in each_gene.get('href'):

                        gene_name = each_gene.text
                        gene_link = 'https://www.genome.jp%s' % each_gene.get('href')
                        gene_soup = BeautifulSoup(get(gene_link).text, features="lxml")

                        for gene_soup_row in gene_soup.find_all('tr'):

                            try:
                                gene_soup_row_name = gene_soup_row.th.text

                                if gene_soup_row_name == 'AA seq':

                                    # remove button text
                                    for button in gene_soup_row.find_all('button'):
                                        button.clear()

                                    to_write_out = '>%s %s\n' % (gene_name, gene_soup_row.td.text)

                                    file_out_handle.write(to_write_out)

                            except:
                                pass
        except:
            pass

    file_out_handle.close()


def main_code(args):

    interested_ko_file  = args['ko']
    output_folder       = args['o']
    force_overwrite     = args['force']
    num_threads         = args['t']

    # create output folder
    if (os.path.isdir(output_folder) is True) and (force_overwrite is False):
        print('Output folder already exists, program exited!')
        exit()
    else:
        force_create_folder(output_folder)

    # store interested ko in list
    interested_ko_list = []
    for interested_ko in open(interested_ko_file):
        interested_ko_list.append(interested_ko.strip())

    list_for_multiple_arguments_download_ko_genes = []
    for interested_ko in interested_ko_list:
        interested_ko_fasta_output = '%s/%s.fasta' % (output_folder, interested_ko)
        list_for_multiple_arguments_download_ko_genes.append([interested_ko, interested_ko_fasta_output])

    pool = mp.Pool(processes=num_threads)
    pool.map(download_ko_genes, list_for_multiple_arguments_download_ko_genes)
    pool.close()
    pool.join()


if __name__ == '__main__':

    download_ko_genes_parser = argparse.ArgumentParser()

    # arguments for COG_parser
    download_ko_genes_parser.add_argument('-ko',              required=True,                              help='interested KOs at level D, one id per line')
    download_ko_genes_parser.add_argument('-o',               required=False,                             help='output folder')
    download_ko_genes_parser.add_argument('-force',           required=False, action='store_true',        help='force overwrite')
    download_ko_genes_parser.add_argument('-t',               required=False, type=int, default=1,        help='number of threads')

    args = vars(download_ko_genes_parser.parse_args())
    main_code(args)
