import os
import glob
import argparse
from Bio import SeqIO


get_MAG_reads_long_usage = '''
================================== get_MAG_reads_long example commands ==================================

BioSAK get_MAG_reads_long -r LongReads.fastq -s assembly.sam -m MAG_1.fa -o MAG_1_reads -l 3000
BioSAK get_MAG_reads_long -r LongReads.fastq -s assembly.sam -m MAG_dir -x fa -o extracted_reads -l 0
BioSAK get_MAG_reads_long -r LongReads.fastq -s assembly.sam -m MAG_dir -x fa -o extracted_reads -l 5000

# Note
The input sam file is obtained by mapping quality-filtered reads to all contigs derived from these reads.
minimap2 -d assembly.mmi assembly.fasta
minimap2 --secondary=no -t 12 -a -o assembly.sam assembly.mmi long_reads.fastq

=========================================================================================================
'''


def get_MAG_reads_long(args):

    reads_file          = args['r']
    sam_file            = args['s']
    mag_folder          = args['m']
    mag_ext             = args['x']
    min_read_len        = args['l']
    extracted_read_dir  = args['o']

    # create folder
    if os.path.isdir(extracted_read_dir) is True:
        print('Output folder detected, program exited!')
        print(extracted_read_dir)
        exit()

    mag_list = []
    if os.path.isfile(mag_folder) is True:
        mag_list = [mag_folder]
    elif os.path.isdir(mag_folder) is True:
        mag_re = '%s/*.%s' % (mag_folder, mag_ext)
        mag_list = glob.glob(mag_re)
        if len(mag_list) == 0:
            print('No MAG detected, program exited!')
            exit()
    else:
        print('No MAG file/dir detected, program exited!')
        exit()

    print('getting ctg_to_gnm dict')
    ctg_to_gnm_dict = {}
    ctgs_in_multi_mags_set = set()
    for each_gnm in mag_list:
        gnm_id = '.'.join(os.path.basename(each_gnm).strip().split('.')[:-1])
        for each_seq in SeqIO.parse(each_gnm, 'fasta'):
            ctg_id = each_seq.id
            if ctg_id not in ctg_to_gnm_dict:
                ctg_to_gnm_dict[ctg_id] = {gnm_id}
            else:
                ctg_to_gnm_dict[ctg_id].add(gnm_id)
                ctgs_in_multi_mags_set.add(ctg_id)

    # check if there is contig in multiple MAGs
    if len(ctgs_in_multi_mags_set) > 0:
        print('The following contigs found in multiple MAGs, proceed anyway!')
        print(','.join(sorted([i for i in ctgs_in_multi_mags_set])))

    print('getting read_to_gnm_dict')
    read_to_gnm_dict = {}
    for each_line in open(sam_file):
        if not each_line.startswith('@'):
            each_line_split = each_line.strip().split('\t')
            read_id = each_line_split[0]
            ref_id = each_line_split[2]
            ref_gnm_set = ctg_to_gnm_dict.get(ref_id, set())
            for each_ref_gnm in ref_gnm_set:
                if read_id not in read_to_gnm_dict:
                    read_to_gnm_dict[read_id] = {each_ref_gnm}
                else:
                    read_to_gnm_dict[read_id].add(each_ref_gnm)

    os.system('mkdir %s' % extracted_read_dir)

    print('Extracting reads')
    for each_read in SeqIO.parse(reads_file, 'fastq'):
        read_id = each_read.id
        read_len = len(each_read.seq)
        read_gnms = read_to_gnm_dict.get(read_id, {})
        if read_len >= min_read_len:
            if read_gnms != {}:
                for each_gnm in read_gnms:
                    pwd_read_file = '%s/%s_reads_min%sbp.fastq' % (extracted_read_dir, each_gnm, min_read_len)
                    pwd_read_file_handle = open(pwd_read_file, 'a')
                    SeqIO.write(each_read, pwd_read_file_handle, 'fastq')
                    pwd_read_file_handle.close()

    # print flye cmds
    print('Suggested flye commands:')
    for each_gnm in sorted(mag_list):
        gnm_id = '.'.join(os.path.basename(each_gnm).strip().split('.')[:-1])
        flye_cmd = 'flye --meta --threads 12 --out-dir flye_reassemble_wd_%s --nano-raw %s_reads_min%sbp.fastq' % (gnm_id, gnm_id, min_read_len)
        print(flye_cmd)

    print('Done!')


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-r',    required=True,                             help='reads file')
    parser.add_argument('-s',    required=True,                             help='sam file')
    parser.add_argument('-m',    required=True,                             help='MAG file/folder')
    parser.add_argument('-x',    required=False, default='fasta',           help='MAG extension, default: fasta')
    parser.add_argument('-l',    required=False, type=int, default=5000,    help='minimal read length (in bp), default: 5000')
    parser.add_argument('-o',    required=True,                             help='Output folder')
    args = vars(parser.parse_args())
    get_MAG_reads_long(args)
