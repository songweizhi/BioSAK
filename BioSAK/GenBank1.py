import os
import glob
import argparse


GenBank1_usage = '''
====================== GenBank1 example commands ======================

BioSAK GenBank1 -i accession.txt -o op_dir -f -seq 
BioSAK GenBank1 -i accession.txt -o op_dir -f -seq -organism
BioSAK GenBank1 -i accession.txt -o op_dir -f -seq -organism -voucher

=======================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def combine_esearch_op_organism(op_dir, op_txt):

    file_re = '%s/*_organism.txt' % op_dir
    file_list = glob.glob(file_re)

    op_txt_handle = open(op_txt, 'w')
    for each_file in file_list:
        f_name, f_path, f_base, f_ext = sep_path_basename_ext(each_file)
        accession_id = f_base.split('_organism')[0]

        organism_info = ''
        with open(each_file) as f:
            organism_info = f.readline().replace('ORGANISM', '').strip()
        op_txt_handle.write('%s\t%s\n' % (accession_id, organism_info))
    op_txt_handle.close()


def combine_esearch_op_voucher(op_dir, op_txt):

    file_re = '%s/*_voucher.txt' % op_dir
    file_list = glob.glob(file_re)

    op_txt_handle = open(op_txt, 'w')
    for each_file in file_list:
        f_name, f_path, f_base, f_ext = sep_path_basename_ext(each_file)
        accession_id = f_base.split('_voucher')[0]

        voucher_info = ''
        with open(each_file) as f:
            #voucher_info = f.readline().replace('ORGANISM', '').strip()
            voucher_info = f.readline().strip()
            if 'specimen_voucher=' in voucher_info:
                voucher_info = voucher_info.replace('/specimen_voucher="', '')[:-1]
        if voucher_info != '':
            op_txt_handle.write('%s\t%s\n' % (accession_id, voucher_info))
    op_txt_handle.close()


def GenBank1(args):

    accession_txt       = args['i']
    op_dir              = args['o']
    get_sequence        = args['seq']
    get_organism_info   = args['organism']
    get_voucher_info    = args['voucher']
    force_overwrite     = args['f']

    tmp_dir             = '%s/tmp'                       % op_dir
    cmd_txt_fasta       = '%s/cmds_get_fasta.sh'         % op_dir
    cmd_txt_organism    = '%s/cmds_get_organism_info.sh' % op_dir
    cmd_txt_voucher     = '%s/cmds_get_voucher_info.sh'  % op_dir

    if os.path.isdir(op_dir):
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('%s already exists, program exited!' % op_dir)
            exit()
    os.mkdir(op_dir)
    os.mkdir(tmp_dir)

    cmd_list_get_fasta = []
    cmd_list_get_organism = []
    cmd_list_get_voucher = []
    for accession_id in open(accession_txt):
        accession_id     = accession_id.strip().split()[0]
        get_fasta_cmd    = 'esearch -db nucleotide -query %s | efetch -format fasta > %s/%s.fasta'                              % (accession_id, 'tmp', accession_id)
        get_organism_cmd = 'esearch -db nucleotide -query %s | efetch -format gb | grep "ORGANISM" > %s/%s_organism.txt'        % (accession_id, 'tmp', accession_id)
        get_voucher_cmd  = 'esearch -db nucleotide -query %s | efetch -format gb | grep "specimen_voucher" > %s/%s_voucher.txt' % (accession_id, 'tmp', accession_id)
        cmd_list_get_fasta.append(get_fasta_cmd)
        cmd_list_get_organism.append(get_organism_cmd)
        cmd_list_get_voucher.append(get_voucher_cmd)

    if get_sequence is True:
        with open(cmd_txt_fasta, 'w') as f:
            f.write('\n'.join(cmd_list_get_fasta))

    if get_organism_info is True:
        with open(cmd_txt_organism, 'w') as f:
            f.write('\n'.join(cmd_list_get_organism))

    if get_voucher_info is True:
        with open(cmd_txt_voucher, 'w') as f:
            f.write('\n'.join(cmd_list_get_voucher))

    # report
    if get_sequence is True:
        print('Commands for getting sequence info exported to:\t%s' % cmd_txt_fasta)

    if get_organism_info is True:
        print('Commands for getting organism info exported to:\t%s' % cmd_txt_organism)

    if get_voucher_info is True:
        print('Commands for getting voucher info exported to:\t%s'  % cmd_txt_voucher)


if __name__ == '__main__':

    GenBank1_parser = argparse.ArgumentParser(usage=GenBank1_usage)
    GenBank1_parser.add_argument('-i',          required=True,                       help='input txt containing accession id')
    GenBank1_parser.add_argument('-o',          required=True,                       help='output dir')
    GenBank1_parser.add_argument('-seq',        required=False, action="store_true", help='get sequences')
    GenBank1_parser.add_argument('-organism',   required=False, action="store_true", help='get organism info')
    GenBank1_parser.add_argument('-voucher',    required=False, action="store_true", help='get voucher info')
    GenBank1_parser.add_argument('-f',          required=False, action="store_true", help='force overwrite')
    args = vars(GenBank1_parser.parse_args())
    GenBank1(args)
