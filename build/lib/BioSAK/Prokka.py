import os
import glob
import argparse
import multiprocessing as mp


prokka_usage = '''
======================== prokka example commands ========================

BioSAK prokka -i gnm_folder -x fa -t 12 -d genome_domain.txt
BioSAK prokka -i mag_folder -x fa -t 12 -m genome_type.txt

# gnm_domain.txt (no file extension, tab separated)
# genomes not included in gnm_domain.txt will be treated as "Bacteria"
gnm_1   Archaea
gnm_2   Bacteria
   
# genome_type.txt (no file extension, tab separated)
# genomes not included in gnm_type.txt will be treated as "nonMAG"
gnm_1   MAG
gnm_2   nonMAG

=========================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)
    return f_path, f_base, f_ext


def prokka(args):

    gnm_dir         = args['i']
    gnm_ext         = args['x']
    gnm_domain_txt  = args['d']
    meta_mode_txt   = args['m']
    execute_cmd     = args['e']
    num_threads     = args['t']

    prokka_cmds_txt = '%s_prokka_cmds.txt' % gnm_dir

    # read in gnm domain info
    gnm_domain_dict = dict()
    if os.path.isfile(gnm_domain_txt) is True:
        for each_line in open(gnm_domain_txt):
            each_line_split = each_line.strip().split('\t')
            gnm_id = each_line_split[0]
            gnm_domain = each_line_split[1]
            if gnm_domain not in ['Archaea', 'Bacteria']:
                print('please specify either Archaea or Bacteria')
                exit()
            gnm_domain_dict[gnm_id] = gnm_domain

    # read in gnm type info (MAG or not)
    gnm_type_dict = dict()
    if os.path.isfile(meta_mode_txt) is True:
        for each_line in open(meta_mode_txt):
            each_line_split = each_line.strip().split('\t')
            gnm_id = each_line_split[0]
            gnm_type = each_line_split[1]
            if gnm_type not in ['MAG', 'nonMAG']:
                print('please specify either MAG or nonMAG')
                exit()
            gnm_type_dict[gnm_id] = gnm_type

    gnm_file_re   = '%s/*.%s' % (gnm_dir, gnm_ext)
    gnm_file_list = glob.glob(gnm_file_re)

    prokka_cmd_list = []
    for each_gnm in gnm_file_list:
        _, gnm_id, _ = sep_path_basename_ext(each_gnm)
        gnm_domain = gnm_domain_dict.get(gnm_id, 'Bacteria')
        gnm_type = gnm_type_dict.get(gnm_id, 'nonMAG')

        prokka_cmd = 'prokka --force --compliant --cpus 1 --kingdom %s --prefix %s --locustag %s --strain %s --outdir %s_prokka_wd %s' % (gnm_domain, gnm_id, gnm_id, gnm_id, gnm_id, each_gnm)
        if gnm_type == 'MAG':
            prokka_cmd = 'prokka --force --compliant --metagenome --cpus 1 --kingdom %s --prefix %s --locustag %s --strain %s --outdir %s_prokka_wd %s' % (gnm_domain, gnm_id, gnm_id, gnm_id, gnm_id, each_gnm)
        prokka_cmd_list.append(prokka_cmd)

    # write out command
    prokka_cmds_txt_handle = open(prokka_cmds_txt, 'w')
    for each_cmd in prokka_cmd_list:
        prokka_cmds_txt_handle.write(each_cmd + '\n')
    prokka_cmds_txt_handle.close()

    # run prokka with mp
    if execute_cmd is True:
        print('Annotating %s genomes with %s cores' % (len(prokka_cmd_list), num_threads))
        pool = mp.Pool(processes=num_threads)
        pool.map(os.system, prokka_cmd_list)
        pool.close()
        pool.join()
        print('Done!')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage=prokka_usage)
    parser.add_argument('-i',          required=True,                          help='genome folder')
    parser.add_argument('-x',          required=True, default='fna',           help='file extension, deafult: fna')
    parser.add_argument('-d',          required=True,                          help='genome domain, Bacteria or Archaea')
    parser.add_argument('-m',          required=False, action="store_true",    help='annotate MAG')
    parser.add_argument('-e',          required=False, action="store_true",    help='execute commands')
    parser.add_argument('-t',          required=False, type=int, default=1,    help='number of threads')
    args = vars(parser.parse_args())
    prokka(args)
