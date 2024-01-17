import os
import argparse
import subprocess
import multiprocessing as mp

MetaBiosample_usage = '''
==================== MetaBiosample example commands ====================

BioSAK metaBiosample -i SAMEA5126984 -a 'isolation_source,host,geo_loc_name,lat_lon' -o op_dir
BioSAK metaBiosample -i biosamples.txt -a 'isolation_source,host,geo_loc_name' -o op_dir -t 6

Dependencies: xtract, esearch and esummary
# https://www.ncbi.nlm.nih.gov/books/NBK179288/

=======================================================================
'''

def check_executables(program_list):

    not_detected_programs = []
    for needed_program in program_list:
        if subprocess.call(['which', needed_program], stdout=open(os.devnull, 'wb')) != 0:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not detected, program exited!' % ','.join(not_detected_programs))
        exit()


def MetaBiosample(args):

    file_in         = args['i']
    att_str         = args['a']
    op_dir          = args['o']
    thread_num      = args['t']
    force_overwrite = args['f']
    execute_cmd     = args['exe']

    check_executables(['xtract', 'esearch', 'esummary'])

    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('specified output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s/tmp' % op_dir)

    attributes_to_extract = att_str.split(',')
    cmd_att_part = ''
    for each_attr_to_extract in attributes_to_extract:
        cmd_att_part += '-group Attribute -if Attribute@harmonized_name -equals "%s" -element Attribute ' % each_attr_to_extract

    accession_set = set()
    if os.path.isfile(file_in) is True:
        for each_id in open(file_in):
            accession_set.add(each_id.strip())
    else:
        accession_set.add(file_in)

    esearch_cmd_list = []
    for each_id in accession_set:
        biosample_id = each_id.strip()
        op_meta_txt = '%s/tmp/%s_metadata.txt' % (op_dir, biosample_id)
        xtract_cmd  = 'xtract -pattern DocumentSummary -element Accession -element Date -first Title %s' % cmd_att_part
        esearch_cmd = 'esearch -db biosample -query %s | esummary | %s> %s' % (biosample_id, xtract_cmd, op_meta_txt)
        esearch_cmd_list.append(esearch_cmd)

    cat_cmd = 'cat %s/tmp/*_metadata.txt > %s/combined_metadata.txt' % (op_dir, op_dir)

    cmd_txt = '%s/cmd.txt' % op_dir
    with open(cmd_txt, 'w') as cmd_txt_handle:
        cmd_txt_handle.write('\n'.join(sorted(esearch_cmd_list)) + '\n')
        cmd_txt_handle.write(cat_cmd)

    if execute_cmd is True:
        pool = mp.Pool(processes=thread_num)
        pool.map(os.system, esearch_cmd_list)
        pool.close()
        pool.join()
        os.system(cat_cmd)


if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser(usage=MetaBiosample_usage)
    arg_parser.add_argument('-i',   required=True,                        help='biosample id file, one id per line')
    arg_parser.add_argument('-a',   required=True,                        help='attributes to extract')
    arg_parser.add_argument('-o',   required=True,                        help='output folder')
    arg_parser.add_argument('-t',   required=False, type=int, default=1,  help='number of threads, default: 1')
    arg_parser.add_argument('-f',   required=False, action="store_true",  help='Force overwrite existing results')
    arg_parser.add_argument('-exe', required=False, action="store_true",  help='execute cmds')
    args = vars(arg_parser.parse_args())
    MetaBiosample(args)
