import os
import argparse


mmseqs_usage = '''
============================== mmseqs example command ==============================

BioSAK mmseqs -i contigs.fa -o op_dir -db /scratch/PI/ocessongwz/DB/mmseqs/GTDB -f

====================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def mmseqs(args):

    seq_file        = args['i']
    op_dir          = args['o']
    db_dir          = args['db']
    num_threads     = args['t']
    force_overwrite = args['f']

    # create op dir
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output directory detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    _, _, seq_file_base, _ = sep_path_basename_ext(seq_file)

    cmd_1 = 'mmseqs createdb %s %s/%s.DB'                                                               % (seq_file, op_dir, seq_file_base)
    cmd_2 = 'mmseqs taxonomy %s/%s.DB %s %s/taxonomyResult %s/tmp --threads %s'                         % (op_dir, seq_file_base, db_dir, op_dir, op_dir, num_threads)
    cmd_3 = 'mmseqs createtsv %s/%s.DB %s/taxonomyResult %s/taxonomyResult.tsv'                         % (op_dir, seq_file_base, op_dir, op_dir)
    cmd_4 = 'mmseqs taxonomyreport %s %s/taxonomyResult %s/taxonomyResult_report'                       % (db_dir, op_dir, op_dir)
    cmd_5 = 'mmseqs taxonomyreport %s %s/taxonomyResult %s/taxonomyResult_report.html --report-mode 1'  % (db_dir, op_dir, op_dir)

    os.system(cmd_1)
    os.system(cmd_2)
    os.system(cmd_3)
    os.system(cmd_4)
    os.system(cmd_5)

    print('Classification results exported to:')
    print('Kraken format\t%s/taxonomyResult_report'     % op_dir)
    print('Krona format\t%s/taxonomyResult_report.html' % op_dir)
    print('For the visualization of the Kraken report using Pavian: https://github.com/soedinglab/mmseqs2/wiki#taxonomy-assignment')


if __name__ == '__main__':

    mmseqs_parser = argparse.ArgumentParser(usage=mmseqs_usage)
    mmseqs_parser.add_argument('-i',    required=True,                          help='sequence file')
    mmseqs_parser.add_argument('-o',    required=True,                          help='output directory')
    mmseqs_parser.add_argument('-db',   required=True,                          help='db dir')
    mmseqs_parser.add_argument('-t',    required=False, type=int, default=1,    help='number of threads, default is 1')
    mmseqs_parser.add_argument('-f',    required=False, action="store_true",    help='force overwrite')
    args = vars(mmseqs_parser.parse_args())
    mmseqs(args)
