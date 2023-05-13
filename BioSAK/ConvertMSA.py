import os
import glob
import argparse
from Bio import SeqIO
from Bio import AlignIO


ConvertMSA_usage = '''
================================= ConvertMSA example commands =================================

# phylip to fasta
TreeSAK ConvertMSA -i concatenated.phy -fi phylip-relaxed -o concatenated.fasta -fo fasta
TreeSAK ConvertMSA -i phy_files -fi phylip-relaxed -xi phy -o MSA_in_fasta -fo fasta -xo fa 

# examples of alignment format (https://biopython.org/wiki/AlignIO): 
fasta, phylip, phylip-relaxed, phylip-sequential, clustal

===============================================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def ConvertMSA(args):

    aln_in            = args['i']
    aln_in_ext        = args['xi']
    aln_in_format     = args['fi']
    aln_out           = args['o']
    aln_out_ext       = args['xo']
    aln_out_format    = args['fo']
    one_line          = args['oneline']
    no_gap            = args['nogap']
    force_overwriting = args['f']

    if ((one_line is True) and (aln_out_format != 'fasta')) or ((no_gap is True) and (aln_out_format != 'fasta')):
        print('Please provide "-oneline" and/or "-nogap" only if "-fo" is fasta')
        exit()

    if os.path.isfile(aln_in) is True:
        if (one_line is False) and (no_gap is False):
            AlignIO.convert(aln_in, aln_in_format, aln_out, aln_out_format)
        else:
            aln_out_tmp = aln_out + '.tmp'
            AlignIO.convert(aln_in, aln_in_format, aln_out_tmp, aln_out_format)
            pwd_aln_out_handle = open(aln_out, 'w')
            for each_seq in SeqIO.parse(aln_out_tmp, 'fasta'):
                seq_id = each_seq.id
                seq_sequence = str(each_seq.seq)
                if no_gap is False:
                    pwd_aln_out_handle.write('>%s\n' % seq_id)
                    pwd_aln_out_handle.write('%s\n' % seq_sequence)
                else:
                    pwd_aln_out_handle.write('>%s\n' % seq_id)
                    pwd_aln_out_handle.write('%s\n' % seq_sequence.replace('-', ''))
            pwd_aln_out_handle.close()
            os.system('rm %s' % aln_out_tmp)

        print('Done!')

    elif os.path.isdir(aln_in) is True:
        aln_in_re = '%s/*.%s' % (aln_in, aln_in_ext)
        aln_in_list = [os.path.basename(file_name) for file_name in glob.glob(aln_in_re)]

        # check input
        if len(aln_in_list) == 0:
            print('Input file not detected, program exited!')
            exit()

        # check output folder
        if os.path.isdir(aln_out) is True:
            if force_overwriting is True:
                os.system('rm -r %s' % aln_out)
            else:
                print('Output folder already exist, program exited!')
                exit()
        os.system('mkdir %s' % aln_out)

        # convert
        for each_aln_in in aln_in_list:

            aln_in_path, aln_in_basename, aln_in_ext = sep_path_basename_ext(each_aln_in)
            pwd_aln_in      = '%s/%s'        % (aln_in, each_aln_in)
            pwd_aln_out     = '%s/%s.%s'     % (aln_out, aln_in_basename, aln_out_ext)
            pwd_aln_out_tmp = '%s/%s_tmp.%s' % (aln_out, aln_in_basename, aln_out_ext)

            if (one_line is False) and (no_gap is False):
                AlignIO.convert(pwd_aln_in, aln_in_format, pwd_aln_out, aln_out_format)
            else:
                AlignIO.convert(pwd_aln_in, aln_in_format, pwd_aln_out_tmp, aln_out_format)
                pwd_aln_out_handle = open(pwd_aln_out, 'w')
                for each_seq in SeqIO.parse(pwd_aln_out_tmp, 'fasta'):
                    seq_id = each_seq.id
                    seq_sequence = str(each_seq.seq)
                    if no_gap is False:
                        pwd_aln_out_handle.write('>%s\n' % seq_id)
                        pwd_aln_out_handle.write('%s\n' % seq_sequence)
                    else:
                        sequence_no_gap = seq_sequence.replace('-', '')
                        if len(sequence_no_gap) > 0:
                            pwd_aln_out_handle.write('>%s\n' % seq_id)
                            pwd_aln_out_handle.write('%s\n' % sequence_no_gap)
                pwd_aln_out_handle.close()
                os.system('rm %s' % pwd_aln_out_tmp)
        print('Done!')
    else:
        print('Input file not found, program exited!')


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',       required=True,                       help='input alignment')
    parser.add_argument('-xi',      required=False, default='aln',       help='input alignment extension')
    parser.add_argument('-fi',      required=True,                       help='input alignment format, e.g., fasta, phylip')
    parser.add_argument('-o',       required=True,                       help='output alignment')
    parser.add_argument('-xo',      required=False, default='aln',       help='output alignment extension')
    parser.add_argument('-fo',      required=True,                       help='output alignment format, e.g., fasta, phylip')
    parser.add_argument('-oneline', required=False, action="store_true", help='put sequence in single line, available if -fo is fasta')
    parser.add_argument('-nogap',   required=False, action="store_true", help='remove gaps from alignment, available if -fo is fasta')
    parser.add_argument('-f',       required=False, action="store_true", help='force overwrite existing output folder')
    args = vars(parser.parse_args())
    ConvertMSA(args)
