import os
import argparse
from Bio import AlignIO


SliceMSA_usage = '''
========================= SliceMSA example commands =========================

TreeSAK SliceMSA -i 16S_aln.fasta -s 200-300 -o 16S_aln_200-300.fasta
TreeSAK SliceMSA -i 16S_aln.phylip -fi phylip-relaxed -s sections.txt -o SliceMSA_op -fo phylip-relaxed

# example
200-300     select columns 200-300
-100        select columns 1-300
500-        select columns from 500 to the end

# Example of sections.txt (one section per line):
200-300
-100
500-

# Examples of alignment format (https://biopython.org/wiki/AlignIO): 
fasta, phylip, phylip-relaxed, phylip-sequential, clustal

=============================================================================
'''


def msa2fasta(msa_object, fasta_out):

    with open(fasta_out, 'w') as fasta_out_handle:
        for each_seq in msa_object:
            fasta_out_handle.write('>%s\n' % each_seq.id)
            fasta_out_handle.write('%s\n' % str(each_seq.seq))


def msa2phylip(msa_object, phylip_out):

    max_seq_id_len = 0
    for each_seq in msa_object:
        seq_id_len = len(each_seq.id)
        if seq_id_len > max_seq_id_len:
            max_seq_id_len = seq_id_len

    with open(phylip_out, 'w') as phylip_out_handle:
        phylip_out_handle.write('%s %s\n' % (len(msa_object), msa_object.get_alignment_length()))
        for each_seq in msa_object:
            seq_id = each_seq.id
            seq_id_with_space = '%s%s' % (seq_id, ' ' * (max_seq_id_len + 2 - len(seq_id)))
            phylip_out_handle.write('%s%s\n' % (seq_id_with_space, str(each_seq.seq)))


def SliceMSA(args):

    msa_in_file         = args['i']
    aln_in_format       = args['fi']
    col_to_select_txt   = args['s']
    op_dir              = args['o']
    aln_out_format      = args['fo']
    force_overwriting   = args['force']

    aln_out_ext = 'fasta'
    if aln_out_format == 'phylip-relaxed':
        aln_out_ext = 'phylip'

    if os.path.isfile(msa_in_file) is False:
        print('Input MSA not found, program exited!')
        exit()

    # read in msa
    msa_in = AlignIO.read(msa_in_file, aln_in_format)

    # parse provided sections
    section_to_select_list = []
    if os.path.isfile(col_to_select_txt) is False:
        col_to_select_txt_split = col_to_select_txt.strip().split('-')
        if col_to_select_txt == '-':
            section_to_select_list.append(['1', str(msa_in.get_alignment_length())])
        elif col_to_select_txt.startswith('-'):
            section_to_select_list.append(['1', col_to_select_txt_split[1]])
        elif col_to_select_txt.endswith('-'):
            section_to_select_list.append([col_to_select_txt_split[0], str(msa_in.get_alignment_length())])
        else:
            section_to_select_list.append(col_to_select_txt_split)
    else:
        for each_section in open(col_to_select_txt):
            each_section = each_section.strip()
            each_section_split = each_section.strip().split('-')
            if each_section == '-':
                section_to_select_list.append(['1', str(msa_in.get_alignment_length())])
            elif each_section.startswith('-'):
                section_to_select_list.append(['1', each_section_split[1]])
            elif each_section.endswith('-'):
                section_to_select_list.append([each_section_split[0], str(msa_in.get_alignment_length())])
            else:
                section_to_select_list.append(each_section_split)

    # check output folder
    if len(section_to_select_list) > 1:
        if os.path.isdir(op_dir) is True:
            if force_overwriting is True:
                os.system('rm -r %s' % op_dir)
            else:
                print('Output folder already exist, program exited!')
                exit()
        os.system('mkdir %s' % op_dir)

    # write out sections
    if len(section_to_select_list) == 1:
        current_section = msa_in[:, (int(section_to_select_list[0][0]) - 1):(int(section_to_select_list[0][1]))]
        if aln_out_ext == 'fasta':
            msa2fasta(current_section, op_dir)
        if aln_out_ext == 'phylip':
            msa2phylip(current_section, op_dir)
    else:
        for each_section in section_to_select_list:

            pwd_op_file     = '%s/%s.%s' % (op_dir, '-'.join(each_section), aln_out_ext)
            current_section = msa_in[:, (int(each_section[0])-1):(int(each_section[1]))]

            # write out
            if aln_out_ext == 'fasta':
                msa2fasta(current_section, pwd_op_file)
            if aln_out_ext == 'phylip':
                msa2phylip(current_section, pwd_op_file)

    print('MSA subset(s) exported to %s, Done!' % op_dir)


if __name__ == '__main__':

    # arguments for rename_seq_parser
    SliceMSA_parser = argparse.ArgumentParser()
    SliceMSA_parser.add_argument('-i',      required=True,                         help='input MSA in fasta format')
    SliceMSA_parser.add_argument('-fi',     required=False, default='fasta',       help='format (NOT file extension) of input MSA, default: fasta')
    SliceMSA_parser.add_argument('-s',      required=True,                         help='columns to export, e.g. 200-300, -100, 50-')
    SliceMSA_parser.add_argument('-o',      required=True,                         help='output file or folder')
    SliceMSA_parser.add_argument('-fo',     required=False, default='fasta',       help='format of output MSA, select from fasta and phylip-relaxed, default: fasta')
    SliceMSA_parser.add_argument('-force',  required=False, action="store_true",   help='force overwrite existing output folder')
    args = vars(SliceMSA_parser.parse_args())
    SliceMSA(args)

