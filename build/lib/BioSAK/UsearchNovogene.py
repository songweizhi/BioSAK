import os
import glob
import argparse
from Bio import SeqIO


UsearchNovogene_usage = '''
======================== UsearchNovogene example commands ========================

BioSAK UsearchNovogene -i CleanData -x fna -r SILVA_138.2.fa -o op_dir -t 12 -f

# SILVA reference sequences on Mac
/Users/songweizhi/DB/SILVA/138.2/SILVA_138.2_SSURef_NR99_tax_silva.fasta

cd /Users/songweizhi/Desktop/SpongeMicrobiomeProject
BioSAK UsearchNovogene -i s01_CleanData -x fna -o demo_op_dir -t 10 -f -r /Users/songweizhi/DB/SILVA/138.2/SILVA_138.2_SSURef_NR99_tax_silva.fasta

# This is a wrapper for the following steps:


==================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def best_hit(file_in, file_out):

    file_out_handle = open(file_out, 'w')
    best_hit_line = ''
    best_hit_query_id = ''
    best_hit_score = 0
    for blast_hit in open(file_in):
        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]
        bit_score = float(blast_hit_split[11])

        if best_hit_query_id == '':
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

        elif (query_id == best_hit_query_id) and (bit_score > best_hit_score):
            best_hit_score = bit_score
            best_hit_line = blast_hit

        elif query_id != best_hit_query_id:
            file_out_handle.write(best_hit_line)
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

    file_out_handle.write(best_hit_line)
    file_out_handle.close()


def UsearchNovogene(args):

    clean_data_dir      = args['i']
    clean_data_ext      = args['x']
    silva_ref_seq       = args['r']
    op_dir              = args['o']
    run_blca            = args['blca']
    blca_ref_seq        = args['ref']
    blca_ref_tax        = args['tax']
    num_threads         = args['t']
    force_create_op_dir = args['f']

    # define file name
    dir_clean_data_renamed  = '%s/s01_CleanData_renamed'                        % op_dir
    dir_DereplicatedData    = '%s/s02_DereplicatedData'                         % op_dir
    s07_OtuTable            = '%s/s07_AllSamples_unoise_otu_table1.txt'         % op_dir
    s09_FinalOtuTable       = '%s/s09_AllSamples_unoise_otu_table_final.txt'    % op_dir
    s10_blca_op_dir         = '%s/s10_AllSamples_BLCA_classifications'          % op_dir

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % dir_clean_data_renamed)
    os.system('mkdir %s' % dir_DereplicatedData)

    clean_data_re = '%s/*.%s' % (clean_data_dir, clean_data_ext)
    clean_data_list = glob.glob(clean_data_re)

    for each_file in clean_data_list:

        fa_name, fa_path, fa_base, fa_ext = sep_path_basename_ext(each_file)

        pwd_fa_renamed = '%s/%s'                % (dir_clean_data_renamed, fa_name)
        pwd_fa_uniques = '%s/%s_uniques.fasta'  % (dir_DereplicatedData, fa_base)

        # add sample id to contig id
        pwd_fa_renamed_handle = open(pwd_fa_renamed, 'w')
        for each_seq in SeqIO.parse(each_file, 'fasta'):
            seq_id = each_seq.id
            pwd_fa_renamed_handle.write('>%s;sample=%s;\n' % (seq_id, fa_base))
            pwd_fa_renamed_handle.write(str(each_seq.seq) + '\n')
        pwd_fa_renamed_handle.close()

        # run dereplicate
        cmd_fastx_uniques = 'usearch -fastx_uniques %s -fastaout %s -sizeout' % (pwd_fa_renamed, pwd_fa_uniques)
        print(cmd_fastx_uniques)
        os.system(cmd_fastx_uniques)

    # combine sequences from all samples
    cmd_2 = 'cat %s/*.fasta > %s/s03_AllSamples.fasta' % (dir_DereplicatedData, op_dir)
    print(cmd_2)
    os.system(cmd_2)

    # dereplicate AllSamples.fasta
    cmd_3  = 'usearch -fastx_uniques %s/s03_AllSamples.fasta -fastaout %s/s04_AllSamples_uniques.fasta -sizein -sizeout -strand both' % (op_dir, op_dir)
    print(cmd_3)
    os.system(cmd_3)

    # Generating unique sequences using UNOISE
    cmd_4  = 'usearch -unoise3 %s/s04_AllSamples_uniques.fasta -zotus %s/s05_AllSamples_denoised.fasta' % (op_dir, op_dir)
    print(cmd_4)
    os.system(cmd_4)

    # Chimera Removal
    cmd_5  = 'usearch -uchime2_ref %s/s05_AllSamples_denoised.fasta -db %s -strand plus -mode high_confidence -notmatched %s/s06_AllSamples_unoise_nc.fasta' % (op_dir, silva_ref_seq, op_dir)
    print(cmd_5)
    os.system(cmd_5)

    # generate OTU table, an OTU table is made by the otutab command
    cmd_6_1 = 'usearch -otutab %s/s03_AllSamples.fasta -db %s/s06_AllSamples_unoise_nc.fasta -id 0.97 -otutabout %s/s07_AllSamples_unoise_otu_table1.txt'                               % (op_dir, op_dir, op_dir)
    cmd_6_2 = 'usearch -otutab %s/s03_AllSamples.fasta -db %s/s06_AllSamples_unoise_nc.fasta -id 0.97 -otutabout %s/s07_AllSamples_unoise_otu_table2.txt -maxaccepts 0 -maxrejects 0'   % (op_dir, op_dir, op_dir)
    print(cmd_6_1)
    os.system(cmd_6_1)
    print(cmd_6_2)
    os.system(cmd_6_2)

    # Mapping of OTUs on Reference Database
    if run_blca is False:

        cmd_7  = 'blastn -query %s/s06_AllSamples_unoise_nc.fasta -outfmt 6 -out %s/s08_AllSamples_unoise_nc.txt -db %s -evalue 1e-20 -num_threads %s' % (op_dir, op_dir, silva_ref_seq, num_threads)
        print(cmd_7)
        os.system(cmd_7)

        # keep best hit
        blast_op          = '%s/s08_AllSamples_unoise_nc.txt'   % op_dir
        blast_op_best_hit = '%s/s08_AllSamples_unoise_nc2.txt'  % op_dir
        best_hit(blast_op, blast_op_best_hit)

        # Merge Table and Taxonomy
        silva_ref_to_tax_dict = dict()
        for each_seq in SeqIO.parse(silva_ref_seq, 'fasta'):
            seq_id = each_seq.id
            seq_tax = ' '.join(each_seq.description.split(' ')[1:])
            silva_ref_to_tax_dict[seq_id] = seq_tax

        otu_best_hit_dict = dict()
        for each in open(blast_op_best_hit):
            each_split = each.strip().split('\t')
            info_list = [each_split[1], each_split[2], each_split[3], each_split[10]]
            otu_best_hit_dict[each_split[0]] = info_list

        s09_FinalOtuTable_handle = open(s09_FinalOtuTable, 'w')
        for each_line in open(s07_OtuTable):
            if each_line.startswith('#'):
                s09_FinalOtuTable_handle.write('%s\tSilvaId\tIdentity\tAlignment_Length\tEvalue\ttaxonomy\n' % each_line.strip())
            else:
                each_line_split = each_line.strip().split('\t')
                otu_id = each_line_split[0]
                best_hit_info_list = otu_best_hit_dict.get(otu_id, ['na', 'na', 'na', 'na'])
                best_hit_info_str = '\t'.join(best_hit_info_list)
                ref_id = best_hit_info_list[0]
                ref_tax = silva_ref_to_tax_dict.get(ref_id, 'na')
                s09_FinalOtuTable_handle.write('%s\t%s\t%s\n' % (each_line.strip(), best_hit_info_str, ref_tax))
        s09_FinalOtuTable_handle.close()

    else:
        print('classify sequences with BLCA')
        from blca import blca
        blca_arg_dict = dict()
        blca_arg_dict['i']   = op_dir
        blca_arg_dict['x']   = '_unoise_nc.fasta'
        blca_arg_dict['o']   = s10_blca_op_dir
        blca_arg_dict['ref'] = blca_ref_seq
        blca_arg_dict['tax'] = blca_ref_tax
        blca_arg_dict['t']   = num_threads
        blca_arg_dict['f']   = True
        blca(blca_arg_dict)


if __name__ == '__main__':

    UsearchNovogene_parser = argparse.ArgumentParser(usage=UsearchNovogene_usage)
    UsearchNovogene_parser.add_argument('-i',       required=True,                        help='path to input sequences')
    UsearchNovogene_parser.add_argument('-x',       required=True,                        help='file extension')
    UsearchNovogene_parser.add_argument('-r',       required=True,                        help='SSU references, e.g., SILVA_138.2_SSURef_NR99_tax_silva.fasta')
    UsearchNovogene_parser.add_argument('-o',       required=True,                        help='output directory')
    UsearchNovogene_parser.add_argument('-blca',    required=False, action="store_true",  help='perform classification with BLCA')
    UsearchNovogene_parser.add_argument('-ref',     required=False,                       help='BLCA reference sequences')
    UsearchNovogene_parser.add_argument('-tax',     required=False,                       help='BLCA reference taxonomy')
    UsearchNovogene_parser.add_argument('-t',       required=False, type=int, default=1,  help='number of threads, default is 1')
    UsearchNovogene_parser.add_argument('-f',       required=False, action="store_true",  help='force overwrite')
    args = vars(UsearchNovogene_parser.parse_args())
    UsearchNovogene(args)
