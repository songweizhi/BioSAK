import os
import glob
import argparse
from Bio import SeqIO


GenBank_usage = '''
====================== GenBank example commands ======================

BioSAK GenBank -f -i accession.txt -o op_dir -tax ncbi_taxonomy.txt

-tax /Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt

======================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def gbk2faa(gbk_in, faa_out):

    _, _, faa_base, _ = sep_path_basename_ext(faa_out)

    faa_out_handle = open(faa_out, 'w')
    for seq_record in SeqIO.parse(gbk_in, 'genbank'):
        for feature in seq_record.features:
            if feature.type == 'CDS':
                feature_translation = ''
                if 'translation' in feature.qualifiers:
                    feature_translation = feature.qualifiers['translation'][0]
                if feature_translation != '':
                    faa_out_handle.write('>%s__%s\n' % (faa_base, feature.qualifiers['protein_id'][0]))
                    faa_out_handle.write('%s\n' % feature_translation)
    faa_out_handle.close()


def get_origin_seq_from_gbk(gbk_file):

    after_origin_line = False
    before_slash_line = True
    concatenated_seq = ''
    for each_line in open(gbk_file):
        each_line = each_line.strip()
        if each_line == 'ORIGIN':
            after_origin_line = True
        if each_line == '//':
            before_slash_line = False
        if (after_origin_line is True) and (before_slash_line is True):
            if each_line != 'ORIGIN':
                each_line_split = each_line.split(' ')
                seq_str = ''.join(each_line_split[1:])
                seq_str = seq_str.upper()
                concatenated_seq += seq_str

    return concatenated_seq


def combine_esearch_op_organism(op_dir, tax_lineage_dict, op_txt):

    file_re = '%s/*_organism.txt' % op_dir
    file_list = glob.glob(file_re)

    op_txt_handle = open(op_txt, 'w')
    for each_file in file_list:
        f_name, f_path, f_base, f_ext = sep_path_basename_ext(each_file)
        accession_id = f_base.split('_organism')[0]

        organism_info = ''
        full_lineage_str = ''
        with open(each_file) as f:
            organism_info = f.readline().replace('ORGANISM', '').strip()
            if organism_info != '':
                organism_g = organism_info.split()[0]
                g_full_lineage = tax_lineage_dict.get(organism_g, ('g__' + organism_g))
                full_lineage_str = '%s;s__%s' % (g_full_lineage, organism_info)

        if len(tax_lineage_dict) == 0:
            if organism_info != '':
                op_txt_handle.write('%s\t%s\n' % (accession_id, organism_info))
        else:
            if full_lineage_str != '':
                op_txt_handle.write('%s\t%s\n' % (accession_id, full_lineage_str))
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


def GenBank(args):

    accession_txt       = args['i']
    op_dir              = args['o']
    tax_file            = args['tax']
    get_aa              = args['aa']
    skip_esearch        = args['skip_esearch']
    force_overwrite     = args['f']

    tmp_dir                     = '%s/tmp'                      % op_dir
    cmd_txt_gbk                 = '%s/cmds_get_gbk.sh'          % op_dir
    combined_fa                 = '%s/combined.fasta'           % op_dir
    combined_faa                = '%s/combined.faa'             % op_dir
    combined_organism_info_txt  = '%s/combined_organism.txt'    % op_dir
    combined_voucher_info_txt   = '%s/combined_voucher.txt'     % op_dir

    if skip_esearch is True:
        if force_overwrite is True:
            print('"-skip_esearch" is provided, please remove "-f"' % op_dir)
            exit()
    else:
        if os.path.isdir(op_dir):
            if force_overwrite is True:
                os.system('rm -r %s' % op_dir)
            else:
                print('%s already exists, program exited!' % op_dir)
                exit()
        os.mkdir(op_dir)
        os.mkdir(tmp_dir)

    if skip_esearch is False:
        cmd_list_get_gbk = []
        for accession_id in open(accession_txt):

            accession_id = accession_id.strip().split()[0]
            gbk_file     = '%s/%s.gbk'          % (tmp_dir, accession_id)
            fa_file      = '%s/%s.fasta'        % (tmp_dir, accession_id)
            faa_file      = '%s/%s.faa'         % (tmp_dir, accession_id)
            voucher_txt  = '%s/%s_voucher.txt'  % (tmp_dir, accession_id)
            organism_txt = '%s/%s_organism.txt' % (tmp_dir, accession_id)

            get_gbk_cmd  = 'esearch -db nucleotide -query %s | efetch -format gb > %s' % (accession_id, gbk_file)
            os.system(get_gbk_cmd)
            cmd_list_get_gbk.append(get_gbk_cmd)

            # get organism and voucher info
            os.system('cat %s | grep "ORGANISM" > %s'         % (gbk_file, organism_txt))
            os.system('cat %s | grep "specimen_voucher" > %s' % (gbk_file, voucher_txt))

            # get sequence file
            if os.path.isfile(gbk_file):
                gbk_record = SeqIO.read(gbk_file, "genbank")
                seq_id   = gbk_record.id
                seq_desc = gbk_record.description
                seq_str  = get_origin_seq_from_gbk(gbk_file)
                with open(fa_file, 'w') as f1:
                    f1.write('>%s %s\n' % (seq_id, seq_desc))
                    f1.write('%s\n' % seq_str)

                if get_aa is True:
                    gbk2faa(gbk_file, faa_file)

        # write out command
        with open(cmd_txt_gbk, 'a') as f2:
            f2.write('\n'.join(cmd_list_get_gbk))

    ###################################### combine the results for all accessions ######################################

    tax_lineage_dict = dict()
    if tax_file is not None:
        for each_line in open(tax_file):
            each_line_split = each_line.strip().split('\t')
            tax_lineage_dict[each_line_split[0]] = each_line_split[1]

    os.system('cat %s/*.fasta > %s' % (tmp_dir, combined_fa))
    os.system('cat %s/*.faa > %s'   % (tmp_dir, combined_faa))
    combine_esearch_op_organism(tmp_dir, tax_lineage_dict, combined_organism_info_txt)
    combine_esearch_op_voucher(tmp_dir, combined_voucher_info_txt)


if __name__ == '__main__':

    GenBank_parser = argparse.ArgumentParser(usage=GenBank_usage)
    GenBank_parser.add_argument('-i',               required=True,                       help='input txt containing accession id')
    GenBank_parser.add_argument('-o',               required=True,                       help='output dir')
    GenBank_parser.add_argument('-tax',             required=False, default=None,        help='NCBI taxonomy file, specify to get the full taxon lineage')
    GenBank_parser.add_argument('-aa',              required=False, action="store_true", help='export amino acid sequence')
    GenBank_parser.add_argument('-f',               required=False, action="store_true", help='force overwrite')
    GenBank_parser.add_argument('-skip_esearch',    required=False, action="store_true", help='skip esearch')
    args = vars(GenBank_parser.parse_args())
    GenBank(args)
