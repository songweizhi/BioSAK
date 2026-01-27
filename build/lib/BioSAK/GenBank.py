import os
import glob
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


GenBank_usage = '''
====================== GenBank example commands ======================

BioSAK GenBank -f -i accession.txt -o op_dir -tax ncbi_taxonomy.txt

-tax /Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt

Note: the "-tax" file can be prepared with BioSAK's taxdump module

======================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def gbk_to_ffn_faa(gbk_in, ffn_out, faa_out):

    _, _, gbk_base, _ = sep_path_basename_ext(gbk_in)
    faa_out_handle = open(faa_out, 'w')
    ffn_out_handle = open(ffn_out, 'w')
    for seq_record in SeqIO.parse(gbk_in, 'genbank'):
        record_sequence = str(seq_record.seq)
        for feature in seq_record.features:
            if feature.type == 'CDS':
                feature_pos = '%s_%s' % (feature.location.start, feature.location.end)
                feature_seq_nt = record_sequence[feature.location.start:feature.location.end]
                if feature.location.strand in ['-1', -1]:
                    feature_pos = '%s_%s' % (feature.location.end, feature.location.start)
                    feature_seq_nt = str(Seq(feature_seq_nt).reverse_complement())
                feature_seq_aa = ''
                if 'translation' in feature.qualifiers:
                    feature_seq_aa = feature.qualifiers['translation'][0]

                protein_id = ''
                if 'protein_id' in feature.qualifiers:
                    protein_id = feature.qualifiers['protein_id'][0]

                cds_id_to_used = feature_pos.replace('<', '').replace('>', '')
                if protein_id != '':
                    cds_id_to_used = protein_id

                feature_product = ''
                if 'protein_id' in feature.qualifiers:
                    feature_product = feature.qualifiers['product'][0]

                gene_name = ''
                if 'gene' in feature.qualifiers:
                    gene_name = feature.qualifiers['gene'][0]

                feature_note = ''
                if 'note' in feature.qualifiers:
                    feature_note = feature.qualifiers['note'][0]

                desc_to_use = ''
                if feature_product != '':
                    desc_to_use = feature_product
                elif feature_note != '':
                    desc_to_use = feature_note
                elif gene_name != '':
                    desc_to_use = gene_name

                # write out nt sequence
                ffn_out_handle.write('>%s__%s %s\n' % (gbk_base, cds_id_to_used, desc_to_use))
                ffn_out_handle.write('%s\n' % feature_seq_nt)

                # write out aa sequence
                if feature_seq_aa != '':
                    faa_out_handle.write('>%s__%s %s\n' % (gbk_base, cds_id_to_used, desc_to_use))
                    faa_out_handle.write('%s\n' % feature_seq_aa)
    faa_out_handle.close()
    ffn_out_handle.close()


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
    include_itol_header = args['itol']

    tmp_dir                     = '%s/tmp'                      % op_dir
    cmd_txt_gbk                 = '%s/cmds_get_gbk.sh'          % op_dir
    combined_fa                 = '%s/combined.fasta'           % op_dir
    combined_ffn                = '%s/combined.ffn'             % op_dir
    combined_faa                = '%s/combined.faa'             % op_dir
    combined_organism_info_txt  = '%s/combined_organism.txt'    % op_dir
    combined_voucher_info_txt   = '%s/combined_voucher.txt'     % op_dir
    combined_fa_desc            = '%s/combined_fasta_desc.txt'  % op_dir
    combined_ffn_desc           = '%s/combined_ffn_desc.txt'    % op_dir
    combined_faa_desc           = '%s/combined_faa_desc.txt'    % op_dir

    if skip_esearch is True:
        if force_overwrite is True:
            print('"-skip_esearch" is provided, please remove "-f"')
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

        accession_id_set = set()
        for each_id in open(accession_txt):
            each_id = each_id.strip().split()[0]
            accession_id_set.add(each_id)

        for accession_id in sorted(list(accession_id_set)):
            gbk_file     = '%s/%s.gbk'                      % (tmp_dir, accession_id)
            fa_file      = '%s/%s.fasta'                    % (tmp_dir, accession_id)
            ffn_file      = '%s/%s.ffn'                     % (tmp_dir, accession_id)
            faa_file      = '%s/%s.faa'                     % (tmp_dir, accession_id)
            voucher_txt  = '%s/%s_voucher.txt'              % (tmp_dir, accession_id)
            organism_txt = '%s/%s_organism.txt'             % (tmp_dir, accession_id)

            get_gbk_cmd  = 'esearch -db nucleotide -query %s | efetch -format gb > %s' % (accession_id, gbk_file)
            os.system(get_gbk_cmd)
            cmd_list_get_gbk.append(get_gbk_cmd)

            # get organism and voucher info
            os.system('cat %s | grep "ORGANISM" > %s'         % (gbk_file, organism_txt))
            os.system('cat %s | grep "specimen_voucher" > %s' % (gbk_file, voucher_txt))

            # get sequence file
            if (os.path.isfile(gbk_file) is True) and (os.stat(gbk_file).st_size > 0):
                gbk_record = SeqIO.read(gbk_file, "genbank")
                seq_id   = gbk_record.id
                seq_desc = gbk_record.description
                seq_str  = get_origin_seq_from_gbk(gbk_file)
                with open(fa_file, 'w') as f1:
                    f1.write('>%s %s\n' % (seq_id, seq_desc))
                    f1.write('%s\n' % seq_str)

                if get_aa is True:
                    gbk_to_ffn_faa(gbk_file, ffn_file, faa_file)

        # write out command
        with open(cmd_txt_gbk, 'a') as f2:
            f2.write('\n'.join(cmd_list_get_gbk))

    ###################################### combine the results for all accessions ######################################

    tax_lineage_dict = dict()
    if tax_file is not None:
        for each_line in open(tax_file):
            each_line_split = each_line.strip().split('\t')
            tax_lineage_dict[each_line_split[0]] = each_line_split[1]

    if os.path.isfile(combined_fa) is True:
        os.system('rm %s' % combined_fa)

    os.system('cat %s/*.fasta > %s' % (tmp_dir, combined_fa))
    combined_fa_desc_handle = open(combined_fa_desc, 'w')
    if include_itol_header is True:
        combined_fa_desc_handle.write('LABELS\nSEPARATOR TAB\nDATA\n')
    for each_fa in SeqIO.parse(combined_fa, 'fasta'):
        combined_fa_desc_handle.write('%s\t%s\n' % (each_fa.id, each_fa.description.replace(' ', '_')))
    combined_fa_desc_handle.close()

    if get_aa is True:
        # get combined_faa
        if os.path.isfile(combined_faa) is True:
            os.system('rm %s' % combined_faa)
        os.system('cat %s/*.faa > %s'   % (tmp_dir, combined_faa))
        combined_faa_desc_handle = open(combined_faa_desc, 'w')
        if include_itol_header is True:
            combined_faa_desc_handle.write('LABELS\nSEPARATOR TAB\nDATA\n')
        for each_faa in SeqIO.parse(combined_faa, 'fasta'):
            combined_faa_desc_handle.write('%s\t%s\n' % (each_faa.id, each_faa.description.replace(' ', '_')))
        combined_faa_desc_handle.close()

        # get combined_ffn
        if os.path.isfile(combined_ffn) is True:
            os.system('rm %s' % combined_ffn)
        os.system('cat %s/*.ffn > %s'   % (tmp_dir, combined_ffn))
        combined_ffn_desc_handle = open(combined_ffn_desc, 'w')
        if include_itol_header is True:
            combined_ffn_desc_handle.write('LABELS\nSEPARATOR TAB\nDATA\n')
        for each_ffn in SeqIO.parse(combined_ffn, 'fasta'):
            combined_ffn_desc_handle.write('%s\t%s\n' % (each_ffn.id, each_ffn.description.replace(' ', '_')))
        combined_ffn_desc_handle.close()

    combine_esearch_op_organism(tmp_dir, tax_lineage_dict, combined_organism_info_txt)
    combine_esearch_op_voucher(tmp_dir, combined_voucher_info_txt)


if __name__ == '__main__':

    GenBank_parser = argparse.ArgumentParser(usage=GenBank_usage)
    GenBank_parser.add_argument('-i',               required=True,                       help='input txt containing accession id')
    GenBank_parser.add_argument('-o',               required=True,                       help='output dir')
    GenBank_parser.add_argument('-tax',             required=False, default=None,        help='NCBI taxonomy file, specify to get the full taxon lineage')
    GenBank_parser.add_argument('-aa',              required=False, action="store_true", help='export amino acid sequence')
    GenBank_parser.add_argument('-skip_esearch',    required=False, action="store_true", help='skip esearch')
    GenBank_parser.add_argument('-itol',            required=False, action="store_true", help='include iTOL header')
    GenBank_parser.add_argument('-f',               required=False, action="store_true", help='force overwrite')
    args = vars(GenBank_parser.parse_args())
    GenBank(args)
