import os
import argparse


sra2biosample_usage = '''
================= sra2biosample example commands =================

BioSAK sra2biosample -i SRR18249231 -o SRR18249231_biosample.txt
BioSAK sra2biosample -i sra.txt -o biosample.txt

==================================================================
'''


def sep_path_basename_ext(file_in):
    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]
    return f_name, f_path, f_base, f_ext


def sra_to_biosample_bioproject(sra_id_set, tmp_dir):

    # create tmp_dir
    if os.path.isdir(tmp_dir) is False:
        os.mkdir(tmp_dir)

    sra2biosample_dict = dict()
    sra2bioproject_dict = dict()
    for each_sra in sra_id_set:
        op_txt_biosample = '%s/%s_Biosample.txt' % (tmp_dir, each_sra)
        sra_to_biosample_cmd = 'esearch -db sra -query %s | efetch -format docsum | xtract -pattern DocumentSummary -element Biosample > %s' % (each_sra, op_txt_biosample)
        # print(sra_to_biosample_cmd)
        os.system(sra_to_biosample_cmd)
        with open(op_txt_biosample) as f:
            biosample_id = f.readline().strip('\n')
        sra2biosample_dict[each_sra] = biosample_id

        op_txt_bioproject = '%s/%s_Bioproject.txt' % (tmp_dir, each_sra)
        sra_to_bioproject_cmd = 'esearch -db sra -query %s | efetch -format docsum | xtract -pattern DocumentSummary -element Bioproject > %s' % (each_sra, op_txt_bioproject)
        # print(sra_to_bioproject_cmd)
        os.system(sra_to_bioproject_cmd)
        with open(op_txt_bioproject) as f:
            bioproject_id = f.readline().strip('\n')
        sra2bioproject_dict[each_sra] = bioproject_id

    return sra2biosample_dict, sra2bioproject_dict


def xml_to_txt(xml_file, txt_file):

    txt_file_handle = open(txt_file, 'w')
    for each_line in open(xml_file):
        each_line = each_line.strip()
        if (each_line.count('<') == 1) and (each_line.count('>') == 1):
            pass
        else:
            each_line = each_line.split('</')[0]
            title_str = each_line.split('>')[0][1:]
            value_str = each_line.split('>')[-1]
            if 'display_name' in title_str:
                title_str = title_str.split('display_name')[1].replace('=','').replace('"','')
            elif title_str.count('=') == 1:
                title_str = title_str.split('=')[1].replace('"','')
            elif title_str.count('=') == 0:
                title_str = title_str
            elif 'is_primary' in title_str:
                title_str = title_str.split('is_primary')[0].split('=')[1].replace('"','')
            txt_file_handle.write('%s\t%s\n' % (title_str, value_str))
    txt_file_handle.close()


def sra2biosample(args):

    file_in = args['i']
    op_txt  = args['o']

    accession_set = set()
    if os.path.isfile(file_in) is True:
        for each_id in open(file_in):
            accession_set.add(each_id.strip())
    else:
        accession_set.add(file_in)


    op_f_name, op_f_path, op_f_base, op_f_ext = sep_path_basename_ext(op_txt)

    tmp_dir = '%s/%s_sra2biosample_tmp_dir' % (op_f_path, op_f_base)

    # get sra_to_biosample_dict and sra_to_bioproject_dict
    sra_to_biosample_dict, sra_to_bioproject_dict = sra_to_biosample_bioproject(accession_set, tmp_dir)

    # write out
    op_txt_handle = open(op_txt, 'w')
    op_txt_handle.write('SRA\tBiosample\tBioproject\n')
    for each_sra_id in sorted(list(accession_set)):
        op_txt_handle.write('%s\t%s\t%s\n' % (each_sra_id, sra_to_biosample_dict.get(each_sra_id, 'na'), sra_to_bioproject_dict.get(each_sra_id, 'na')))
    op_txt_handle.close()


if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser(usage=sra2biosample_usage)
    arg_parser.add_argument('-i',   required=True,                        help='input SRA id')
    arg_parser.add_argument('-o',   required=True,                        help='output file')
    args = vars(arg_parser.parse_args())
    sra2biosample(args)
