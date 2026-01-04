import os
import argparse
import subprocess


metaBiosample_usage = '''
=========== metaBiosample example commands ===========

Dependency: efetch

BioSAK metaBiosample -f -i SAMEA5126984 -o op_dir
BioSAK metaBiosample -f -i biosamples.txt -o op_dir

=====================================================
'''

def check_executables(program_list):

    not_detected_programs = []
    for needed_program in program_list:
        if subprocess.call(['which', needed_program], stdout=open(os.devnull, 'wb')) != 0:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not detected, program exited!' % ','.join(not_detected_programs))
        exit()


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def efetch_op_to_dict(efetch_op_txt):

    metadata_dict = dict()
    for each_line in open(efetch_op_txt):
        each_line = each_line.strip()
        if not each_line.endswith(':'):
            if each_line.startswith('/'):
                each_line = each_line[1:]
                each_line_split = each_line.split('=')
                attribute_name = each_line_split[0]
                attribute_value = each_line_split[1][1:-1]
                metadata_dict[attribute_name] = attribute_value
            elif each_line.startswith('Identifiers: '):
                each_line = each_line[len('Identifiers: '):]
                each_line_split = each_line.split(';')
                for each_identifier in each_line_split:
                    each_identifier = each_identifier.strip()
                    each_identifier_split = each_identifier.split(': ')
                    metadata_dict[each_identifier_split[0]] = each_identifier_split[1]
            elif each_line.startswith('1: '):
                desc_line = each_line[len('1: '):]
                metadata_dict['sample_description'] = desc_line
            elif each_line.startswith('Organism: '):
                desc_line = each_line[len('Organism: '):]
                metadata_dict['organism'] = desc_line
            else:
                pass

    return metadata_dict


def metaBiosample(args):

    file_in         = args['i']
    op_dir          = args['o']
    force_overwrite = args['f']

    check_executables(['efetch'])

    tmp_dir               = '%s/tmp'            % op_dir
    combined_metadata_txt = '%s/metadata.txt'   % op_dir

    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('specified output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % tmp_dir)

    accession_set = set()
    if os.path.isfile(file_in) is True:
        for each_id in open(file_in):
            accession_set.add(each_id.strip())
    else:
        accession_set.add(file_in)

    efetch_cmd_list = []
    op_file_set = set()
    for each_id in accession_set:
        biosample_id = each_id.strip()
        op_meta_txt = '%s/%s.txt' % (tmp_dir, biosample_id)
        efetch_cmd = 'efetch -db biosample -id %s -format abstract > %s' % (biosample_id, op_meta_txt)
        efetch_cmd_list.append(efetch_cmd)
        op_file_set.add(op_meta_txt)

    cmd_txt = '%s/cmd.txt' % op_dir
    with open(cmd_txt, 'w') as cmd_txt_handle:
        cmd_txt_handle.write('\n'.join(sorted(efetch_cmd_list)) + '\n')

    for each_cmd in efetch_cmd_list:
        os.system(each_cmd)

    all_attr_set = set()
    metadata_dod = dict()
    for each_file in op_file_set:
        f_name, f_path, f_base, f_ext = sep_path_basename_ext(each_file)
        current_metadata_dict = efetch_op_to_dict(each_file)
        for i in current_metadata_dict:
            all_attr_set.add(i)
        metadata_dod[f_base] = current_metadata_dict

    all_attr_list_sorted = sorted(list(all_attr_set))

    combined_metadata_txt_handle = open(combined_metadata_txt, 'w')
    combined_metadata_txt_handle.write('Biosample\t%s\n' % ('\t'.join(all_attr_list_sorted)))
    for each_biosample in sorted(list(metadata_dod.keys())):
        current_biosample_attr_list = [each_biosample]
        for each_attr in all_attr_list_sorted:
            current_biosample_attr_list.append(metadata_dod[each_biosample].get(each_attr, 'na'))
        combined_metadata_txt_handle.write('\t'.join(current_biosample_attr_list) + '\n')
    combined_metadata_txt_handle.close()


if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser(usage=metaBiosample_usage)
    arg_parser.add_argument('-i',   required=True,                        help='biosample id file, one id per line')
    arg_parser.add_argument('-o',   required=True,                        help='output folder')
    arg_parser.add_argument('-f',   required=False, action="store_true",  help='Force overwrite existing results')
    args = vars(arg_parser.parse_args())
    metaBiosample(args)
