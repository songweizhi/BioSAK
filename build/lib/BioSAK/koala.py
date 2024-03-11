import os
import argparse


koala_usage = '''
============= koala example commands =============

BioSAK koala -i user_ko.txt -o separated_ko_dir

==================================================
'''


def koala(args):

    file_in             = args['i']
    op_dir              = args['o']
    force_create_op_dir = args['f']

    # create op_dir
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('output directory exist, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    # read in user_ko.txt
    gnm_to_annotation_dict = dict()
    for each in open(file_in):
        each_split = each.strip().split('\t')
        gnm_id = '_'.join(each_split[0].split('_')[:-1])
        if gnm_id not in gnm_to_annotation_dict:
            gnm_to_annotation_dict[gnm_id] = set()
        else:
            gnm_to_annotation_dict[gnm_id].add(each.strip())

    # write out
    for each_gnm in gnm_to_annotation_dict:
        pwd_annotation_file = '%s/%s.txt' % (op_dir, each_gnm)
        annotation_list = sorted(list(gnm_to_annotation_dict[each_gnm]))
        pwd_annotation_file_handle = open(pwd_annotation_file, 'w')
        pwd_annotation_file_handle.write('\n'.join(annotation_list))
        pwd_annotation_file_handle.close()

    print('Done!')


if __name__ == '__main__':

    koala_parser = argparse.ArgumentParser(usage=koala_usage)
    koala_parser.add_argument('-i', required=True,                          help='user_ko.txt from BlastKOALA or GhostKOALA')
    koala_parser.add_argument('-o', required=True,                          help='output directory')
    koala_parser.add_argument('-f', required=False, action="store_true",    help='force overwrite')
    args = vars(koala_parser.parse_args())
    koala(args)
