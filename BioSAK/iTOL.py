import os
import sys
import argparse
from itolapi import Itol, ItolExport


iTOL_usage = '''
=========================================== iTOL example commands ===========================================

BioSAK iTOL -tree Demo.tree -label Demo_label.txt -color Demo_color.txt -collapse Demo_collapse.txt

# Note:
1. iTOL module needs internet connection for working.
2. Some example input files can be found from the tutorial folder in BioSAK GitHub repository.
3. iTOL module doesn't show collapsed nodes. You should go to the provided link for better visualization. 

=============================================================================================================
'''

'''
# vertical_shift_factor
# horizontal_scale_factor
# leaf_sorting (disable Leaf sorting, important!!!)
'''

def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def iTOL(args):

    # read in arguments
    tree_file =         args['tree']
    collapse_file =     args['collapse']
    label_file =        args['label']
    color_file =        args['color']

    # define output file name
    tree_file_path, tree_file_basename, tree_file_extension = sep_path_basename_ext(tree_file)
    pwd_output_pdf = '%s/%s.pdf' % (tree_file_path, tree_file_basename)

    # Create the Itol class and add tree files
    my_job = Itol()
    my_job.add_file(tree_file)

    if color_file is not None:
        my_job.add_file(color_file)

    if label_file is not None:
        my_job.add_file(label_file)

    if collapse_file is not None:
        my_job.add_file(collapse_file)

    # upload tree
    print('Uploading tree. This may take some time depending on how large the tree is and how much load there is on the itol server')
    good_upload = my_job.upload()
    if not good_upload:
        print('There was an error:' + my_job.comm.upload_output)
        sys.exit(1)

    print('Tree Web Page URL: %s'     % my_job.get_webpage())

    # Export the tree above to pdf
    itol_exporter = my_job.get_itol_export()
    itol_exporter.set_export_param_value('format', 'pdf')
    itol_exporter.set_export_param_value('datasetList', 'dataset1')
    itol_exporter.export(pwd_output_pdf)

    print('Tree plot exported to %s.pdf' % tree_file_basename)


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()

    parser.add_argument('-tree',      required=True,  help='tree file in newick format')
    parser.add_argument('-label',     required=True,  help='label')
    parser.add_argument('-color',     required=False, help='color')
    parser.add_argument('-collapse',  required=False, help='collapse)')

    args = vars(parser.parse_args())

    iTOL(args)
