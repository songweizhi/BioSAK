import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-i', required=True, help='file in')
parser.add_argument('-o', required=True, help='file out')
parser.add_argument('-s', required=True, help='read strand, 1 or 2')

args = vars(parser.parse_args())
file_in     = args['i']
file_out    = args['o']
read_strand = args['s']


file_out_handle = open(file_out, 'w')
for each in open(file_in):

    if each.startswith('>'):
        each_split = each.strip().split(' ')
        print(each_split)
        file_out_handle.write('%s/%s\n' % ('_'.join(each_split[0].split(':')), read_strand))
    else:
        file_out_handle.write(each)
file_out_handle.close()

