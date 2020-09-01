import argparse


def get_ref_mean_depth(sam_depth_file):
    ref_id_list = []
    ref_len_dict = {}
    ref_depth_dict = {}
    for each_line in open(sam_depth_file):

        each_line_split = each_line.strip().split('\t')
        ref_id = each_line_split[0]
        pos = int(each_line_split[1])
        pos_depth = int(each_line_split[2])

        # get ref_id_list
        if ref_id not in ref_id_list:
            ref_id_list.append(ref_id)

        # get ref_len_dict
        if ref_id not in ref_len_dict:
            ref_len_dict[ref_id] = 1
        else:
            ref_len_dict[ref_id] += 1

        # get ref_depth_dict
        if ref_id not in ref_depth_dict:
            ref_depth_dict[ref_id] = pos_depth
        else:
            ref_depth_dict[ref_id] += pos_depth

    mean_depth_dict = {}
    for ref in sorted(ref_id_list):
        ref_len = ref_len_dict[ref]
        ref_depth = ref_depth_dict[ref]
        ref_depth_mean = float("{0:.2f}".format(ref_depth / ref_len))
        mean_depth_dict[ref] = ref_depth_mean

    return mean_depth_dict


parser = argparse.ArgumentParser()
parser.add_argument('-d', required=True, help='depth file')
args = vars(parser.parse_args())
sam_depth_file = args['d']

mean_depth_dict = get_ref_mean_depth(sam_depth_file)

for i in mean_depth_dict:
    print('%s\t%s\t%s' % (sam_depth_file, i, mean_depth_dict[i]))
