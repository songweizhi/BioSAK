
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


sam_depth_file = '/Users/songweizhi/Desktop/OneHGT_50_to_50.depth'
sam_depth_file = '/Users/songweizhi/Desktop/OneHGT_75_to_25.depth'
sam_depth_file = '/Users/songweizhi/Desktop/OneHGT_90_to_10.depth'
sam_depth_file = '/Users/songweizhi/Desktop/OneHGT_95_to_05.depth'
sam_depth_file = '/Users/songweizhi/Desktop/OneHGT_99_to_01.depth'

print(get_ref_mean_depth(sam_depth_file))

# 50_to_50：{'2.10_chromosome': 55.29, 'D2_c': 45.67}
# 75_to_25：{'2.10_chromosome': 82.91, 'D2_c': 22.84}
# 90_to_10：{'2.10_chromosome': 99.47, 'D2_c': 9.14}
# 95_to_05: {'2.10_chromosome': 105.0, 'D2_c': 4.62}
# 99_to_01: {'2.10_chromosome': 109.41, 'D2_c': 1.53}
