

def scale_str_to_size_list(scale_str):

    scale_list = scale_str.split('-')
    scale_list = [float(i) for i in scale_list]

    shape_size_list = []
    if scale_list[0] == 0:
        shape_size_list = [0]
        for each_value in scale_list[1:-1]:
            current_size = each_value/scale_list[-1]
            shape_size_list.append(current_size)
        shape_size_list.append(1)

    if scale_list[0] != 0:
        shape_size_list = [0.1]
        interval_num = len(scale_list) - 1
        interval_value = (1 - 0.1)/interval_num
        n = 1
        for each_value in scale_list[1:-1]:
            shape_size_list.append(interval_value * n + 0.1)
            n += 1
        shape_size_list.append(1)

    return shape_size_list

