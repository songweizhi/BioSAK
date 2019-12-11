
def split_with_size(input_list, size_cutoff):
    new_list = []
    temp_list = []
    sum = 0
    for each in input_list:
        if each >= size_cutoff:
            if temp_list == []:
                # append current element
                temp_list.append(each)
                new_list.append(temp_list)
                # reset temp_list
                temp_list = []
            else:
                # if temp_list != [], append temp_list first
                new_list.append(temp_list)
                # reset temp_list
                temp_list = []
                # then, append new element
                temp_list.append(each)
                new_list.append(temp_list)
                # reset temp_list and sum
                temp_list = []
                sum = 0
        else:
            sum += each
            if sum <= size_cutoff:
                temp_list.append(each)
                pass
            else:
                new_list.append(temp_list)
                temp_list = []
                temp_list.append(each)
                sum = each
    new_list.append(temp_list)
    return new_list


num_list = [1, 200, 100, 100, 190, 3, 4, 600, 90, 90, 100, 30, 40, 1, 500, 200, 100, 100, 100, 200, 1, 2, 3, 4, 5]

print(split_with_size(num_list, 200))
