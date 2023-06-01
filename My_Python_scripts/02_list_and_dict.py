
######################################################### list #########################################################

import random

num_list = [0,1,2,3,4,5,6,7,8,9]
num_list_randomized = random.sample(num_list, len(num_list))
print(num_list_randomized)
# [7, 8, 0, 5, 2, 4, 1, 9, 3, 6]


#################### unique_list_elements ####################

def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output

print(unique_list_elements([1,2,3,3,3,3,3,4,5, 'a', 'a']))
# [1, 2, 3, 4, 5, 'a']


#################### get_intersection ####################

a = [1,2,3,4,'a',5,6]
b = [3,4,5,6,'a', 'b',7,8]
c = set(a).intersection(b)
print(c)
# {3, 4, 5, 6, 'a'}


#################### enumerate ####################

for i, j in enumerate(['a', 'b', 'c', 'd']):
    print(i, j)
# 0 a
# 1 b
# 2 c
# 3 d


#################### for loop in one line ####################

colors = ['black', 'white']
sizes = ['S', 'M', 'L']
t_shirts = [(color, size) for color in colors for size in sizes]
print(t_shirts)
# [('black', 'S'), ('black', 'M'), ('black', 'L'), ('white', 'S'), ('white', 'M'), ('white', 'L')]


######################################################### dict #########################################################

#################### merge_two_dict ####################

def merge_two_dict(dict_a, dict_b):
    dict_c = dict_a.copy()
    dict_c.update(dict_b)
    return dict_c


dict_a = {'A': 'aa', 'B': 'bb'}
dict_b = {'C': 'cc', 'D': 'dd'}
dict_c = merge_two_dict(dict_a, dict_b)
print(dict_c)
# {'A': 'aa', 'B': 'bb', 'C': 'cc', 'D': 'dd'}



