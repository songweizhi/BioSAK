import random


# num_list = [0,1,2,3,4,5,6,7,8,9]
# print(num_list)
#
# num_list_shuffle = random.shuffle(num_list)
#
#
# print(num_list)
# print(num_list_shuffle)




def get_color_list(color_num):

    color_list_to_return_sorted = [1,2,3,4,5]
    random.shuffle(color_list_to_return_sorted)
    return color_list_to_return_sorted



a = get_color_list(5)
print(a)

