
########## use * to grab excess items ##########

a1, b1, *rest1 = [1, 2, 3, 4, 5]
a2, *rest2, b2 = [1, 2, 3, 4, 5]
a3, *rest3, b3 = [1, 2, 3]
a4, *rest4, b4 = [1, 2]

print(rest1)    # [3, 4, 5]
print(rest2)    # [2, 3, 4]
print(rest3)    # [2]
print(rest4)    # []


########## swapping the values of variables without using a temporary variable ##########

a = 'apple'
b = 'orange'
a, b = b, a

print(a)    # orange
print(b)    # apple



