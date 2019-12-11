
wd = '/Users/songweizhi/Desktop/new'

transfer_file = open('%s/donor2recip.txt' % wd)
negatives = open('%s/HGT_candidates_ng_with_direction.txt' % wd)
#predicted = open('%s/HGT_candidates_100.txt' % wd)
predicted = open('%s/HGT_candidates_100_with_direction.txt' % wd)


donor_list = []
for each in transfer_file:
    each_split = each.strip().split('\t')
    donor = each_split[0]
    if donor != 'donor_gene':
        donor_list.append(donor)

print('donor_list(%s):' % len(donor_list))
print('\t'.join(donor_list))


neg_A_to_B_list = []
neg_B_to_A_list = []
for each_neg in negatives:
    #print(each_neg)
    each_neg_split = each_neg.strip().split('\t')
    if each_neg.startswith('A'):
        neg_B_to_A_list.append(each_neg_split[0])
    if each_neg.startswith('B'):
        neg_A_to_B_list.append(each_neg_split[1])

# neg_A_to_B_list_new = []
# for each in neg_A_to_B_list:
#     if each not in set(neg_A_to_B_list).intersection(neg_B_to_A_list):
#         neg_A_to_B_list_new.append(each)
#
# neg_B_to_A_list_new = []
# for each in neg_B_to_A_list:
#     if each not in set(neg_A_to_B_list).intersection(neg_B_to_A_list):
#         neg_B_to_A_list_new.append(each)


print('\nneg_A_to_B_list(%s):' % len(neg_A_to_B_list))
print('\t'.join(neg_A_to_B_list))

# print('\nneg_A_to_B_list_new(%s):' % len(neg_A_to_B_list_new))
# print('\t'.join(neg_A_to_B_list_new))

print('\nneg_B_to_A_list(%s):' % len(neg_B_to_A_list))
print('\t'.join(neg_B_to_A_list))

# print('\nneg_B_to_A_list_new(%s):' % len(neg_B_to_A_list_new))
# print('\t'.join(neg_B_to_A_list_new))

print('\nintersection negative group(%s)' % len(set(neg_A_to_B_list).intersection(neg_B_to_A_list)))
print('\t'.join(set(neg_A_to_B_list).intersection(neg_B_to_A_list)))

predicted_A_to_B_list = []
predicted_B_to_A_list = []
for each_hgt in predicted:
    each_hgt_split = each_hgt.strip().split('\t')
    if each_hgt.startswith('A'):
        predicted_B_to_A_list.append(each_hgt_split[0])
    if each_hgt.startswith('B'):
        predicted_A_to_B_list.append(each_hgt_split[1])
print('\npredicted_A_to_B_list(%s):' % len(predicted_A_to_B_list))
print('\t'.join(predicted_A_to_B_list))
print('\npredicted_B_to_A_list(%s):' % len(predicted_B_to_A_list))
print('\t'.join(predicted_B_to_A_list))
print('\nintersection predicted group(%s)' % len(set(predicted_A_to_B_list).intersection(predicted_B_to_A_list)))
print('\t'.join(set(predicted_A_to_B_list).intersection(predicted_B_to_A_list)))

validated_A_to_B_list = []
for each in predicted_A_to_B_list:
    if each in donor_list:
        validated_A_to_B_list.append(each)
print('\nvalidated_A_to_B_list(%s):' % len(validated_A_to_B_list))
print('\t'.join(validated_A_to_B_list))


validated_B_to_A_list = []
for each in predicted_B_to_A_list:
    if each in donor_list:
        validated_B_to_A_list.append(each)
print('\nvalidated_B_to_A_list(%s):' % len(validated_B_to_A_list))
print('\t'.join(validated_B_to_A_list))


detected_in_both_direction = []
not_predicted_list = []
for each in donor_list:
    if (each in validated_A_to_B_list) and (each in validated_B_to_A_list):
        detected_in_both_direction.append(each)
    if (each not in predicted_A_to_B_list) and (each not in predicted_B_to_A_list):
        not_predicted_list.append(each)

print('\nvalidated_in_both_direction(%s):' % len(detected_in_both_direction))
print('\t'.join(detected_in_both_direction))
print('\nnot_predicted_list(%s):' % len(not_predicted_list))
print('\t'.join(not_predicted_list))