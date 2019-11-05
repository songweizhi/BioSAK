# CAZyDB_fam_activities = '/Users/songweizhi/Desktop/CAZyDB.07312019.fam-activities.txt'
#
# fam_to_activities_dict = {}
# for each_fam in open(CAZyDB_fam_activities):
#     each_fam_split = each_fam.strip().split('	  ')
#     if len(each_fam_split) == 2:
#         fam_id = each_fam_split[0]
#         fam_activities = each_fam_split[1]
#         fam_to_activities_dict[fam_id] = fam_activities
#
#
# print('Query\tFamily\tActivities')
# hmm_to_num_dict = {}
# for hmm_hit in open('/Users/songweizhi/Desktop/Refined_3.out.dm.ps.stringent'):
#     hmm_hit_split = hmm_hit.strip().split('\t')
#     query_id = hmm_hit_split[2]
#     matched_hmm = hmm_hit_split[0]
#
#     matched_hmm_id = matched_hmm.split('.hmm')[0]
#     # if '_' in matched_hmm_id:
#     #     matched_hmm_id = matched_hmm_id.split('_')[0]
#
#     matched_hmm_activities = 'NA'
#     if matched_hmm_id in fam_to_activities_dict:
#         matched_hmm_activities = fam_to_activities_dict[matched_hmm_id]
#
#     if matched_hmm not in hmm_to_num_dict:
#         hmm_to_num_dict[matched_hmm] = 1
#     else:
#         hmm_to_num_dict[matched_hmm] += 1
#
#
#     print('%s\t%s\t%s' % (query_id, matched_hmm, matched_hmm_activities))
#
#



list_1 = ['A','B','C','D' ]
print(list_1)
list_2 = list_1[::-1]
#print(list_2)
print(list_1)