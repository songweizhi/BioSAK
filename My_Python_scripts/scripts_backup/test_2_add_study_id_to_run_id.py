import os
# get run to study dict

selected = 'antarctic_metagenome_output_selected.txt'
run2study_dict = {}
for each in open(selected):
    #print(each)
    each_split = each.strip().split('\t')
    study_id = each_split[0]
    run_id = each_split[1]
    run2study_dict[run_id] = study_id

print(run2study_dict)

for each_run in open('antarctic_metagenome_output_run_ids.txt'):
    each_run = each_run.strip()
    #print(each_run)
    #print('mv %s_Q30.fasta %s_%s.fasta' % (each_run, run2study_dict[each_run], each_run))
    os.system('mv %s_Q30.fasta %s_%s.fasta' % (each_run, run2study_dict[each_run], each_run))


