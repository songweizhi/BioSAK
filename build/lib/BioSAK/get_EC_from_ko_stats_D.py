import os

ec_list_file =          '/Users/songweizhi/MetaCyc_demo/Ecoli_ko_stats_D_GeneNumber.txt'
ec_list_file_out_tmp =  '/Users/songweizhi/MetaCyc_demo/Ecoli_ec_tmp.txt'
ec_list_file_out =      '/Users/songweizhi/MetaCyc_demo/Ecoli_ec.txt'

identified_ec_list = set()
ec_list_file_out_handle = open(ec_list_file_out_tmp, 'w')
for each_line in open(ec_list_file):
    if '[EC:' in each_line:
        ec_str = each_line.strip().split('[EC:')[1].split(']')[0]

        ec_list = [ec_str]
        if ' ' in ec_str:
            ec_list = ec_str.split(' ')
        for ec in ec_list:
            print(ec)
            ec_list_file_out_handle.write(ec + '\n')
ec_list_file_out_handle.close()

os.system('cat %s | sort | uniq > %s' % (ec_list_file_out_tmp, ec_list_file_out))
os.system('rm %s' % ec_list_file_out_tmp)