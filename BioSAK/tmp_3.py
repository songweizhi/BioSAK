

# for i in range(1, 73):
#
#     cmd = 'srun -n 1 python3 00_DataNeeded/BLCA/2.blca_main.py -x -p 1 -r 00_DataNeeded/db_GTDB214.1/ssu_all_r214_BLCAparsed.taxonomy -q 00_DataNeeded/db_GTDB214.1/ssu_all_r214_BLCAparsed.fasta -i b_rep_set_%s.fa -o b_rep_set_%s_BLCA_out.1.txt &' % (i, i)
#     print(cmd)


def ko_stats_to_txt(ko_stats_dict, ko_desc_dict, op_txt, op_stats_txt):
    op_txt_handle = open(op_txt, 'w')
    op_stats_txt_handle = open(op_stats_txt, 'w')
    for ko_high in sorted(list(ko_stats_dict.keys())):
        ko_d_set = ko_stats_dict[ko_high]
        ko_high_desc = ko_desc_dict[ko_high]
        op_stats_txt_handle.write('%s\t%s\t%s\n' % (ko_high, len(ko_d_set), ko_high_desc))
        for ko_d in sorted(list(ko_d_set)):
            ko_d_desc = ko_desc_dict[ko_d]
            op_txt_handle.write('%s\t%s\t%s\t%s\n' % (ko_high, ko_high_desc, ko_d, ko_d_desc))
    op_txt_handle.close()
    op_stats_txt_handle.close()

print('aaa')

