import io
import sys
import pandas as pd


########################################################################################################################

rpkm_file       = '/Users/songweizhi/Desktop/cov_rpkm_stat_files/Carteriospongia1.sorted_filtered.rpkm'
allsample_stat  = '/Users/songweizhi/Desktop/allsample.stat'

'''
python3 /scratch/PI/ocessongwz/Sponge_r220/7_symbiont_abundance/get_abd_scripts/bin_rpkm_calculation.py Cliona1.sorted_filtered.rpkm allsample.stat Cliona1.cal_bin.rpkm
'''

########################################################################################################################

read_num_dict = dict()
for each_line in open('/Users/songweizhi/Desktop/allsample.stat'):
    if not each_line.startswith('file	num_seqs'):
        each_line_split = each_line.strip().split('\t')
        sample_id = each_line_split[0].split('_1.fastq')[0].split('/')[-1]
        read_num = int(each_line_split[1].replace(',', ''))
        read_num_dict[sample_id] = read_num

rpkm_out1 = pd.read_csv(rpkm_file, sep='\t', header=4)

rpkm_out2 = rpkm_out1[['#Name', 'Length', 'Reads']]
rpkm_out2['bin_name'] = rpkm_out2['#Name'].str.rsplit('_', n=1, expand=True)[0]  # contig name 'bin1_1', bin name 'bin1'
rpkm_out2['Length'] = rpkm_out2['Length'].astype(int)
rpkm_out2['Reads'] = rpkm_out2['Reads'].astype(int)

rpkm_out3 = rpkm_out2.groupby(["bin_name"])[['Length', 'Reads']].sum()

sample_name = rpkm_file.split('/')[-1].split('.sorted_filtered.rpkm')[0]
sample_read_num = read_num_dict[sample_name]

# Calulate RPKM
rpkm_out3[sample_name + '_rpkm'] = (rpkm_out3['Reads'] * 1000000000) / (rpkm_out3['Length'] * sample_read_num * 2)
print(rpkm_out3)
