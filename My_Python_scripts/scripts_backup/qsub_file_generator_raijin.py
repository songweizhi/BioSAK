import os
import glob
import shutil


###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 14
memory = 120
walltime_needed = '47:59:00'
email = 'wythe1987@163.com'
#email = 'weizhi.song@student.unsw.edu.au'
modules_needed = [/3.9.0']

wd = '/Users/songweizhi/Desktop'
outputs_folder = 'qsub_6_idba_ud'

###########################################################################################

os.chdir(wd)

# create outputs folder
if not os.path.exists(outputs_folder):
    os.makedirs(outputs_folder)
else:
    shutil.rmtree(outputs_folder)
    os.makedirs(outputs_folder)


# Prepare header
line_1 = '#!/bin/bash\n'
line_2 = '#PBS -P du5\n'
line_3 = '#PBS -q hugemem\n'
line_4 = '#PBS -l ncpus=14,walltime=11:59:00,mem=300GB\n'
line_5 = '#PBS -l wd\n'
line_6 = '#PBS -M wythe1987@163.com\n'
line_7 = '#PBS -m ae\n\n'
header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7

# Prepare module lines
module_lines = ''
for module in modules_needed:
    module_lines += 'module load ' + module + '\n'

print(header)
print(module_lines)


################################################################

for each_run in open('tests.txt'):
    each_run = each_run.strip()
    output_file_handle = open('%s/qsub_%s.sh' % (outputs_folder, each_run), 'w')
    #print(each_run)
    # if len(each_run) == 9:
    #     each_run_no_last_3_4 = each_run[:-3]
    # if len(each_run) == 10:
    #     each_run_no_last_3_4 = each_run[:-4]

    #print(each_run_no_last_3)
    #print(header)
    #print(module_lines)
    #print('cd /srv/scratch/z5039045/Antarctica')
    #print('wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/%s/%s/%s.sra' % (each_run_no_last_3_4, each_run, each_run))
    output_file_handle.write(header)
    #output_file_handle.write(module_lines)
    output_file_handle.write('cd /short/du5/wzs561/MetaCHIP/combined_idba_reads\n')
    #output_file_handle.write('fq2fa %s_Q30.fastq %s_Q30.fasta\n' % (each_run, each_run))
    #output_file_handle.write('spades.py --meta --only-assembler -t %s -s combined_%s_reads.fasta -o combined_%s_k21-127 -k 21,33,55,77,99,127' % (ppn_number, each_run, each_run))
    output_file_handle.write('/short/du5/wzs561/idba-master/bin/idba_ud --pre_correction --num_threads 14 --mink 20 --maxk 124 --step 20 --read %s --out combined_%s_k20-124' % (each_run, each_run))



    #output_file_handle.write('wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/%s/%s/%s.sra\n' % (each_run_no_last_3_4, each_run, each_run))

    output_file_handle.close()




# spades.py --meta --only-assembler -t 6 -s combined_reads_m5_l100i250.fasta -o combined_reads_m5_l100i250_k21-125 -k 21,33,55,77,99,127

# /short/du5/wzs561/idba-master/bin/idba_ud --pre_correction --num_threads 14 --mink 20 --maxk 124 --step 20 --read GemSIM_100bp_m0.fasta --out GemSIM_100bp_m0_k20-124
