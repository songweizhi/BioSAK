import os
import shutil


##################################################### CONFIGURATION ####################################################

# parameters for generate_qsub_file
nodes_number = 1
ppn_number = 3
memory = 20
walltime_needed = '11:59:00'
email = 'songwz03@163.com'
modules_needed = ['python/3.5.2', 'blast+/2.6.0', 'diamond/0.9.24']
qsub_file_folder = '/Users/songweizhi/Desktop/qsub_files_KEGG'
wd_on_katana = '/srv/scratch/z5039045/MetaCHIP_TT_90MGs/GoodBins_0.5_0.05_MetaCHIP_wd/GoodBins_0.5_0.05_all_faa_files'


################################################## generate_qsub_file ##################################################


# Create qsub file folder
if os.path.isdir(qsub_file_folder):
    shutil.rmtree(qsub_file_folder)
    if os.path.isdir(qsub_file_folder):
        shutil.rmtree(qsub_file_folder)
        if os.path.isdir(qsub_file_folder):
            shutil.rmtree(qsub_file_folder)
os.system('mkdir %s' % qsub_file_folder)


# Prepare header
line_1 = '#!/bin/bash'
line_2 = '#PBS -l nodes=%s:ppn=%s' % (str(nodes_number), str(ppn_number))
line_3 = '#PBS -l vmem=%sgb' % str(memory)
line_4 = '#PBS -l walltime=%s' % walltime_needed
line_5 = '#PBS -j oe'
line_6 = '#PBS -M %s' % email
line_7 = '#PBS -m ae'
header = '%s\n%s\n%s\n%s\n%s\n%s\n%s\n' % (line_1, line_2, line_3, line_4, line_5, line_6, line_7)


# Prepare module lines
module_lines = ''
for module in modules_needed:
    module_lines += 'module load %s\n' % module


# write to qsub files


for each_faa in open('/Users/songweizhi/Desktop/faa_file_list.txt'):

    each_faa_basename = '.'.join(each_faa.strip().split('.')[:-1])
    print(each_faa_basename)


    current_qsub_file = '%s/qsub_KEGG_%s.sh' % (qsub_file_folder, each_faa_basename)
    handle = open(current_qsub_file, 'w')
    handle.write('\n' + header + '\n')
    handle.write(module_lines)
    handle.write('\ncd %s\n' % wd_on_katana)
    KEGG_cmd = 'python3 /srv/scratch/z5039045/Scripts/KEGG_wrapper.py -seq_in %s -p %s -t %s -DB_dir /srv/scratch/z5039045/KEGG_DB -diamond' % (each_faa.strip(), each_faa_basename, ppn_number)
    handle.write(KEGG_cmd)
    handle.close()




