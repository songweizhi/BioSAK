import os
import shutil

def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


########################################################################################################################

combined_checkm_quality_file = '/Users/songweizhi/Desktop/Mariana_Trench_checkm.txt'
completeness_cutoff = 50
contamination_cutoff = 5
bin_folder_in = '/Users/songweizhi/Desktop/Mariana_Trench_refined_bins_all'


########################################################################################################################


infile_path, infile_basename, infile_extension = sep_path_basename_ext(combined_checkm_quality_file)
file_out = '%s/%s_cmplt%s_contam%s%s' % (infile_path, infile_basename, str(completeness_cutoff), str(contamination_cutoff), infile_extension)


bin_folder_out = '%s_%s_%s' % (bin_folder_in, completeness_cutoff, contamination_cutoff)
force_create_folder(bin_folder_out)


filtered_checkm_quality_handle = open(file_out, 'w')
filtered_checkm_quality_handle.write('Genome_ID\tCompleteness\tContamination\n')
qualified_bin_num = 0
for each in open(combined_checkm_quality_file):

    if not ((each.startswith('---')) or (each.startswith('  Bin Id'))):
        quality_split = each.strip().split(' ')
        quality_split_new = []
        for each in quality_split:
            if each != '':
                quality_split_new.append(each)

        genome_id = quality_split_new[0]
        completeness = float(quality_split_new[12])
        contamination = float(quality_split_new[13])

        if (completeness >= completeness_cutoff) and (contamination <= contamination_cutoff):
            filtered_checkm_quality_handle.write('%s\t%s\t%s\n' % (genome_id, completeness, contamination))
            #os.system('cp %s/%s.fasta %s/' % (bin_folder_in, genome_id, bin_folder_out))
            qualified_bin_num += 1

filtered_checkm_quality_handle.close()

print('The number of bins with completeness >= %s and contamination <= %s: %s' % (completeness_cutoff, contamination_cutoff, qualified_bin_num))

