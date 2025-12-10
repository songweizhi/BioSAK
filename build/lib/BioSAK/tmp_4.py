import os
import glob


def summarise_metadata(file_dir, file_ext, summary_txt):


    file_re = '%s/*.%s' % (file_dir, file_ext)
    file_list = glob.glob(file_re)

    for each_file in file_list:
        file_name = os.path.basename(each_file)
        sample_id = file_name[:-(len(file_ext) + 1)]
        current_isolate = 'na'
        current_host = 'na'
        for each_line in open(each_file):
            if '/host="' in each_line:
                current_host = each_line.split('"')[1]
            elif '/isolate="' in each_line:
                current_isolate = each_line.split('"')[1]

        print('%s\t%s\t%s' % (sample_id, current_isolate, current_host))



file_dir    = '/Users/songweizhi/Desktop/biosamples_metadata/tmp'
file_ext    = 'txt'
summary_txt = '/Users/songweizhi/Desktop/biosamples_metadata/tmp.txt'


summarise_metadata(file_dir, file_ext, summary_txt)