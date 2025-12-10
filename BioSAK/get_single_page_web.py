import argparse
import os

get_single_page_web_usage = '''
================================ get_single_page_web example commands ================================

BioSAK get_single_page_web -d doc/images -m doc/metadata.txt -w doc/template.html -o single_page_web.html
BioSAK get_single_page_web -d doc/images -m doc/metadata.txt -w doc/template_with_barchart.html -o single_page_web_with_barchart.html

# Note
There should be no parentheses in metadata

======================================================================================================
'''


def get_single_page_web(args):

    sample_metadata_txt    = args['m']
    web_page_file_template = args['w']
    photo_dir              = args['d']
    web_page_file_out      = args['o']

    # read in sample metadata
    photo_section_str = ''
    col_index = dict()
    line_num_index = 0
    for each_line in open(sample_metadata_txt):
        line_num_index += 1
        line_split = each_line.strip().split('\t')
        if line_num_index == 1:
            col_index = {key: i for i, key in enumerate(line_split)}
        else:
            sample_accession    = line_split[col_index['Accession']]
            sample_image        = line_split[col_index['Image']]
            sample_taxon        = line_split[col_index['Taxon']]
            sample_location     = ''
            sample_location_eng = line_split[col_index['Location_Eng']]
            sample_location_chi = line_split[col_index['Location_Chi']]
            sample_latitude     = line_split[col_index['Latitude']]
            sample_longitude    = line_split[col_index['Longitude']]
            sample_source       = line_split[col_index['Source']]
            sample_desc         = line_split[col_index['Description']]
            sample_date         = line_split[col_index['Date']]
            sample_link         = line_split[col_index['Link']]
            loci_with_coord     = '%s (%s/%s)' % (sample_location, sample_latitude, sample_longitude)
            sample_coord        = '%s/%s' % (sample_latitude, sample_longitude)

            pwd_photo = sample_image
            if photo_dir is not None:
                pwd_photo = '%s/%s' % (photo_dir, sample_image)

            photo_attribute_str = '      {accession: "%s", photo: "%s", taxon: "%s", loci: "%s", loci_eng: "%s", loci_chi: "%s", lat: %s, lng: %s, desc: "%s", coord: "%s", source: "%s", date: "%s", link: "%s"},' % (sample_accession, pwd_photo, sample_taxon, sample_location, sample_location_eng, sample_location_chi, sample_latitude, sample_longitude, sample_desc, sample_coord, sample_source, sample_date, sample_link)
            photo_section_str += photo_attribute_str
            photo_section_str += '\n'

    # write out
    web_page_file_out_handle = open(web_page_file_out, 'w')
    for each_line in open(web_page_file_template):
        if '// add photo attributes here' not in each_line:
            web_page_file_out_handle.write(each_line)
        else:
            web_page_file_out_handle.write(photo_section_str)
    web_page_file_out_handle.close()


if __name__ == '__main__':
    get_single_page_web_parser = argparse.ArgumentParser(usage=get_single_page_web_usage)
    get_single_page_web_parser.add_argument('-m', required=True,                help='sample metadata')
    get_single_page_web_parser.add_argument('-w', required=True,                help='webpage template')
    get_single_page_web_parser.add_argument('-d', required=False, default=None, help='photo directory')
    get_single_page_web_parser.add_argument('-o', required=True,                help='output webpage')
    args = vars(get_single_page_web_parser.parse_args())
    get_single_page_web(args)
