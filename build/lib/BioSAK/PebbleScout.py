import os
import glob
import csv
import xml.sax
import argparse
import requests
import xmltodict
import pandas as pd
from Bio import SeqIO
import multiprocessing as mp


PebbleScout_usage = '''
======================== PebbleScout example commands ========================

BioSAK PebbleScout -f -m -t 12 -db meta,meta_vol2 -i gnm_dir -x fna -o op_dir

==============================================================================
'''

########################################################################################################################
############################################### biosample section start ################################################
########################################################################################################################

class BioSamplesParser(xml.sax.ContentHandler):
    '''
    SAX parser to extract information from BioSamples XML file and store it in a dictionary.
    '''
    def __init__(self, content_dict) -> None:
        super().__init__()
        self.cur_dict = content_dict
        self.biosample_id = ''

        self.is_db = False
        self.cur_db = ''
        self.is_title = False
        self.paragraph_string = ''
        self.is_paragraph = False
        self.is_owner = False
        self.is_package = False
        self.attribute_name = ''
        self.is_link = False
        self.link_name = ''
        self.owner_xml_string = '<Owner>'

    def unpack_owner_dict(self, xml_dict: dict, parent_key: str = '', output_list: list = None) -> list:
            '''
            Unpack a given owner dictionary into a single list that has all values from the dictionary.
            The list contains tuples of the form (key, value).
            If the dictionary has #text as a key, store the other attributes in the key as well.
            '''
            if output_list is None:
                output_list = []
            for k, v in xml_dict.items():
                key = parent_key + '/' + k
                if type(v) == dict:
                    self.unpack_owner_dict(v, key, output_list)
                elif type(v) == list and len(v) > 0 and type(v[0]) == dict:
                    key_copy = key
                    for item in v:
                        if '#text' in item.keys():
                            key = key_copy
                            for k2, v2 in item.items():
                                if k2 != '#text':
                                    key = key + '/' + k2 + '=' + v2
                            output_list.append((key, item['#text']))
                        else:
                            self.unpack_owner_dict(item, key, output_list)
                else:
                    output_list.append((key, v))

            return output_list

    def startElement(self, name, attrs):
        if name == 'BioSample':
            self.biosample_id == attrs['accession']
            self.cur_dict['submission_date'] = attrs['submission_date'] if 'submission_date' in attrs else ''
            self.cur_dict['last_update'] = attrs['last_update'] if 'last_update' in attrs else ''
            self.cur_dict['publication_date'] = attrs['publication_date'] if 'publication_date' in attrs else ''
            self.cur_dict['access'] = attrs['access'] if 'access' in attrs else ''
        elif name == 'Id':
            self.is_db = True
            try:
                self.cur_db = attrs['db']
            except KeyError:
                self.cur_db = attrs['db_label']
        elif name == 'Title':
            self.is_title = True
        elif name == 'Organism':
            self.cur_dict['organism_taxonomy_id'] = attrs['taxonomy_id']
            self.cur_dict['organism_taxonomy_name'] = attrs['taxonomy_name']
        elif name == 'Paragraph':
            self.is_paragraph = True
        elif name == 'Owner':
            self.is_owner = True
        elif self.is_owner:
            self.owner_xml_string += '<{}'.format(name)
            if len(attrs) > 0:
                for attr in attrs.getNames():
                    self.owner_xml_string += ' {}="{}"'.format(attr, attrs[attr])
            self.owner_xml_string += '>'
        elif name == 'Package':
            self.is_package = True
        elif name == 'Attribute':
            try:
                self.attribute_name = attrs['harmonized_name']
            except KeyError:
                self.attribute_name = attrs['attribute_name']
        elif name == 'Link':
            # assumes format for <Link> tag has either url or entrez attribute
            if 'type' in attrs and attrs['type'] == 'url' and 'label' in attrs:
                # doesn't work well for GEO sample labels
                self.link_name = attrs['label']
            elif 'type' in attrs and attrs['type'] == 'entrez':
                self.link_name = attrs['target']
            self.is_link = True
        elif name == 'Status':
            self.cur_dict['status'] = attrs['status']
            self.cur_dict['when'] = attrs['when']

    def characters(self, content):
        if self.is_title:
            self.cur_dict['title'] = content
            self.is_title = False
        elif self.is_db:
            self.cur_dict[self.cur_db + '_db'] = content
            self.is_db = False
        elif self.is_paragraph:
            self.paragraph_string += content
            self.is_paragraph = False
        elif self.is_owner:
            self.owner_xml_string += content
        elif self.is_package:
            self.cur_dict['package'] = content
            self.is_package = False
        elif self.attribute_name:
            self.cur_dict[self.attribute_name] = content
            self.attribute_name = ''
        elif self.is_link:
            self.cur_dict[self.link_name] = content
            self.link_name = ''
            self.is_link = False

    def endElement(self, name):
        if name == 'Owner':
            self.owner_xml_string += '</Owner>'
            # replace ampersand with 'and' to avoid xml parsing error
            self.owner_xml_string = self.owner_xml_string.replace('&', 'and')
            owner_dict = self.unpack_owner_dict(xml_dict=xmltodict.parse(self.owner_xml_string))
            for k, v in owner_dict:
                self.cur_dict[k[7:]] = v
            self.is_owner = False
        elif self.is_owner:
            self.owner_xml_string += '</{}>'.format(name)
        if name == 'BioSample':
            self.cur_dict['paragraph'] = self.paragraph_string
            self.paragraph_string = ''

    def endDocument(self):
        pass

def biosample(args):

    input_file = args['i']
    output_file = args['o']
    values_file = None

    # extract the values from the values file if provided
    csv_headers = ['biosample_id']
    set_values = False
    if values_file:
        set_values = True
        with open(values_file, 'r') as f:
            values = f.read().splitlines()
            csv_headers.extend(values)

    with open(input_file, 'r') as f:
        # results_dict contains a dictionary for each biosample id with the metadata
        results_dict = {}
        biosample_ids = f.read().splitlines()
        n = 1
        for accession in biosample_ids:
            print('Getting BioSample metadata %s/%s: %s' % (n, len(biosample_ids), accession))
            n += 1

            # content_dict contains the metadata for a single biosample id
            content_dict = {}
            handler = BioSamplesParser(content_dict)
            url = 'https://serratus-biosamples.s3.us-east-1.amazonaws.com/biosamples_split/' + accession + '.xml'
            headers = {'Host': 'serratus-biosamples.s3.amazonaws.com'}
            # get the xml file from the serratous-biosamples s3 bucket
            response = requests.get(url, headers=headers)
            if response.status_code != 200:
                print('Error: status code {} for {}'.format(response.status_code, accession))
                continue
            result_string = response.content.decode('utf-8')
            # parse the response stirng with xml.sax
            try:
                xml.sax.parseString(result_string, handler)
            except xml.sax.SAXParseException:
                print('Error: SAXParseException for {}'.format(accession))
                continue
            # store the metadata for the biosample id in results_dict
            # filter the content_dict to only include the values in values_file if provided
            if set_values:
                filtered_dict = {k: content_dict[k] for k in content_dict if k in values}
                results_dict[accession] = filtered_dict
            else:
                results_dict[accession] = content_dict

    # add headers to the csv_headers if there is no values file
    if csv_headers == ['biosample_id']:
        for _, data_dict in results_dict.items():
            for header, _ in data_dict.items():
                if header not in csv_headers:
                    csv_headers.append(header)

    # write results to csv file
    with open(output_file, 'w', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(csv_headers)
        for accession, data_dict in results_dict.items():
            row = [accession]
            for header in csv_headers[1:]:
                if header in data_dict:
                    row.append(data_dict[header])
                else:
                    row.append('na')
            writer.writerow(row)

########################################################################################################################
################################################ biosample section end #################################################
########################################################################################################################


def sep_path_basename_ext(file_in):
    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]
    return f_name, f_path, f_base, f_ext


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):

    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


def pebblescout_results_to_data_matrix(file_base_list, db_list, pebblescout_op_dir, min_cov_pct, op_txt, op_txt_t, op_txt_t_with_metadata, metadata_col_list, metadata_dict):

    subject_to_biosample_dict = dict()
    biosample_title_dict = dict()
    subject_id_set_all = set()
    biosample_subject_id_set_all = set()
    pebblescout_op_dod = dict()
    for f_base in file_base_list:
        current_gnm_dict = dict()
        for each_db in db_list:
            op_file = '%s/%s_%s.txt' % (pebblescout_op_dir, f_base, each_db)
            if os.path.isfile(op_file):
                for each_line in open(op_file):
                    if not each_line.startswith('QueryID\tSubjectID'):
                        each_line_split = each_line.strip().split('\t')
                        subject_id      = each_line_split[1]
                        pct_coverage    = each_line_split[3]
                        bioSample_id    = each_line_split[5]
                        bioSample_title = each_line_split[6]
                        if (bioSample_id != '""') and (bioSample_title != '""'):
                            if float(pct_coverage) >= min_cov_pct:
                                current_gnm_dict[subject_id] = pct_coverage
                                subject_to_biosample_dict[subject_id] = bioSample_id
                                biosample_title_dict[bioSample_id] = bioSample_title
                                subject_id_set_all.add(subject_id)
                                biosample_subject_id_set_all.add('%s_%s' % (bioSample_id, subject_id))
        pebblescout_op_dod[f_base] = current_gnm_dict

    if len(subject_id_set_all) > 0:

        subject_id_set_all_with_biosample_sorted = sorted(list(biosample_subject_id_set_all))
        subject_desc_list_sorted = []
        for each_biosample_subject in subject_id_set_all_with_biosample_sorted:
            biosample_id    = each_biosample_subject.split('_')[0]
            biosample_title = biosample_title_dict.get(biosample_id, 'NA')
            subject_desc_list_sorted.append(biosample_title)

        # write out to file
        op_txt_handle = open(op_txt, 'w')
        op_txt_handle.write('Genome\t%s\n' % ('\t'.join(subject_id_set_all_with_biosample_sorted)))
        for each_gnm in sorted(list(pebblescout_op_dod.keys())):
            current_gnm_dict = pebblescout_op_dod[each_gnm]
            value_list = [each_gnm]
            for each_biosample_subject in subject_id_set_all_with_biosample_sorted:
                subject_id = each_biosample_subject.split('_')[1]
                value_list.append(current_gnm_dict.get(subject_id, '0'))
            op_txt_handle.write('%s\n' % ('\t'.join(value_list)))
        op_txt_handle.write('Description\t%s\n' % ('\t'.join(subject_desc_list_sorted)))
        op_txt_handle.close()

        transpose_csv(op_txt, op_txt_t, '\t', 0, 0)
        os.system('rm -r %s' % op_txt)

        biosample_to_keep = set()
        n = 0
        for each_line in open(op_txt_t, 'r'):
            if n > 0:
                biosample_id = each_line.strip().split('\t')[0].split('_')[0]
                biosample_to_keep.add(biosample_id)
            n += 1

        col_index_to_value_dict = dict()
        for each_sample in metadata_dict:
            if each_sample in biosample_to_keep:
                index = 0
                for each_element in metadata_dict[each_sample]:
                    if index not in col_index_to_value_dict:
                        col_index_to_value_dict[index] = set()
                    col_index_to_value_dict[index].add(each_element)
                    index += 1

        cols_to_ignore_set = set()
        for each_col in col_index_to_value_dict:
            if (len(col_index_to_value_dict[each_col]) == 1) and ('na' in col_index_to_value_dict[each_col]):
                cols_to_ignore_set.add(each_col)

        if len(metadata_dict) > 0:
            op_txt_t_with_metadata_handle = open(op_txt_t_with_metadata, 'w')
            n = 0
            for each_line in open(op_txt_t, 'r'):
                if n == 0:
                    list_to_use  = metadata_col_list
                else:
                    biosample_id = each_line.strip().split('\t')[0].split('_')[0]
                    list_to_use  = metadata_dict[biosample_id]

                list_to_use_removed_na_col = []
                col_index = 0
                for each_col in list_to_use:
                    if col_index not in cols_to_ignore_set:
                        list_to_use_removed_na_col.append(each_col)
                    col_index += 1

                if n == 0:
                    op_txt_t_with_metadata_handle.write('Dataset\t%s\t%s\n' % (each_line.strip(), '\t'.join(list_to_use_removed_na_col)))
                else:
                    op_txt_t_with_metadata_handle.write('%s\t%s\n' % (each_line.strip(), '\t'.join(list_to_use_removed_na_col)))
                n += 1
            op_txt_t_with_metadata_handle.close()

            os.system('rm %s' % op_txt_t)


def PebbleScout(args):

    file_dir            = args['i']
    file_ext            = args['x']
    db_str              = args['db']
    # min_cov_pct         = args['min_cov']
    num_threads         = args['t']
    force_create_op_dir = args['f']
    op_dir              = args['o']
    gnm_id_txt          = args['id']
    op_prefix           = args['p']
    get_metadata        = args['m']

    if op_prefix is None:
        prefix_str = ''
    else:
        prefix_str = '%s_' % op_prefix

    db_list             = db_str.split(',')
    cov_pct_list        = [0, 1, 3, 5, 10, 20, 30, 50]

    # define file name
    gnm_dir                 = '%s/concatenated_genome'              % op_dir
    pebblescout_op_dir      = '%s/pebblescout_output'               % op_dir
    cmd_txt                 = '%s/commands.txt'                     % op_dir
    biosample_txt           = '%s/BioSample.txt'                    % op_dir
    biosample_metadata_txt  = '%s/BioSample_metadata.txt'           % op_dir

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % gnm_dir)
    os.system('mkdir %s' % pebblescout_op_dir)

    gnm_id_set = set()
    if gnm_id_txt is not None:
        if os.path.isfile(gnm_id_txt) is True:
            for each_gnm in open(gnm_id_txt):
                gnm_id_set.add(each_gnm.strip().split()[0])
        else:
            print('%s not found, program exited!' % gnm_id_txt)
            exit()

    # concatenate sequences in input genomes
    file_re = '%s/*.%s' % (file_dir, file_ext)
    file_list = glob.glob(file_re)
    file_base_list_to_process = []
    for each_gnm in file_list:
        f_name, f_path, f_base, f_ext = sep_path_basename_ext(each_gnm)

        to_process = ''
        if len(gnm_id_set) == 0:
            to_process = True
        else:
            if (f_base in gnm_id_set) or (f_name in gnm_id_set):
                to_process = True

        if to_process is True:
            file_base_list_to_process.append(f_base)
            current_gnm_after_cat = '%s/%s.%s' % (gnm_dir, f_base, f_ext)
            concatenates_seq_str = ''
            seq_index = 0
            for each_seq in SeqIO.parse(each_gnm, 'fasta'):
                if seq_index == 0:
                    concatenates_seq_str = concatenates_seq_str + str(each_seq.seq).strip()
                else:
                    concatenates_seq_str = concatenates_seq_str + ('NNNNNNNNNNNNNNNNNNNNNNNNN%s' % str(each_seq.seq).strip())
                seq_index += 1

            current_gnm_after_cat_handle = open(current_gnm_after_cat, 'w')
            current_gnm_after_cat_handle.write('>%s\n' % f_base)
            current_gnm_after_cat_handle.write(concatenates_seq_str + '\n')
            current_gnm_after_cat_handle.close()

    # get pebblescout commands
    cmd_txt_handle = open(cmd_txt, 'w')
    cmd_set = set()
    pebblescout_op_txt_list = []
    for f_base in file_base_list_to_process:
        for each_db in db_list:
            seq_file           = '%s/%s.%s'                                                                                                % (gnm_dir, f_base, file_ext)
            pebblescout_op_txt = '%s/%s_%s.txt'                                                                                            % (pebblescout_op_dir, f_base, each_db)
            url_str            = 'https://pebblescout.ncbi.nlm.nih.gov/sra-cl-be/sra-cl-be.cgi?db=%s&m=2&rettype=pebblescout&download=yes' % each_db
            pebblescout_cmd    = 'curl -s -F "fasta=@%s" "%s" > %s'                                                                        % (seq_file, url_str, pebblescout_op_txt)
            cmd_set.add(pebblescout_cmd)
            pebblescout_op_txt_list.append(pebblescout_op_txt)
            cmd_txt_handle.write(pebblescout_cmd + '\n')
    cmd_txt_handle.close()

    # run pebblescout with multi-processing
    print('Running %s commands with %s cores' % (len(cmd_set), num_threads))
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, sorted([i for i in cmd_set]))
    pool.close()
    pool.join()

    # get biosample set
    biosample_set = set()
    for each_pebblescout_op_txt in pebblescout_op_txt_list:
        for each_line in open(each_pebblescout_op_txt):
            if not each_line.startswith('QueryID\tSubjectID'):
                each_line_split = each_line.strip().split('\t')
                bioSample_id = each_line_split[5]
                if bioSample_id != '""':
                    biosample_set.add(bioSample_id)

    # get biosample metadata
    biosample_metadata_col_list = []
    biosample_metadata_dict = dict()
    if get_metadata is True:
        with open(biosample_txt, 'w') as f:
            f.write('\n'.join(sorted(list(biosample_set))))
        biosample({'i': biosample_txt, 'o': biosample_metadata_txt})
        os.system('rm %s' % biosample_txt)

        n = 0
        col_num = 0
        for each_line in open(biosample_metadata_txt):
            each_line_split     = each_line.strip().split('\t')
            if n == 0:
                biosample_metadata_col_list = each_line_split[1:]
                col_num = len(each_line_split)
            else:
                biosample_id = each_line_split[0]
                biosample_meta_list = each_line_split[1:]
                if len(each_line_split) == col_num:
                    biosample_metadata_dict[biosample_id] = biosample_meta_list
                else:
                    biosample_metadata_dict[biosample_id] = ['na']*(col_num - 1)
            n +=1

    # get data matrix
    for each_cov_pct in cov_pct_list:
        df_txt_tmp      = '%s/%sPebbleScout_all_tmp.txt'                % (op_dir, prefix_str)
        df_txt          = '%s/%sPebbleScout_all.txt'                    % (op_dir, prefix_str)
        df_txt_meta     = '%s/%sPebbleScout_all_with_metadata.txt'      % (op_dir, prefix_str)
        if each_cov_pct > 0:
            df_txt_tmp  = '%s/%sPebbleScout_cov%s_tmp.txt'              % (op_dir, prefix_str, each_cov_pct)
            df_txt      = '%s/%sPebbleScout_cov%s.txt'                  % (op_dir, prefix_str, each_cov_pct)
            df_txt_meta = '%s/%sPebbleScout_cov%s_with_metadata.txt'    % (op_dir, prefix_str, each_cov_pct)
        pebblescout_results_to_data_matrix(file_base_list_to_process, db_list, pebblescout_op_dir, each_cov_pct, df_txt_tmp, df_txt, df_txt_meta, biosample_metadata_col_list, biosample_metadata_dict)

    print('Done!')


if __name__ == '__main__':

    PebbleScout_parser = argparse.ArgumentParser(usage=PebbleScout_usage)
    PebbleScout_parser.add_argument('-p',           required=False, default=None,                   help='output prefix')
    PebbleScout_parser.add_argument('-i',           required=True,                                  help='input fasta file/dir')
    PebbleScout_parser.add_argument('-x',           required=False, default=None,                   help='file extension')
    PebbleScout_parser.add_argument('-id',          required=False, default=None,                   help='id of genomes to process')
    PebbleScout_parser.add_argument('-db',          required=False, default='meta,meta_vol2',       help='query database, default is meta and meta_vol2')
    PebbleScout_parser.add_argument('-m',           required=False, action="store_true",            help='retrieve biosample metadata from NCBI')
    # PebbleScout_parser.add_argument('-min_cov',     required=False, type=str, default='1,3,5,10',   help='minimum %coverage to include in the datamatrix, default is 1,3,5,10')
    PebbleScout_parser.add_argument('-t',           required=False, type=int, default=1,            help='number of core, default is 1')
    PebbleScout_parser.add_argument('-f',           required=False, action="store_true",            help='force overwrite')
    PebbleScout_parser.add_argument('-o',           required=True,                                  help='output directory')
    args = vars(PebbleScout_parser.parse_args())
    PebbleScout(args)
