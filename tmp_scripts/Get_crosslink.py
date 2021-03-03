import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-q',
                    required=True,
                    help='query sequences (multi-fasta format)')

parser.add_argument('-r',
                    required=True,
                    help='reference sequences folder')

parser.add_argument('-x',
                    required=True,
                    help='extension of reference file(s)')

parser.add_argument('-i',
                    required=True,
                    type=float,
                    help='identity cutoff')

parser.add_argument('-a',
                    required=True,
                    type=int,
                    help='alignment length cutoff')

parser.add_argument('-makeblastdb',
                    required=False,
                    default='makeblastdb',
                    help='makeblastdb')

parser.add_argument('-blastn',
                    required=False,
                    default='blastn',
                    help='blastn')


args = vars(parser.parse_args())

query_contigs = args['q']
references = args['r']
reference_extension = args['x']
identity_cutoff = args['i']
alignment_length_cutoff = args['a']
pwd_makeblastdb_exe = args['makeblastdb']
pwd_blastn_exe = args['blastn']

if references[-1] == '/':
    input_bin_folder_1 = references[:-1]

wd = os.getcwd()


reference_folder = '%s/%s' % (wd, references)
os .system('cat references/*.%s > combined_references.fasta' % (reference_extension))

pwd_combined_references = '%s/combined_references.fasta' % wd
pwd_query_contigs = '%s/%s' % (wd, query_contigs)
pwd_blast_result = '%s/blastn_results.tab' % wd
pwd_blast_result_filtered = '%s/blastn_results_filtered.tab' % wd
pwd_blast_result_filtered_sorted = '%s/blastn_results_filtered_sorted.tab' % wd
pwd_googleVis_input = '%s/googleVis_input.csv' % wd
pwd_googleVis_output = '%s/googleVis_output.html' % wd

# make blast database
print('Making Blastdb with combined_references.fasta')
os.system('%s -in %s -dbtype nucl -parse_seqids' % (pwd_makeblastdb_exe, pwd_combined_references))

# run blastn
outfmt = '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"'
blast_parameters = '-evalue 1e-5 -perc_identity %s -outfmt 6 -task blastn' % identity_cutoff
print('\nRunning Blast, be patient...')
os.system('%s -query %s -db %s -out %s -perc_identity 100 %s' % (pwd_blastn_exe, pwd_query_contigs, pwd_combined_references, pwd_blast_result, outfmt))
os.system('%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, pwd_query_contigs, pwd_combined_references, pwd_blast_result, blast_parameters))

# process blast results
pwd_blast_result_filtered_handle = open(pwd_blast_result_filtered, 'w')
for each in open(pwd_blast_result):
    #print(each.strip())
    each_split = each.strip().split('\t')
    query = each_split[0]
    subject = each_split[1]
    identity = float(each_split[2])
    alignment_length = int(each_split[3])
    if (identity >= identity_cutoff) and (alignment_length >= alignment_length_cutoff):
        if '_ctg' in subject:
            subject_split = subject.split('_')
            subject_new = subject_split[0]
            pwd_blast_result_filtered_handle.write('%s\t%s\t%s\n' % (subject_new, query, alignment_length))
        else:
            pwd_blast_result_filtered_handle.write('%s\t%s\t%s\n' % (subject, query, alignment_length))
pwd_blast_result_filtered_handle.close()

os.system('cat %s | sort > %s' % (pwd_blast_result_filtered, pwd_blast_result_filtered_sorted))

# get googleVis input
current_match = ''
current_match_lengths = []
pwd_googleVis_input_handle = open(pwd_googleVis_input, 'w')
pwd_googleVis_input_handle.write('C1,C2,Length\n')

for each_match in open(pwd_blast_result_filtered_sorted):
    each_match_split = each_match.strip().split('\t')
    each_match_identifier = '%s,%s' % (each_match_split[0], each_match_split[1])
    each_match_alignment_length = int(each_match_split[2])
    if current_match == '':
        current_match = each_match_identifier
        current_match_lengths.append(each_match_alignment_length)
    elif current_match == each_match_identifier:
        current_match_lengths.append(each_match_alignment_length)
    elif current_match != each_match_identifier:
        pwd_googleVis_input_handle.write('%s,%s\n' % (current_match, sum(current_match_lengths)))
        current_match = each_match_identifier
        current_match_lengths = []
        current_match_lengths.append(each_match_alignment_length)


pwd_googleVis_input_handle.write('%s,%s\n' % (current_match, sum(current_match_lengths)))
pwd_googleVis_input_handle.close()

import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages

def GoogleVis_Sankey_plotter(input_csv, output_html, height):
    out = open(output_html, 'w')
    utils = rpackages.importr('googleVis')
    packages_needed = ['googleVis']
    for each_package in packages_needed :
        if not rpackages.isinstalled(each_package):
            utils.install_packages(each_package)
        else:
            pass

    df = robjects.DataFrame.from_csvfile(input_csv)
    sankey_plot = robjects.r['gvisSankey'](df,
                                           option=robjects.r['list'](
                                               sankey="{node : {colorMode: 'unique', labelPadding: 10}, link:{colorMode: 'source'}}",
                                               height=height,
                                               width=600))
    out.write(str(sankey_plot))
    out.close()


GoogleVis_Sankey_plotter(pwd_googleVis_input, pwd_googleVis_output, 900)
