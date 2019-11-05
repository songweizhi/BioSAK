import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

usage = '''
python3 /Users/songweizhi/PycharmProjects/python_scripts/DNA2protein.py -in input_nc.fna -out output_aa.faa
'''

parser = argparse.ArgumentParser()

parser.add_argument('-in',
                    required=True,
                    help='input nucleotide sequences file')

parser.add_argument('-out',
                    required=False,
                    help='output amino acid sequences file')

args = vars(parser.parse_args())
input_file = args['in']
output_file = args['out']

output_handle = open(output_file, 'w')
for nc_record in SeqIO.parse(input_file, 'fasta'):
    aa_object = nc_record.seq.translate()
    aa_record = SeqRecord(aa_object)
    aa_record.id = nc_record.id
    aa_record.description = nc_record.description
    SeqIO.write(aa_record, output_handle, 'fasta')
output_handle.close()
