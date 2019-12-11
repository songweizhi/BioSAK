import sys
import os
from Bio import SeqIO

records = SeqIO.parse('/Users/songweizhi/Desktop/sequence.gbk', 'genbank')
for record in records:
    for feature in record.features:
        if 'locus_tag' in feature.qualifiers:
            print(record)
            if feature.qualifiers['locus_tag'] == 'PTD2_03541':
                print(record)
            # if feature.locus_tag == 'PTD2_03541':
            #     print(feature)