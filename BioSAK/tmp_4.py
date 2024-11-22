#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import io, sys

fn1 = sys.argv[1] #GSD-20.sorted_filtered.rpkm
fn2 = sys.argv[2] #allsample.stat
fn3 = sys.argv[3] #test_out.txt

def get_sample_reads(sample_name):
    with io.open (fn2,'r', encoding = 'utf-8' ) as f:
        lines = f.readlines()
        for line in lines:
            if str(sample_name) == line.split('\t')[0].split('_1.fastq')[0]: #sample name GSD-20_1.fastq
                sample_reads = line.split('\t')[1]
                sample_reads_fnl = int(sample_reads.replace('\n','').replace(',', ''))
                return  sample_reads_fnl

if __name__ == '__main__':

    rpkm_out1 = pd.read_csv(fn1, sep = '\t',header = 4)

    rpkm_out2 = rpkm_out1[['#Name', 'Length', 'Reads']]
    rpkm_out2['bin_name'] = rpkm_out2['#Name'].str.rsplit('_', n=1, expand=True)[0] #contig name 'bin1_1', bin name 'bin1'
    rpkm_out2['Length'] =  rpkm_out2['Length'].astype(int)
    rpkm_out2['Reads'] =  rpkm_out2['Reads'].astype(int)

    rpkm_out3 = rpkm_out2.groupby(["bin_name"])[['Length','Reads']].sum()

    sample_name = str(fn1).split('.')[0]
    sample_reads_fnl = get_sample_reads(sample_name)
    rpkm_out3[sample_name+'_rpkm']=(rpkm_out3['Reads']*1000000000)/(rpkm_out3['Length']*sample_reads_fnl*2) #calulate RPKM
    rpkm_out3.to_csv(fn3,sep='\t')


