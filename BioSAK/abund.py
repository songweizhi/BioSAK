import os
import glob
import argparse
from Bio import SeqIO


abund_usage = '''
========================== mask example commands ==========================

BioSAK mask -h
BioSAK abund -h

cd /Users/songweizhi/Desktop
python3 /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/abund.py -i symbiont_genomes -x fna -o symbiont_genomes_abund -f 

python3 abund.py -i symbiont_genomes -x fna -o symbiont_genomes_abund -f -t 6

===========================================================================
'''

def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def mask(args):

    gnm_dir         = args['i']
    file_ext        = args['x']

def mapping(fq_r1, fq_r2, ref_seq, op_prefix, num_threads):

    bwa_cmd             = 'bwa mem -5SP -t %s %s %s %s | samblaster > %s.sam'                                                           % (num_threads, ref_seq, fq_r1, fq_r2, op_prefix)
    samtools_view_cmd   = 'samtools view -@ 32 -bS -h -b %s.sam > %s.bam'                                                               % (op_prefix, op_prefix)
    samtools_sort_cmd   = 'samtools sort -@ 32 %s.bam -o %s.sorted.bam'                                                                 % (op_prefix, op_prefix)
    bamm_filter_cmd     = 'bamm filter -b %s.sorted.bam --percentage_id 0.99 --percentage_aln 0.9'                                      % (op_prefix)
    coverm_filter_cmd   = 'coverm filter -b %s.sorted.bam --min-read-aligned-percent 0.9 --min-read-percent-identity 0.99'              % (op_prefix)
    pileup_sh_cmd       = 'pileup.sh in=%s.sorted_filtered.bam out=%s.sorted_filtered.cov rpkm=%s.sorted_filtered.rpkm overwrite=true'  % (op_prefix, op_prefix, op_prefix)
    seqkit_stat_cmd     = 'seqkit stat %s.fastq > %s.stat'                                                                              % (op_prefix, op_prefix)

    os.system(bwa_cmd)
    os.system(samtools_view_cmd)
    os.system(samtools_sort_cmd)
    os.system(coverm_filter_cmd)
    os.system(pileup_sh_cmd)
    os.system(seqkit_stat_cmd)


def get_rpkm():
    pass


def abund(args):

    gnm_dir             = args['i']
    file_ext            = args['x']
    output_folder       = args['o']
    num_threads         = args['t']
    force_create_op_dir = args['f']

    gnm_file_re = '%s/*.%s' % (gnm_dir, file_ext)
    input_gnm_list = glob.glob(gnm_file_re)

    # create output folder
    if os.path.isdir(output_folder) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % output_folder)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % output_folder)

    concatenated_renamed_gnm = '%s/combined.fasta' % output_folder

    concatenated_renamed_fasta_handle = open(concatenated_renamed_gnm, 'w')
    gnm_name_dict_old_to_new = dict()
    gnm_name_dict_new_to_old = dict()
    for each_gnm in input_gnm_list:
        gnm_name, gnm_path, gnm_base, gnm_ext = sep_path_basename_ext(each_gnm)
        gnm_name_new = gnm_base
        gnm_name_new = gnm_name_new.replace('_', '').replace('-', '')
        gnm_name_dict_old_to_new[gnm_base] = gnm_name_new
        gnm_name_dict_new_to_old[gnm_name_new] = gnm_base
        seq_index = 1
        for each_seq in SeqIO.parse(each_gnm, 'fasta'):
            concatenated_renamed_fasta_handle.write('>%s_%s\n' % (gnm_name_new, seq_index))
            concatenated_renamed_fasta_handle.write('%s\n' % each_seq.seq)
            seq_index += 1
    concatenated_renamed_fasta_handle.close()

    metaxa2_ssu_cmd = 'metaxa2 --plus --mode m --cpu 12 --multi_thread T --table T -g ssu --not_found T -i %s -o %s.metaxa2_ssu' % (concatenated_renamed_gnm, concatenated_renamed_gnm)
    metaxa2_lsu_cmd = 'metaxa2 --plus --mode m --cpu 12 --multi_thread T --table T -g lsu --not_found T -i %s -o %s.metaxa2_lsu' % (concatenated_renamed_gnm, concatenated_renamed_gnm)
    os.system(metaxa2_ssu_cmd)
    os.system(metaxa2_lsu_cmd)




if __name__ == '__main__':

    abund_parser = argparse.ArgumentParser(usage=abund_usage)
    abund_parser.add_argument('-i', required=True,                          help='input genome directory')
    abund_parser.add_argument('-x', required=True,                          help='file extension')
    abund_parser.add_argument('-o', required=True,                          help='output directory')
    abund_parser.add_argument('-f', required=False, action="store_true",    help='force overwrite')
    abund_parser.add_argument('-t', required=False, type=int, default=1,    help='number of threads')
    args = vars(abund_parser.parse_args())
    abund(args)


'''
cd /scratch/PI/ocessongwz/Sponge_r220/7_symbiont_abundance/symbiont_genomes
metaxa2 --plus --mode m --cpu 12 --multi_thread T --table T -g ssu --not_found T -i VYC.fna -o VYC.metaxa2_ssu
metaxa2 --plus --mode m --cpu 12 --multi_thread T --table T -g lsu --not_found T -i VYC.fna -o VYC.metaxa2_lsu
cut -f 1,9,10 VYC.metaxa2_ssu.extraction.results > masked_metaxa_ssu.bed
cut -f 1,9,10 VYC.metaxa2_lsu.extraction.results > masked_metaxa_lsu.bed
cat masked_metaxa_ssu.bed masked_metaxa_lsu.bed > masked_metaxa.bed
barrnap --kingdom bac --threads 32 --reject 0.3 VYC.fna > out_bac.gff
barrnap --kingdom arc --threads 32 --reject 0.3 VYC.fna > out_arc.gff
barrnap --kingdom euk --threads 32 --reject 0.3 VYC.fna > out_euk.gff
cut -f 1,4,5 out_bac.gff >> masked_barrnap.bed
cut -f 1,4,5 out_arc.gff >> masked_barrnap.bed
cut -f 1,4,5 out_euk.gff >> masked_barrnap.bed
cat masked_barrnap.bed masked_metaxa.bed > masked_rrna.bed

# run tRNAscan-SE
tRNAscan-SE -A -b VYC.fna.ar.bedformat --thread 1 VYC.fna
tRNAscan-SE -B -b VYC.fna.bac.bedformat --thread 1 VYC.fna
cat VYC.fna.ar.bedformat VYC.fna.bac.bedformat > VYC.fna.bedformat
cut -f 1,2,3 VYC.fna.bedformat > masked_trnascan.bed

# Find duplicate areas
dustmasker -in VYC.fna -out VYC.fna.lowcom.out -outfmt acclist
cut -d '>' -f 2 VYC.fna.lowcom.out > masked_dustmasker.bed 

# Cover the above area with NNNN
cat masked_rrna.bed masked_trnascan.bed masked_dustmasker.bed > VYC.fna_mask.bed
awk -F'\t' '$2>=0' VYC.fna_mask.bed > VYC.fna_mask_0.bed
bedtools maskfasta -fi VYC.fna -bed VYC.fna_mask_0.bed -fo VYC.fna_masked.fasta
'''
