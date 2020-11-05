import argparse
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment


SubsetAlnCols_usage = '''
========================= SubsetAlnCols example commands =========================

BioSAK SubsetAlnCols -in 16S.aln -r 200-300 -out 16S_v1.aln
BioSAK SubsetAlnCols -in 16S.aln -r 200-300 -pct 80 -out 16S_v1.aln

===============================================================================
'''


def SubsetAlnCols(args):

    aln_in_fasta      = args['in']
    col_range         = args['r']
    min_non_empty_pct = args['pct']
    align_subset_out  = args['out']


    # get range
    col_start        = int(col_range.split('-')[0])
    col_end          = int(col_range.split('-')[1])
    align_subset_len = col_end - col_start + 1

    # subset
    align_in         = AlignIO.read(aln_in_fasta, 'fasta')
    align_in_subset  = align_in[:, (col_start - 1):col_end]

    if min_non_empty_pct is None:

        # write out
        align_file_out_handle = open(align_subset_out, 'w')
        AlignIO.write(align_in_subset, align_file_out_handle, 'fasta')
        align_file_out_handle.close()

    else:

        # filter
        align_in_subset_filtered = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
        for seq_record in align_in_subset:
            non_empty_col_pct = (align_subset_len - seq_record.seq.count('-'))*100/align_subset_len
            if non_empty_col_pct >= float(min_non_empty_pct):
                align_in_subset_filtered.add_sequence(seq_record.id, str(seq_record.seq))

        # write out
        align_file_out_handle = open(align_subset_out, 'w')
        AlignIO.write(align_in_subset_filtered, align_file_out_handle, 'fasta')
        align_file_out_handle.close()


if __name__ == '__main__':

    SubsetAlnCols_parser = argparse.ArgumentParser()

    # arguments for rename_seq_parser
    SubsetAlnCols_parser.add_argument('-in',    required=True,                help='input MSA in fasta format')
    SubsetAlnCols_parser.add_argument('-r',     required=True,                help='columns to keep, e.g. 200-300, one based')
    SubsetAlnCols_parser.add_argument('-pct',   required=False, default=None, help='minimum percentage of nonempty bases (e.g. 70), default keep all')
    SubsetAlnCols_parser.add_argument('-out',   required=True,                help='output file')

    args = vars(SubsetAlnCols_parser.parse_args())

    SubsetAlnCols(args)

