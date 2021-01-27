from Bio.Blast import NCBIXML


def blast_xml_parser(blast_result_xml, blast_result_with_description):

    # columns in output file
    # query	subject	identity	length	evalue	bitscore	subject_description

    blast_result_with_description_handle = open(blast_result_with_description, 'w')
    blast_result_with_description_handle.write('query\tsubject\tidentity\tlength\tevalue\tbitscore\tsubject_description\n')
    for blast_record in NCBIXML.parse(open(blast_result_xml)):
        query_id = blast_record.query.split(' ')[0]
        for alignment in blast_record.alignments:
            subject_id = alignment.hit_id
            subject_def = ' '.join(alignment.hit_def.split(' ')[1:])
            for hsp in alignment.hsps:
                identity = float("{0:.3f}".format(hsp.identities*100/hsp.align_length))
                aln_len = hsp.align_length
                e_value = hsp.expect
                bit_score = hsp.bits
                blast_result_with_description_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (query_id, subject_id, identity, aln_len, e_value, bit_score, subject_def))
    blast_result_with_description_handle.close()


blast_result_xml         = '/Users/songweizhi/Desktop/query_blastp_xml.tab'
blast_result_xml_summary = '/Users/songweizhi/Desktop/query_blastp_xml_summary.tab'
blast_xml_parser(blast_result_xml, blast_result_xml_summary)


