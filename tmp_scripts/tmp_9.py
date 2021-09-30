
for each in open('/Users/songweizhi/Desktop/bin_list_with_path.txt'):

    pwd_bin = each.strip()
    bin_id = pwd_bin.split('/')[-1].split('.')[0]

    blastn_cmd = 'blastn -query %s -db combined_hc_refs.fa -out %s_vs_ref.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 12' % (pwd_bin, bin_id)
    print(blastn_cmd)