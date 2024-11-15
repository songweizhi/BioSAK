
fastq_r1 = 'R1.fastq.gz'
file_ext = 'fastq.gz'

r1_base = fastq_r1[:-(len(file_ext)+1)]

print(r1_base)
