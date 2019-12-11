from Bio import SeqIO

input_fastq = '/Users/songweizhi/Desktop/test.fastq'
output_fastq_R1 = '/Users/songweizhi/Desktop/test_R1.fastq'
output_fastq_R2 = '/Users/songweizhi/Desktop/test_R2.fastq'

input_fastq = SeqIO.parse(input_fastq, 'fastq')

R1_handle = open(output_fastq_R1, 'w')
R2_handle = open(output_fastq_R2, 'w')

for read in input_fastq:
    description = read.description
    description_split = description.strip().split(' ')
    if read.id == 'SOLEXA4:34:C13W1ACXX:7:1101:1926:2143':
        pass
    if description_split[1][0] == '1':
        SeqIO.write(read, R1_handle, 'fastq')
    if description_split[1][0] == '2':
        SeqIO.write(read, R2_handle, 'fastq')

R1_handle.close()
R2_handle.close()
