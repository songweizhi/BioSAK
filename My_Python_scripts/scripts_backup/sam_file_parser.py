
import pysam

samfile = pysam.AlignmentFile("/Users/songweizhi/Desktop/sam_file_parser/CC6-YY-74_recircularise_sorted.bam", "rb")

iter = samfile.fetch("Chr1", 10, 20)
n = 0
for x in iter:
    print(str(x))
    n += 1
    print()

print(n)




