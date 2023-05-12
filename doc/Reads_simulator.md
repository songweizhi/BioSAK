# simulate 100000 pairs of reads from ref.fa
BioSAK Reads_simulator -p Test -r ref.fa -l 250 -i 300 -split -n 100000

# simulate reads to a depth of 50X
BioSAK Reads_simulator -p Test -r ref.fa -l 250 -i 300 -split -d 50
