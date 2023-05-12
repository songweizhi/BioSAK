BioSAK MeanMappingDepth -depth ctgs.bam.depth
BioSAK MeanMappingDepth -depth ctgs.bam.depth -T

# get depth file with samtools
samtools depth ctgs.bam > ctgs.bam.depth
