#!/usr/bin/env python
"""
A utility script to split a large FASTA file into smaller ones, with an 
arbitrary number of sequences in each one. 

Takes two arguments:
    1) The large FASTA file that needs to be split
    2) The number of sequences that should be put into each FASTA
"""
#   define a usage message
usage = """Split_FASTA.py - splits a large FASTA file into smaller ones.

Usage:
Split_FASTA.py [FILE] [NUM_SEQS]

Splits [FILE] into a collection of smaller files, with [NUM_SEQS] sequences
in each one."""

#   To take arguments on the command line
import sys
#   To do some math
import math
#   To read/write FASTA. Wrap in a try/except
try:
    from Bio import SeqIO
except ImportError:
    print "This script requires BioPython to be installed!"

#   Simple check to make sure that the correct number of arguments are given
if len(sys.argv) != 3:
    print usage
    exit(1)

#   Set these names so that it is easier to follow them in the code
to_split = sys.argv[1]
num_seqs = sys.argv[2]

#   Make sure that the file provided is readable. 
try:
    open(to_split, 'r')
except IOError:
    print "Error! " + to_split + " does not exist, or is not readable!"
    exit(1)

#   And make sure that the second argument is a positive integer
#   The try-Except block tests for integer, and the if tests for positive
try:
    int(num_seqs)
except ValueError:
    print "Error! Please provide a positive integer!"
    exit(1)
if int(num_seqs) <= 0:
    print "Error! Please provide a positive integer!"
    exit(1)

#   Read in the FASTA file
all_fasta = list(SeqIO.parse(to_split, 'fasta'))
#   How many files to split into?
num_files = int(math.ceil(len(all_fasta)/float(num_seqs)))
#   Print a little message
print "Will split " + to_split + " into " + str(num_files) + " files, with " + str(num_seqs) + " seqs per file."
#   Start splitting!
i = 1
while i <= num_files:
    #   Calculate the start
    start = int(i-1) * int(num_seqs)
    #   and the end
    end = int(i) * int(num_seqs)
    #   Generate the filename
    #   strip off the extension, and get the rest
    filename = to_split.split('.')[:-1][0] + '_' + str(i) + '.fasta'
    #   Write the sequences into the file
    SeqIO.write(all_fasta[start:end], filename, 'fasta')
    #   increment the counter
    i += 1

