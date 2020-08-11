#!/usr/bin/env python

import sys

__author__ = 'weizhisong'

usage = """
####################################################################

  Usage:
  python COG_number_count.py annotation_results.txt output.txt

  COG_annotation_results format:
  G_00005	yebN	S	COG1971	Predicted membrane protein
  G_00009	PA1596	O	COG0326	Molecular chaperone, HSP90 family
  ...

  Author: Weizhi Song

####################################################################
"""

if len(sys.argv) != 3:
    print (usage)
    exit(1)

# Generate COG list for all genes:

input_file = sys.argv[1]
output_file = sys.argv[2]

annotation_results = open(input_file)
out = open(output_file, 'w')

all_cogs = []
for each_gene in annotation_results:
    cog_category = each_gene.strip().split('\t')[2]
    all_cogs.append(cog_category)


all_cogs_temp = []
for item in all_cogs:
    if len(item) == 1:
        all_cogs_temp.append(item)
    else:
        list(item)
        all_cogs_temp += list(item)

Count_A = 'A' + '\t' + str(all_cogs_temp.count('A')) + '\n'
Count_B = 'B' + '\t' + str(all_cogs_temp.count('B')) + '\n'
Count_C = 'C' + '\t' + str(all_cogs_temp.count('C')) + '\n'
Count_D = 'D' + '\t' + str(all_cogs_temp.count('D')) + '\n'
Count_E = 'E' + '\t' + str(all_cogs_temp.count('E')) + '\n'
Count_F = 'F' + '\t' + str(all_cogs_temp.count('F')) + '\n'
Count_G = 'G' + '\t' + str(all_cogs_temp.count('G')) + '\n'
Count_H = 'H' + '\t' + str(all_cogs_temp.count('H')) + '\n'
Count_I = 'I' + '\t' + str(all_cogs_temp.count('I')) + '\n'
Count_J = 'J' + '\t' + str(all_cogs_temp.count('J')) + '\n'
Count_K = 'K' + '\t' + str(all_cogs_temp.count('K')) + '\n'
Count_L = 'L' + '\t' + str(all_cogs_temp.count('L')) + '\n'
Count_M = 'M' + '\t' + str(all_cogs_temp.count('M')) + '\n'
Count_N = 'N' + '\t' + str(all_cogs_temp.count('N')) + '\n'
Count_O = 'O' + '\t' + str(all_cogs_temp.count('O')) + '\n'
Count_P = 'P' + '\t' + str(all_cogs_temp.count('P')) + '\n'
Count_Q = 'Q' + '\t' + str(all_cogs_temp.count('Q')) + '\n'
Count_R = 'R' + '\t' + str(all_cogs_temp.count('R')) + '\n'
Count_S = 'S' + '\t' + str(all_cogs_temp.count('S')) + '\n'
Count_T = 'T' + '\t' + str(all_cogs_temp.count('T')) + '\n'
Count_U = 'U' + '\t' + str(all_cogs_temp.count('U')) + '\n'
Count_V = 'V' + '\t' + str(all_cogs_temp.count('V')) + '\n'
Count_W = 'W' + '\t' + str(all_cogs_temp.count('W')) + '\n'
Count_Y = 'Y' + '\t' + str(all_cogs_temp.count('Y')) + '\n'
Count_Z = 'Z' + '\t' + str(all_cogs_temp.count('Z')) + '\n'

in_sum = Count_A + Count_B + Count_C + Count_D + Count_E + Count_F + Count_G \
         + Count_H + Count_I + Count_J + Count_K + Count_L + Count_M + Count_N \
         + Count_O + Count_P + Count_Q + Count_R + Count_S + Count_T + Count_U \
         + Count_V + Count_W + Count_Y + Count_Z

print(in_sum)
out.write(in_sum)
out.close()
