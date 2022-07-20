library(ape)
library(vegan)
library(optparse)

option_list = list(
  make_option(c("-a", "--tree1"),  type="character", default=NULL, help="the first tree"),
  make_option(c("-b", "--tree2"),  type="character", default=NULL, help="the second tree"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

TREE1 = read.tree(opt$tree1)
TREE2 = read.tree(opt$tree2)

D1 = cophenetic(TREE1)
D1 = D1[order(row.names(D1)),order(row.names(D1))]
D2 = cophenetic(TREE2)
D2 = D2[order(row.names(D2)),order(row.names(D2))]

mantel(xdis = D1, ydis = D2, permutations = 999)

