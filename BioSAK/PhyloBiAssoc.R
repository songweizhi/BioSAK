suppressMessages(suppressWarnings(library("ape")))
suppressMessages(suppressWarnings(library("phytools")))
suppressMessages(suppressWarnings(library("optparse")))

################################################################################

option_list = list(
  make_option(c("-t", "--treefile"), type="character", default=NULL, help="tree file"),
  make_option(c("-d", "--datafile"), type="character", default=NULL, help="data file"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

tree_file = opt$treefile
data_file = opt$datafile

# Rscript PhyloBiAssoc.R -t demo.tre -d demo.txt
# phylosig        7.973475e-26    binaryPGLMM     0.03255813

################################################################################

geotree <- read.tree(tree_file)
geodata <- read.table(data_file, header = TRUE, sep = "\t")

# sort rows in df according to the order of tips in the tree
row.names(geodata) <- geodata[,1]
row.names(geodata) <- geodata$ID
geodata <- geodata[geotree$tip.label,]

for (i in colnames(geodata[, 3:ncol(geodata)])){
  
  # perform phylosig test
  phylosig_test <- phylosig(tree = geotree, x = setNames(geodata[, i], geodata$ID), method = "lambda", test = TRUE)
  phylosig_test_pvalue = phylosig_test$P
  
  # perform binaryPGLMM test if phylosig P-value <= 0.05 (indicating significant phylogenetic signal)
  # perform chi-squared test if phylosig P-value > 0.05  (indicating no phylogenetic signal)
  association_test = ''
  association_p_value = NA
  if (phylosig_test_pvalue <= 0.05) {
    binaryPGLMM_result <- binaryPGLMM(setNames(geodata[, i], geodata$ID) ~ geodata$cate, phy = geotree)
    association_test = 'binaryPGLMM'
    association_p_value = binaryPGLMM_result$B.pvalue[2]
  } else {
    chisq_test <- chisq.test(table(geodata$cate, setNames(geodata[, i], geodata$ID)))
    association_test = 'chisq.test'
    association_p_value = chisq_test$p.value
  }
  
  # print to screen
  cat(i, "phylosig", phylosig_test_pvalue, association_test, association_p_value, '\n', fill=FALSE, sep = "\t")
  
}
