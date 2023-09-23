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
# The header of the first two columns has to be "ID" and "cate".

################################################################################

geotree <- read.tree(tree_file)
geodata <- read.table(data_file, header = TRUE, sep = "\t")

# sort rows in df according to the order of tips in the tree
row.names(geodata) <- geodata[,1]
row.names(geodata) <- geodata$ID
geodata <- geodata[geotree$tip.label,]

cat('ID', "phylosig", "binaryPGLMM", "chisq.test", "coefficient", "significant", '\n', fill=FALSE, sep = "\t")
for (i in colnames(geodata[, 3:ncol(geodata)])){

  # perform phylosig test
  phylosig_test <- phylosig(tree = geotree, x = setNames(geodata[, i], geodata$ID), method = "lambda", test = TRUE)
  phylosig_test_pvalue = phylosig_test$P

  # perform binaryPGLMM test if phylosig P-value <= 0.05 (indicating significant phylogenetic signal)
  # perform chi-squared test if phylosig P-value > 0.05  (indicating no phylogenetic signal)
  # do nothing if phylosig returns NaN

  association_test = ''
  association_p_value = NA
  do_nothing = FALSE
  association_coefficient = 'na'
  significant = 'n'
  if (phylosig_test_pvalue == 'NaN') {
    do_nothing = TRUE
    significant = 'na'
  } else if (phylosig_test_pvalue <= 0.05) {
    binaryPGLMM_result <- binaryPGLMM(setNames(geodata[, i], geodata$ID) ~ geodata$cate, phy = geotree)
    association_test = 'binaryPGLMM'
    association_coefficient = binaryPGLMM_result$B[1]
    association_p_value = binaryPGLMM_result$B.pvalue[2]
    if (association_p_value <= 0.05) {
      significant = 'y'
    }

  } else {
    chisq_test <- chisq.test(table(geodata$cate, setNames(geodata[, i], geodata$ID)))
    association_test = 'chisq.test'
    association_p_value = chisq_test$p.value

    #cor_test <- cor.test(geodata$cate, geodata[, i])
    #association_coefficient = cor_test$estimate

    if (association_p_value <= 0.05) {
      significant = 'y'
    }
  }

  # print to screen
  if (do_nothing == FALSE) {

  if (association_test == 'binaryPGLMM'){
    cat(i, phylosig_test_pvalue, association_p_value, 'na', association_coefficient, significant, '\n', fill=FALSE, sep = "\t")
  }
  if (association_test == 'chisq.test'){
    cat(i, phylosig_test_pvalue, 'na', association_p_value, association_coefficient, significant, '\n', fill=FALSE, sep = "\t")
  }
  }
}
