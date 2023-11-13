
data_file = '/Users/songweizhi/Desktop/PhyloBiAssoc/demo.txt'
tree_file = '/Users/songweizhi/Desktop/PhyloBiAssoc/demo.tre'


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
    
    #binaryPGLMM_result <- binaryPGLMM(geodata[, i] ~ geodata$cate, phy = geotree)
    
    binaryPGLMM_result <- binaryPGLMM(setNames(geodata[, i], geodata$ID) ~ geodata$cate, phy = geotree)
    
    association_test = 'binaryPGLMM'
    association_p_value = binaryPGLMM_result$B.pvalue[2]
  } else {
    chisq_test <- chisq.test(table(geodata$cate, setNames(geodata[, i], geodata$ID)))
    association_test = 'chisq.test'
    association_p_value = chisq_test$p.value
  }
  
  # print to screen
  cat(i, "phylosig", phylosig_test_pvalue, association_test, association_p_value,fill=TRUE, sep = "\t")
  
}



