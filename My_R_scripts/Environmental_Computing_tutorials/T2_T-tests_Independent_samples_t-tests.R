# Independent samples t-tests
# http://environmentalcomputing.net/independent-samples-t-test/

River_pH_data = '/Users/songweizhi/PycharmProjects/BioSAK/My_R_scripts/Environmental_Computing_ tutorials/data/River_pH.csv'
River_pH <- read.csv(River_pH_data, header = T)
River_pH

t.test(pH ~ River_name, data = River_pH, var.equal = TRUE)


############################################################
#                 Assumptions to check                     #
############################################################

# Written 
# As a minimum, the observed t statistic, the P-value and the number of degrees of freedom 
# should be reported. For example, you could write “the pH was significantly higher in River 
# A than River B (independent samples t-test: t = 6.98, df = 18, P < 0.001)”.

############################################################
#                Communicating the results                 #
############################################################

# Visual
# Box plots or column graphs with error bars are effective ways of communicating the variation 
# in a single continuous response variable versus a single categorical predictor variable.

boxplot(pH~River_name, data = River_pH, xlab = "River name", ylab = "pH")



