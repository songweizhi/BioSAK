library(phytools)

?phylosig

tree_file = ''

phylosig(tree_file, x, method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL,
         control=list())

data(mammal.tree)
data(mammal.data)


