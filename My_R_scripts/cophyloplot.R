library(ape)
library(randomcoloR)

tree1 <- rtree(40)
tree2 <- rtree(20)
association <- cbind(tree2$tip.label, tree2$tip.label)
palette<-distinctColorPalette(20)

print(association)
print(palette)

# cophyloplot(tree1, tree2, assoc=association, length.line=3, space=30, gap=3, col=palette, rotate=TRUE)
cophyloplot(tree1, tree2, assoc=association, length.line=3, space=30, gap=3, col=palette)

