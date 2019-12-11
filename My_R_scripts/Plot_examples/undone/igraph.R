install.packages("igraph")
library(igraph)
el1 = read.csv('/Users/weizhisong/Desktop/edges_list.csv', header = FALSE)
el2 = as.matrix(el1)
el2
graph_from_edgelist(el, directed = TRUE)


el <- matrix( c("foo", "bar", "bar", "foobar"), nc = 2, byrow = TRUE)
el
graph_from_edgelist(el, directed = TRUE)



tkplot(el)
from_edgelist(el, directed = TRUE)

graph_from_edgelist(el)
make_ring(2)
ring(el)

... = graph_from_edgelist(el)


?"igraph"

graph_from_edgelist(cbind(1:10, c(2:10, 1)))


str(make_ring(10))
str(make_ring(10, directed = TRUE, mutual = TRUE))


g <- graph_from_literal( Alice-Bob-Cecil-Alice, Daniel-Cecil-Eugene, Cecil-Gordon )
g



