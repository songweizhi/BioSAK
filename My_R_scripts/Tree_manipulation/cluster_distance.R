# Question: How To Do Clustering Based On Cutoff Distance Value In R
# https://www.biostars.org/p/52207/


# Create a sample data set.
y <- matrix(rnorm(50), 10, 5, dimnames=list(paste("g", 1:10, sep=""), paste("t", 1:5, sep="")))

# hierarchical cluster using Euclidean distance to get a range of distances
# for which the integer 2 is relevant
hr <- hclust(dist(y), method = "complete", members=NULL)

# examine the dendrogram
plot(hr)

# cut at a distance of 2, and get cluster memberships
myhcl <- cutree(hr, h=2)

# highlight our clusters on the dendrogram
rect.hclust(hr, h=2)

