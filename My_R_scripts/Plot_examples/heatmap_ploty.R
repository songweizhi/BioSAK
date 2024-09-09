library(plotly)

plot_ly(z = volcano, type = "heatmap")

m <- matrix(rnorm(9), nrow = 3, ncol = 3)
plot_ly(z = m, x = c("a", "b", "c"), y = c("d", "e", "f"), type = "heatmap") # colorscale: "Hot", "Greys", "Greens"

# Custom colorscales
vals <- unique(scales::rescale(c(volcano)))
o <- order(vals, decreasing = FALSE)
cols <- scales::col_numeric("Blues", domain = NULL)(vals)
colz <- setNames(data.frame(vals[o], cols[o]), NULL)
plot_ly(z = volcano, colorscale = colz, type = "heatmap")
