
library(googleVis)

plot_99_200_length = read.csv('/Users/songweizhi/Desktop/googleVis_example_data.csv', header = TRUE)
Sankey_plot_99_200_length <- gvisSankey(plot_99_200_length,
                                        options = list(sankey = "{node:{colorMode:'unique', labelPadding: 10 },link:{colorMode:'source'}}", 
                                                       height = 1200, 
                                                       width = 600))
plot(Sankey_plot_99_200_length)
