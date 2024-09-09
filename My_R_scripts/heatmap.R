
# Use 'scale' to normalize
heatmap(data, scale="column")


data2 = read.table("/Users/songweizhi/Desktop/mtcars.txt", row.names=1, header = TRUE, sep = "\t")
data2_matrix = as.matrix(data2)

heatmap(data2_matrix)

