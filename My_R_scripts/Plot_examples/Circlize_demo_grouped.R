library(circlize)

options(digits = 2)
mat1 = matrix(rnorm(25), nrow = 5)
rownames(mat1) = paste0("A", 1:5)
colnames(mat1) = paste0("B", 1:5)

mat2 = matrix(rnorm(25), nrow = 5)
rownames(mat2) = paste0("A", 1:5)
colnames(mat2) = paste0("C", 1:5)

mat3 = matrix(rnorm(25), nrow = 5)
rownames(mat3) = paste0("B", 1:5)
colnames(mat3) = paste0("C", 1:5)


mat = matrix(0, nrow = 10, ncol = 10)
rownames(mat) = c(rownames(mat2), rownames(mat3))
colnames(mat) = c(colnames(mat1), colnames(mat2))
mat[rownames(mat1), colnames(mat1)] = mat1
mat[rownames(mat2), colnames(mat2)] = mat2
mat[rownames(mat3), colnames(mat3)] = mat3

# Set larger gaps between groups. Here we manually adjust gap.after in circos.par().
circos.par(gap.after = rep(c(rep(1, 4), 8), 3))

chordDiagram(mat, annotationTrack = c("grid", "axis"),
             preAllocateTracks = list(
               track.height = uh(4, "mm"),
               track.margin = c(uh(4, "mm"), 0)
             ))

# highlight.sector(rownames(mat1), track.index = 1, col = "red",   text = "A", cex = 0.8, text.col = "white", niceFacing = TRUE)
# highlight.sector(colnames(mat1), track.index = 1, col = "green", text = "B", cex = 0.8, text.col = "white", niceFacing = TRUE)
# highlight.sector(colnames(mat2), track.index = 1, col = "blue",  text = "C", cex = 0.8, text.col = "white", niceFacing = TRUE)


print(rownames(mat1))
print(colnames(mat1))
print(colnames(mat2))

highlight.sector(c("A1", "A2", "A3", "A4", "A5"), track.index = 1, col = "lightblue", text = "A", cex = 0.8, text.col = "black", niceFacing = TRUE)
highlight.sector(c("B1", "B2", "B3", "B4", "B5"), track.index = 1, col = "lightblue", text = "B", cex = 0.8, text.col = "black", niceFacing = TRUE)
highlight.sector(c("C1", "C2", "C3", "C4", "C5"), track.index = 1, col = "lightblue", text = "C", cex = 0.8, text.col = "black", niceFacing = TRUE)

circos.clear()

rm(list=ls())

